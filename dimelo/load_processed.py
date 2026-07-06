import concurrent.futures
import gzip
import multiprocessing
import os
import warnings
from functools import partial
from multiprocessing import shared_memory
from pathlib import Path

import h5py
import numpy as np
import pysam
from tqdm.auto import tqdm

from . import test_data, utils

# the default chunk size is the number of bp to include per processing chunk in parallelization for loaders.
# 1e6 was empirically determined to be a good default: smaller than 1e5 we see slowdowns due to increased
# parallelization overhead, larger than 1e7 we see slowdowns due to worker utilization decreasing because even
# for whole chromosome processing there aren't always enough chunks to go around. In the 1e5-1e7 range, speed
# on 32 cores is fairly similar, but sitting in the middle of the range should support 10x more cores (beyond
# the reasonable upper bound) and 10x fewer cores (which is about the reasonable lower bound).
DEFAULT_CHUNK_SIZE = 1_000_000
AUTO_PARALLEL_MIN_TASKS = 64
AUTO_PARALLEL_TASKS_PER_CORE = 16
AUTO_PARALLEL_BATCHES_PER_WORKER = 8
AUTO_PARALLEL_CHUNKS_PER_TASK = 8
AUTO_PARALLEL_MEMORY_FRACTION = 0.5
AUTO_PARALLEL_PROCESS_OVERHEAD_BYTES = 64 * 1024 * 1024
_CGROUP_MEMORY_SENTINEL_LIMIT = 1 << 60
_TABIX_CACHE: dict[str, pysam.TabixFile] = {}

################################################################################################################
####                                           Loader wrappers                                              ####
################################################################################################################


def regions_to_list(
    function_handle,
    regions,
    window_size: int | None = None,
    quiet: bool = True,
    cores: int | None = None,
    split_large_regions: bool = False,
    **kwargs,
):
    """
    User-facing function.

    Run any standard load_processed pileup or extract loader loading each region from the region
    specifier into a new element of a list.

    Args:
        function_handle: the loader function you want to run.
        regions: the region specifier. Typically we expect to get many regions for this function, in the form of a list
            of strings or bed file paths. regions_to_list will run across all of these one-by-one returning a separate
            function return for each independent region.
        window_size: window around centers of regions, defaults to None
        quiet: disables progress bars
        cores: CPU cores across which to parallelize processing. Default to None, which means all available.
        split_large_regions: if True, regions will be run sequentially in parallelized chunks. If False,
            each individual region's chunks will be run sequentially but there will be parallelization across
            regions, i.e. each core will be assigned one region at a time by the executor. Set to True if you
            are running a small number of very large regions (e.g. one or two chromosomes), otherwise to to False (default).
        **kwargs: all necessary keyword arguments to pass down to the loader

    Returns:
        List(function_handle return objects per region)
    """
    executor = kwargs.pop("executor", None)
    regions_dict = utils.regions_dict_from_input(
        regions,
        window_size,
    )

    # Flatten regions into a list of (chromosome, start, end, strand)
    region_strings = [
        f"{chromosome}:{start}-{end},{strand}"
        for chromosome, region_list in regions_dict.items()
        for start, end, strand in region_list
    ]

    cores_to_run = _resolve_cores_for_task_count(
        requested_cores=cores,
        task_count=len(region_strings),
    )
    # quiet and cores logic below is driven by the following:
    # If the parallelization is within regions:
    #    (1) progress bars should happen within regions if at all, because we assume regions are
    #        large if they make sense to parallelize
    #    (2) the cores_to_run will be allocated to within-region parallelization, and the top-level
    #        jobs sequence is run sequentially
    # For single-core, avoid process-pool startup and run regions directly.
    if not split_large_regions and cores_to_run == 1:
        process_partial = partial(
            apply_loader_function_to_region,
            function_handle=function_handle,
            quiet=quiet,
            cores=1,
            **kwargs,
        )
        return list(
            tqdm(
                map(process_partial, region_strings),
                total=len(region_strings),
                desc="Loading data",
                disable=quiet,
                leave=False,
            )
        )

    process_partial = partial(
        _apply_loader_function_to_region_batch,
        function_handle=function_handle,
        quiet=quiet or not split_large_regions,
        cores=cores_to_run
        if split_large_regions
        else 1,  # if parallelization is within region
        **kwargs,
    )

    # Keep enough tasks per worker to amortize process overhead on large region lists.
    if split_large_regions:
        batch_size = 1
    else:
        target_batches = max(1, cores_to_run * AUTO_PARALLEL_BATCHES_PER_WORKER)
        batch_size = max(1, len(region_strings) // target_batches)
    region_batches = list(_iter_batches(region_strings, batch_size=batch_size))

    def _collect(mapped_batches) -> list:
        return [
            item
            for batch in tqdm(
                mapped_batches,
                total=len(region_batches),
                desc="Loading data",
                disable=quiet or split_large_regions,
                leave=False,
            )
            for item in batch
        ]

    if executor is not None and not split_large_regions:
        return _collect(executor.map(process_partial, region_batches))

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=1 if split_large_regions else cores_to_run
    ) as local_executor:
        return _collect(local_executor.map(process_partial, region_batches))


def apply_loader_function_to_region(region_string, function_handle, **kwargs):
    """
    apply_loader_function_to_region simply exists to convert position arguments into keyword arguments to make executor.map work

    Args:
        region_string: passed down with regions keyword
        function_handle: function to call with regions and other kwargs
        **kwargs: all keyword arguments passed to regions_to_list. These must be sufficient for whichever load_processed function
            if being referenced by function_handle
    Returns:
        function_handle return value
    """
    return function_handle(regions=region_string, **kwargs)


def _apply_loader_function_to_region_batch(region_batch, function_handle, **kwargs):
    return [
        function_handle(regions=region_string, **kwargs)
        for region_string in region_batch
    ]


def _iter_batches(items: list[str], *, batch_size: int):
    for idx in range(0, len(items), batch_size):
        yield items[idx : idx + batch_size]


def _parallel_chunk_batch_size(*, task_count: int, workers: int) -> int:
    target_tasks = max(1, workers * AUTO_PARALLEL_CHUNKS_PER_TASK)
    return max(1, task_count // target_tasks)


def _get_tabix_file(path: str | Path) -> pysam.TabixFile:
    key = str(path)
    cached = _TABIX_CACHE.get(key)
    if cached is None:
        cached = pysam.TabixFile(key)
        _TABIX_CACHE[key] = cached
    return cached


def _clear_tabix_cache() -> None:
    _TABIX_CACHE.clear()


def _resolve_cores_for_task_count(
    *, requested_cores: int | None, task_count: int
) -> int:
    """
    Resolve worker count for loader-style fanout.

    If cores are explicitly requested, honor that request. For automatic mode (cores=None),
    avoid expensive multiprocessing startup for small task batches.
    """
    resolved = utils.cores_to_run(requested_cores)
    if requested_cores is not None:
        return resolved

    cores_by_tasks = 1
    if task_count >= AUTO_PARALLEL_MIN_TASKS:
        cores_by_tasks = max(1, task_count // AUTO_PARALLEL_TASKS_PER_CORE)

    return max(1, min(resolved, cores_by_tasks))


def _parse_positive_int(value: str | None) -> int | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    digits = []
    for char in text:
        if char.isdigit():
            digits.append(char)
        else:
            break
    if not digits:
        return None
    parsed = int("".join(digits))
    return parsed if parsed > 0 else None


def _cgroup_memory_limit_bytes() -> int | None:
    candidates: list[int] = []
    for candidate_path in (
        "/sys/fs/cgroup/memory.max",
        "/sys/fs/cgroup/memory/memory.limit_in_bytes",
    ):
        try:
            value = Path(candidate_path).read_text(encoding="utf-8").strip()
        except Exception:
            continue
        if not value or value.lower() == "max":
            continue
        parsed = _parse_positive_int(value)
        if parsed is None or parsed >= _CGROUP_MEMORY_SENTINEL_LIMIT:
            continue
        candidates.append(parsed)
    if not candidates:
        return None
    return min(candidates)


def _slurm_memory_limit_bytes() -> int | None:
    mem_per_node_mb = _parse_positive_int(os.environ.get("SLURM_MEM_PER_NODE"))
    mem_per_cpu_mb = _parse_positive_int(os.environ.get("SLURM_MEM_PER_CPU"))
    cpus_per_task = _parse_positive_int(os.environ.get("SLURM_CPUS_PER_TASK"))
    cpus_on_node = _parse_positive_int(os.environ.get("SLURM_CPUS_ON_NODE"))

    candidates: list[int] = []
    if mem_per_node_mb is not None:
        candidates.append(mem_per_node_mb * 1024 * 1024)
    if mem_per_cpu_mb is not None:
        cpus = cpus_per_task or cpus_on_node
        if cpus is not None:
            candidates.append(mem_per_cpu_mb * cpus * 1024 * 1024)

    if not candidates:
        return None
    return min(candidates)


def _available_memory_bytes() -> int | None:
    override = _parse_positive_int(os.environ.get("DIMELO_AVAILABLE_MEMORY_BYTES"))
    if override is not None:
        return override

    candidates: list[int] = []
    try:
        import psutil  # type: ignore

        candidates.append(int(psutil.virtual_memory().available))
    except Exception:
        pass

    cgroup_limit = _cgroup_memory_limit_bytes()
    if cgroup_limit is not None:
        candidates.append(cgroup_limit)

    slurm_limit = _slurm_memory_limit_bytes()
    if slurm_limit is not None:
        candidates.append(slurm_limit)

    if not candidates:
        return None
    return min(candidates)


def _memory_limited_cores(
    *,
    requested_cores: int | None,
    task_count: int,
    bytes_per_task: int,
    extra_shared_bytes: int = 0,
) -> int:
    """
    Resolve worker count using CPU availability, task fanout, and available memory.

    Explicit cores override auto-tuning and are always returned unchanged.
    """
    baseline = _resolve_cores_for_task_count(
        requested_cores=requested_cores,
        task_count=task_count,
    )
    if requested_cores is not None:
        return baseline

    available = _available_memory_bytes()
    if available is None:
        return baseline

    budget = int(available * AUTO_PARALLEL_MEMORY_FRACTION) - max(
        0, int(extra_shared_bytes)
    )
    if budget <= 0:
        return 1

    per_worker = max(1, int(bytes_per_task)) + AUTO_PARALLEL_PROCESS_OVERHEAD_BYTES
    memory_cap = max(1, budget // per_worker)
    return max(1, min(baseline, memory_cap))


def _pileup_counts_from_bedmethyl_single_core(
    *,
    bedmethyl_file: str | Path,
    parsed_motif: utils.ParsedMotif,
    chunks_list: list[dict],
    single_strand: bool,
    quiet: bool,
) -> tuple[int, int]:
    """Single-process fallback used when only one core is requested."""
    source_tabix = _get_tabix_file(bedmethyl_file)
    total_valid = 0
    total_modified = 0

    for chunk in tqdm(
        chunks_list,
        total=len(chunks_list),
        disable=quiet,
        desc="Loading data",
        leave=False,
    ):
        chromosome = chunk["chromosome"]
        subregion_start = chunk["subregion_start"]
        subregion_end = chunk["subregion_end"]
        strand = chunk["strand"]

        if chromosome not in source_tabix.contigs:
            continue

        for row in source_tabix.fetch(
            chromosome, max(subregion_start, 0), subregion_end
        ):
            keep_basemod, _, modified_in_row, valid_in_row = process_pileup_row(
                row=row,
                parsed_motif=parsed_motif,
                region_strand=strand,
                single_strand=single_strand,
            )
            if keep_basemod:
                total_valid += valid_in_row
                total_modified += modified_in_row

    return total_modified, total_valid


def _pileup_vectors_from_bedmethyl_single_core(
    *,
    bedmethyl_file: str | Path,
    parsed_motif: utils.ParsedMotif,
    chunks_list: list[dict],
    region_len: int,
    single_strand: bool,
    regions_5to3prime: bool,
    quiet: bool,
) -> tuple[np.ndarray, np.ndarray]:
    """Single-process fallback used when only one core is requested."""
    source_tabix = _get_tabix_file(bedmethyl_file)
    valid_base_counts = np.zeros(region_len, dtype=np.int32)
    modified_base_counts = np.zeros(region_len, dtype=np.int32)

    for chunk in tqdm(
        chunks_list,
        total=len(chunks_list),
        disable=quiet,
        desc="Loading data",
        leave=False,
    ):
        chromosome = chunk["chromosome"]
        region_start = chunk["region_start"]
        region_end = chunk["region_end"]
        subregion_start = chunk["subregion_start"]
        subregion_end = chunk["subregion_end"]
        strand = chunk["strand"]

        flip_coords = regions_5to3prime and strand == "-"
        if flip_coords:
            subregion_offset = region_end - subregion_end
        else:
            subregion_offset = subregion_start - region_start

        if region_end - region_start > region_len:
            print(
                f"WARNING: You have specified a region at {chromosome}:{region_start}-{region_end} that is longer than the first region; the end of the region will be skipped. To make a profile plot with differently-sized region, consider using the window_size parameter to make a profile across centered windows."
            )

        subregion_len = subregion_end - subregion_start
        valid_base_subregion = np.zeros(subregion_len, dtype=np.int32)
        modified_base_subregion = np.zeros(subregion_len, dtype=np.int32)

        if chromosome in source_tabix.contigs:
            for row in source_tabix.fetch(
                chromosome, max(subregion_start, 0), subregion_end
            ):
                keep_basemod, genomic_coord, modified_in_row, valid_in_row = (
                    process_pileup_row(
                        row=row,
                        parsed_motif=parsed_motif,
                        region_strand=strand,
                        single_strand=single_strand,
                    )
                )
                if keep_basemod:
                    if flip_coords:
                        pileup_coord_in_subregion = subregion_end - genomic_coord - 1
                    else:
                        pileup_coord_in_subregion = genomic_coord - subregion_start
                    if pileup_coord_in_subregion < subregion_len:
                        valid_base_subregion[pileup_coord_in_subregion] += valid_in_row
                        modified_base_subregion[pileup_coord_in_subregion] += (
                            modified_in_row
                        )

        valid_base_counts[subregion_offset : subregion_offset + subregion_len] += (
            valid_base_subregion
        )
        modified_base_counts[subregion_offset : subregion_offset + subregion_len] += (
            modified_base_subregion
        )

    return modified_base_counts, valid_base_counts


################################################################################################################
####                                           Pileup loaders                                               ####
################################################################################################################


def pileup_counts_from_bedmethyl(
    bedmethyl_file: str | Path,
    motif: str,
    regions: str | Path | list[str | Path],
    window_size: int | None = None,
    single_strand: bool = False,
    quiet: bool = False,
    cores: int | None = None,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> tuple[int, int]:
    """
    User-facing function.

    Extract number of modified bases and total number of bases from the given bedmethyl file.
    Called by plotters or by the user.

    This function loops through all the provided regions and pulls those regions up in the input
    sorted and indexed bedmethyl file. For rows within those regions, checks that the motif
    is correct (i.e. sequence context, modified base, mod code, and optionally strand). All
    correct locations are included in the sum counts that get returned.

    If no regions are specified, returns the sum total for the motif of interest across the
    entire bedmethyl file.

    Args:
        bedmethyl_file: Path to bedmethyl file
        regions: Path to bed file specifying regions
        motif: type of modification to extract data for
        window_size: (currently disabled) window around center of region, +-window_size
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        quiet: disables progress bars
        cores: CPU cores across which to parallelize processing. Default to None, which means all available.
        chunk_size: size of genomic subregions to assign out to each process

    Returns:
        tuple containing counts of (modified_bases, total_bases)
    """

    parsed_motif = utils.ParsedMotif(motif)

    regions_dict = utils.regions_dict_from_input(regions, window_size)
    chunks_list = utils.process_chunks_from_regions_dict(
        regions_dict, chunk_size=chunk_size
    )

    cores_to_run = _memory_limited_cores(
        requested_cores=cores,
        task_count=len(chunks_list),
        bytes_per_task=0,
    )

    if cores_to_run == 1:
        return _pileup_counts_from_bedmethyl_single_core(
            bedmethyl_file=bedmethyl_file,
            parsed_motif=parsed_motif,
            chunks_list=chunks_list,
            single_strand=single_strand,
            quiet=quiet,
        )

    batch_size = _parallel_chunk_batch_size(
        task_count=len(chunks_list),
        workers=cores_to_run,
    )
    chunk_batches = list(_iter_batches(chunks_list, batch_size=batch_size))
    modified_base_count = 0
    valid_base_count = 0
    with concurrent.futures.ProcessPoolExecutor(max_workers=cores_to_run) as executor:
        futures = [
            executor.submit(
                _pileup_counts_process_chunk_batch_local,
                bedmethyl_file,
                parsed_motif,
                chunk_batch,
                single_strand,
            )
            for chunk_batch in chunk_batches
        ]
        for future in tqdm(
            concurrent.futures.as_completed(futures),
            total=len(futures),
            disable=quiet,
            desc="Loading data",
            leave=False,
        ):
            try:
                modified_in_batch, valid_in_batch = future.result()
                modified_base_count += int(modified_in_batch)
                valid_base_count += int(valid_in_batch)
            except Exception as err:
                raise RuntimeError("pileup_counts_process_chunk failed.") from err

    return modified_base_count, valid_base_count


def pileup_vectors_from_bedmethyl(
    bedmethyl_file: str | Path,
    motif: str,
    regions: str | Path | list[str | Path],
    window_size: int | None = None,
    single_strand: bool = False,
    regions_5to3prime: bool = False,
    quiet: bool = False,
    cores: int | None = None,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> tuple[np.ndarray, np.ndarray]:
    """
    User-facing function.

    Extract per-position pileup counts at valid motifs across one or more superimposed regions.
    Called by profile plotters, can also be used by a user directly.

    Returns two vectors:
    * Total number of times a modified base in the motif was found at each position
    * Total number of times the motif was found at each position

    This function loops through all the provided regions and fetches those regions from the
    bedmethyl file. For rows within those regions, it checks that the motif
    is correct (i.e. sequence context, modified base, mod code, and optionally strand). It then adds
    to two vectors (mod and valid). By default all regions are assumed to
    be the same size (the size of the first region).

    If regions_5to3prime is set to True, then negative strand regions are flipped to that all regions
    are superimposed along the 5 prime to 3 prime direction, which can be helpful if there is
    directionality to the signal (e.g. upstream v downstream relative to TSSs, TF binding sites, and so on).
    A region must be provided because otherwise there is no way to know what vector to return.
    However, a region can be a whole chromosome if desired.

    Args:
        bedmethyl_file: Path to bedmethyl file
        regions: Path to bed file specifying centered equal-length regions
        motif: type of modification to extract data for
        window_size: the extent in either direction for windows around the center of regions.
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        regions_5to3prime: True means negative strand regions get flipped, False means no flipping
        quiet: disables progress bars
        cores: CPU cores across which to parallelize processing. Default to None, which means all available.
        chunk_size: size of genomic subregions to assign out to each process

    Returns:
        tuple containing (modified_base_counts, valid_base_counts)
    """

    parsed_motif = utils.ParsedMotif(motif)

    regions_dict = utils.regions_dict_from_input(regions, window_size)
    chunks_list = utils.process_chunks_from_regions_dict(
        regions_dict, chunk_size=chunk_size
    )
    # Peek at a region to figure out what size the vectors should be
    first_key = next(iter(regions_dict))
    first_tuple = regions_dict[first_key][0]
    region_len = first_tuple[1] - first_tuple[0]

    max_chunk_span = max(
        (
            max(0, int(chunk["subregion_end"]) - int(chunk["subregion_start"]))
            for chunk in chunks_list
        ),
        default=0,
    )
    per_worker_bytes = 2 * max_chunk_span * np.dtype(np.int32).itemsize
    shared_bytes = 2 * region_len * np.dtype(np.int32).itemsize
    cores_to_run = _memory_limited_cores(
        requested_cores=cores,
        task_count=len(chunks_list),
        bytes_per_task=per_worker_bytes,
        extra_shared_bytes=shared_bytes,
    )

    # Avoid process-pool/shared-memory startup overhead for single-core workloads.
    if cores_to_run == 1:
        return _pileup_vectors_from_bedmethyl_single_core(
            bedmethyl_file=bedmethyl_file,
            parsed_motif=parsed_motif,
            chunks_list=chunks_list,
            region_len=region_len,
            single_strand=single_strand,
            regions_5to3prime=regions_5to3prime,
            quiet=quiet,
        )

    # Initialize shared memory as numpy arrays to make it easy to map to buffer in subprocesses
    shm_valid = shared_memory.SharedMemory(
        create=True, size=(region_len) * np.dtype(np.int32).itemsize
    )
    shm_modified = shared_memory.SharedMemory(
        create=True, size=(region_len) * np.dtype(np.int32).itemsize
    )

    manager = multiprocessing.Manager()
    lock = manager.Lock()

    batch_size = _parallel_chunk_batch_size(
        task_count=len(chunks_list),
        workers=cores_to_run,
    )
    chunk_batches = list(_iter_batches(chunks_list, batch_size=batch_size))
    with concurrent.futures.ProcessPoolExecutor(max_workers=cores_to_run) as executor:
        futures = [
            executor.submit(
                _pileup_vectors_process_chunk_batch,
                bedmethyl_file,
                parsed_motif,
                chunk_batch,
                region_len,
                shm_modified.name,
                shm_valid.name,
                lock,
                single_strand,
                regions_5to3prime,
            )
            for chunk_batch in chunk_batches
        ]
        for future in tqdm(
            concurrent.futures.as_completed(futures),
            total=len(futures),
            disable=quiet,
            desc="Loading data",
            leave=False,
        ):
            try:
                future.result()
            except Exception as err:
                raise RuntimeError("pileup_vectors_process_chunk failed.") from err

    # We need to convert these shared memory buffers to numpy arrays which
    # we then copy, so that they no longer reference the shared memory which
    # will soon be de-allocated
    modified_base_counts = np.copy(
        np.ndarray((region_len,), dtype=np.int32, buffer=shm_modified.buf)
    )
    valid_base_counts = np.copy(
        np.ndarray((region_len,), dtype=np.int32, buffer=shm_valid.buf)
    )
    # Close and unlink shared memory - not fully handled by garbage collection otherwise
    shm_modified.close()
    shm_modified.unlink()
    shm_valid.close()
    shm_valid.unlink()

    return modified_base_counts, valid_base_counts


def counts_from_fake(*args, **kwargs) -> tuple[int, int]:
    """
    Test helper function.

    Generates a fake set of enrichment counts. Ignores all arguments.

    Returns:
        tuple containing counts of (modified_bases, total_bases)
    """
    window_halfsize = 500
    return test_data.fake_peak_enrichment(halfsize=window_halfsize, peak_height=0.15)


def vector_from_fake(window_size: int, *args, **kwargs) -> np.ndarray:
    """
    Test helper function.

    Generates a fake peak trace. Ignores all arguments except window_size.

    Args:
        window_size: halfsize of the window; how far the window stretches on either side of the center point

    Returns:
        vector of fraction modified bases calculated for each position; float values between 0 and 1
    """
    return test_data.fake_peak_enrichment_profile(
        halfsize=window_size, peak_height=0.15
    )


def pileup_vectors_process_chunk(
    bedmethyl_file,
    parsed_motif,
    chunk,
    region_len,
    shm_name_modified,
    shm_name_valid,
    lock,
    single_strand,
    regions_5to3prime,
) -> None:
    _pileup_vectors_process_chunk_batch(
        bedmethyl_file=bedmethyl_file,
        parsed_motif=parsed_motif,
        chunk_batch=[chunk],
        region_len=region_len,
        shm_name_modified=shm_name_modified,
        shm_name_valid=shm_name_valid,
        lock=lock,
        single_strand=single_strand,
        regions_5to3prime=regions_5to3prime,
    )


def _pileup_vectors_local_contribution(
    *,
    source_tabix: pysam.TabixFile,
    parsed_motif: utils.ParsedMotif,
    chunk: dict,
    region_len: int,
    single_strand: bool,
    regions_5to3prime: bool,
) -> tuple[int, np.ndarray, np.ndarray] | None:
    chromosome = chunk["chromosome"]
    region_start = chunk["region_start"]
    region_end = chunk["region_end"]
    subregion_start = chunk["subregion_start"]
    subregion_end = chunk["subregion_end"]
    strand = chunk["strand"]

    flip_coords = regions_5to3prime and strand == "-"
    if flip_coords:
        subregion_offset = region_end - subregion_end
    else:
        subregion_offset = subregion_start - region_start

    if region_end - region_start > region_len:
        print(
            f"WARNING: You have specified a region at {chromosome}:{region_start}-{region_end} that is longer than the first region; the end of the region will be skipped. To make a profile plot with differently-sized region, consider using the window_size parameter to make a profile across centered windows."
        )

    subregion_span = subregion_end - subregion_start
    valid_base_subregion = np.zeros(subregion_span, dtype=np.int32)
    modified_base_subregion = np.zeros(subregion_span, dtype=np.int32)

    if chromosome not in source_tabix.contigs:
        return subregion_offset, valid_base_subregion, modified_base_subregion

    for row in source_tabix.fetch(chromosome, max(subregion_start, 0), subregion_end):
        keep_basemod, genomic_coord, modified_in_row, valid_in_row = process_pileup_row(
            row=row,
            parsed_motif=parsed_motif,
            region_strand=strand,
            single_strand=single_strand,
        )
        if not keep_basemod:
            continue
        if flip_coords:
            pileup_coord_in_subregion = subregion_end - genomic_coord - 1
        else:
            pileup_coord_in_subregion = genomic_coord - subregion_start
        if pileup_coord_in_subregion < subregion_span:
            valid_base_subregion[pileup_coord_in_subregion] += valid_in_row
            modified_base_subregion[pileup_coord_in_subregion] += modified_in_row

    return subregion_offset, valid_base_subregion, modified_base_subregion


def _pileup_vectors_process_chunk_batch(
    bedmethyl_file,
    parsed_motif,
    chunk_batch,
    region_len,
    shm_name_modified,
    shm_name_valid,
    lock,
    single_strand,
    regions_5to3prime,
) -> None:
    """
    Worker helper that processes a batch of chunks and acquires the shared-memory lock once.
    """
    source_tabix = _get_tabix_file(bedmethyl_file)
    existing_valid = shared_memory.SharedMemory(name=shm_name_valid)
    existing_modified = shared_memory.SharedMemory(name=shm_name_modified)
    valid_base_counts = np.ndarray(
        (region_len,), dtype=np.int32, buffer=existing_valid.buf
    )
    modified_base_counts = np.ndarray(
        (region_len,), dtype=np.int32, buffer=existing_modified.buf
    )

    updates: list[tuple[int, np.ndarray, np.ndarray]] = []
    for chunk in chunk_batch:
        contribution = _pileup_vectors_local_contribution(
            source_tabix=source_tabix,
            parsed_motif=parsed_motif,
            chunk=chunk,
            region_len=region_len,
            single_strand=single_strand,
            regions_5to3prime=regions_5to3prime,
        )
        if contribution is not None:
            updates.append(contribution)

    if updates:
        with lock:
            for (
                subregion_offset,
                valid_base_subregion,
                modified_base_subregion,
            ) in updates:
                subregion_span = len(valid_base_subregion)
                valid_base_counts[
                    subregion_offset : subregion_offset + subregion_span
                ] += valid_base_subregion
                modified_base_counts[
                    subregion_offset : subregion_offset + subregion_span
                ] += modified_base_subregion

    existing_modified.close()
    existing_valid.close()


def pileup_counts_process_chunk(
    bedmethyl_file,
    parsed_motif,
    chunk,
    shm_name_modified,
    shm_name_valid,
    lock,
    single_strand,
) -> None:
    modified_subregion_counts, valid_subregion_counts = (
        _pileup_counts_process_chunk_batch_local(
            bedmethyl_file=bedmethyl_file,
            parsed_motif=parsed_motif,
            chunk_batch=[chunk],
            single_strand=single_strand,
        )
    )
    existing_valid = shared_memory.SharedMemory(name=shm_name_valid)
    existing_modified = shared_memory.SharedMemory(name=shm_name_modified)
    valid_base_counts = np.ndarray((1,), dtype=np.int32, buffer=existing_valid.buf)
    modified_base_counts = np.ndarray(
        (1,), dtype=np.int32, buffer=existing_modified.buf
    )
    with lock:
        valid_base_counts[0] += valid_subregion_counts
        modified_base_counts[0] += modified_subregion_counts
    existing_valid.close()
    existing_modified.close()


def _pileup_counts_process_chunk_batch_local(
    bedmethyl_file,
    parsed_motif,
    chunk_batch,
    single_strand,
) -> tuple[int, int]:
    """
    Worker helper for counts mode that accumulates locally and returns one tuple per batch.
    """
    source_tabix = _get_tabix_file(bedmethyl_file)
    valid_total = 0
    modified_total = 0

    for chunk in chunk_batch:
        chromosome = chunk["chromosome"]
        subregion_start = chunk["subregion_start"]
        subregion_end = chunk["subregion_end"]
        strand = chunk["strand"]

        if chromosome not in source_tabix.contigs:
            continue
        for row in source_tabix.fetch(
            chromosome, max(subregion_start, 0), subregion_end
        ):
            keep_basemod, _, modified_in_row, valid_in_row = process_pileup_row(
                row=row,
                parsed_motif=parsed_motif,
                region_strand=strand,
                single_strand=single_strand,
            )
            if keep_basemod:
                valid_total += valid_in_row
                modified_total += modified_in_row

    return modified_total, valid_total


def process_pileup_row(
    row: str,
    parsed_motif: utils.ParsedMotif,
    region_strand: str,
    single_strand: bool = False,
) -> tuple[bool, int, int, int]:
    """
    Helper function designed for pileup_counts_from_bedmethyl via pileup_counts_process_chunk, pileup_vectors_from_bedmethyl
    via pileup_vectors_process_chunk, and export.pileup_to_bigwig; changes to logic here may impact some or all of
    these.

    Process a row from a pileup, determining whether the basemod is relevant and passing back its coordinate,
    modification count, and valid read count.

    Args:
        row: a string row from a bedmethyl file
        parsed_motif: a ParsedMotif object
        region_strand: the strand from the query region
        single_strand: True if only mods on the region_strand are to be kept

    Returns: keep_basemod, genomic_coord, modified_in_row, valid_in_row. Values are provided even if keep_basemod is False.
    """
    tabix_fields = row.split("\t")
    pileup_basemod = tabix_fields[3]

    if single_strand and tabix_fields[5] != region_strand:
        # We are on the wrong strand, can't keep this position.
        return (False, 0, 0, 0)

    pileup_basemod_parts = pileup_basemod.split(",")
    if len(pileup_basemod_parts) == 3:
        pileup_modname, pileup_motif, pileup_mod_coord = pileup_basemod_parts
        keep_basemod = (
            pileup_motif == parsed_motif.motif_seq
            and int(pileup_mod_coord) == parsed_motif.modified_pos
            and pileup_modname in parsed_motif.mod_codes
        )
    elif len(pileup_basemod_parts) == 1:
        keep_basemod = pileup_basemod in parsed_motif.mod_codes
    else:
        raise ValueError(
            f"Unexpected format in bedmethyl file: {row} contains {pileup_basemod} which cannot be parsed."
        )

    if not keep_basemod:
        return (False, 0, 0, 0)

    genomic_coord = int(tabix_fields[1])

    # Modern modkit pileup bedMethyl format (0.6.x) emits full 18-column rows:
    # chrom, start, end, mod, score, strand, thickStart, thickEnd, color,
    # valid_cov, frac_mod, N_mod, N_canonical, N_other_mod, N_delete, N_fail,
    # N_diff, N_nocall.
    # Legacy/older formats may encode summary info in tabix_fields[9] as a
    # whitespace-delimited payload where [0]=valid and [2]=modified.
    try:
        if len(tabix_fields) >= 13:
            valid_in_row = int(float(tabix_fields[9]))
            modified_in_row = int(float(tabix_fields[11]))
        else:
            pileup_info = tabix_fields[9].split()
            valid_in_row = int(float(pileup_info[0]))
            modified_in_row = int(float(pileup_info[2]))
    except Exception:
        if not getattr(process_pileup_row, "_warned_malformed_bedmethyl_row", False):
            warnings.warn(
                "Skipping malformed bedMethyl row while loading pileup data. "
                "This warning is shown once per process.",
                RuntimeWarning,
                stacklevel=2,
            )
            process_pileup_row._warned_malformed_bedmethyl_row = True
        return (False, 0, 0, 0)

    return (True, genomic_coord, modified_in_row, valid_in_row)


################################################################################################################
####                                        Single read loaders                                             ####
################################################################################################################


def _validate_subset_parameters(subset_parameters: dict | None) -> None:
    if subset_parameters is None:
        return
    if not isinstance(subset_parameters, dict):
        raise ValueError("subset_parameters must be provided as a dictionary.")
    if "array" in subset_parameters:
        raise ValueError(
            "subset_parameters cannot include 'array'; this is set internally."
        )
    if "n" not in subset_parameters and "frac" not in subset_parameters:
        raise ValueError("subset_parameters must include at least one of n or frac.")


def _subset_indices(
    relevant_read_indices: np.ndarray,
    subset_parameters: dict | None,
    seed: int | None = None,
) -> np.ndarray:
    if subset_parameters is None or relevant_read_indices.size == 0:
        return relevant_read_indices
    # Reproducibility fix: thread an optional seed into random_sample so subsetting
    # is deterministic when a seed is supplied (None preserves prior behavior).
    return np.sort(
        utils.random_sample(relevant_read_indices, seed=seed, **subset_parameters)
    )


def read_vectors_from_hdf5(
    file: str | Path,
    motifs: list[str],
    regions: str | Path | list[str | Path] | None = None,
    window_size: int | None = None,
    single_strand: bool = False,
    sort_by: str | list[str] = ["chromosome", "region_start", "read_start"],
    calculate_mod_fractions: bool = True,
    quiet: bool = True,  # currently unused; change to default False when pbars are implemented
    cores: int | None = None,  # currently unused
    subset_parameters: dict | None = None,
    span_full_window: bool = False,
    random_sample_n_reads: int | None = None,
    min_read_length_bp: int | None = None,
    random_state: int | None = None,
) -> tuple[list[tuple], list[str], dict | None]:
    """
    User-facing function.

    Pulls a list of read data out of an .h5 file containing processed read vectors, formatted
    for read-by-read vector processing downstream use cases.

    The flow of operation here is we load up the h5 file then loop through our regions and pick
    out reads corresponding to our criteria. Criteria include chromosome, read starts and ends
    (compared to region starts and ends), motif, and strand (if single_strand is True). The indices
    for the desired reads are identified region-by-region, then all the reads for the region (or
    the whole h5, if no region is passed) are loaded using the process_data function and put into
    a list. The bytes are then decoded for the array entries, which are manually compressed because
    h5py wasn't behaving.

    There's some adjustment for the raw probability (no thresh) to match modkit extract outputs.
    Specifically, the 0-255 8bit int has 0.5 added before dividing by 256, giving mod qualities
    between 0.001953 and 0.99805 for bases in valid motifs. (Invalid positions have zeros.)

    After this processing, we calculate modification fractions, sort, and return.

    TODO: Implement progress bars and parallelization as with pileup loaders

    Args:
        file: Path to an hdf5 (.h5) file containing modification data for single reads,
            stored in datasets read_name, chromosome, read_start,
            read_end, base modification motif, mod_vector, and val_vector.
        regions: Single or list of Path objects or strings. Path objects must point to .bed
            files, strings can be .bed paths or region string in the format chrX:XXX-XXX.
            All should all be regions for which your original .bam file had reads extracted,
            although by design this method will not raise an error if any region contains
            zero reads, as this may simply be a matter of low read depth.
            If no regions are specified, the entire .h5 file will be returned. This may cause
            memory issues.
        motifs: types of modification to extract data for. Motifs are specified as
            {DNA_sequence},{position_of_modification}. For example, a methylated adenine is specified
            as 'A,0' and CpG methylation is specified as 'CG,0'.
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        window_size: An optional parameter for creating centered windows for the provided regions.
            If provided, all regions will be adjusted to be the same size and centered. If not provided,
            all regions should already be the same size, or there should be only one.
        sort_by: Read properties by which to sort, either one string or a list of strings. Options
            include chromosome, region_start, region_end, read_start, read_end, and motif. More to
            be added in future.
        quiet: silences progress bars (currently unused)
        cores: cores across which to parallelize processes (currently unused)
        subset_parameters: Parameters to pass to the utils.random_sample() method, to subset the
            reads to be returned. If not None, at least one of n or frac must be provided. The array
            parameter should not be provided here.
        span_full_window: If True, only load reads that fully span the window defined by region_start-region_end
        random_sample_n_reads: Optional global random sample size applied after filtering and
            region/motif selection. Samples unique read names (not tuple rows), then retains all
            motif rows associated with those reads.
        min_read_length_bp: Optional minimum read length filter (bp) applied before loading vectors.
        random_state: Optional integer seed threaded into utils.random_sample for all subsampling
            (both subset_parameters per-region draws and the random_sample_n_reads global draw).
            When None (default), sampling uses the module-level unseeded rng (prior behavior).

    Returns:
        a list of tuples, each tuple containing all datasets corresponding to an individual read that
        was within the specified regions.
        a list of strings, naming the datasets returned.
        a regions_dict, containing lists of (region_start,region_end) coordinates by chromosome/contig.

    TODO: The way the subsetting is implemented is confusing, in that you need to pass all but one of
        the available parameters.
    """
    _validate_subset_parameters(subset_parameters)
    if random_sample_n_reads is not None and int(random_sample_n_reads) <= 0:
        raise ValueError("random_sample_n_reads must be > 0 when provided.")
    if min_read_length_bp is not None and int(min_read_length_bp) < 0:
        raise ValueError("min_read_length_bp must be >= 0 when provided.")
    random_sample_n_reads = (
        int(random_sample_n_reads) if random_sample_n_reads is not None else None
    )
    min_read_length_bp = (
        int(min_read_length_bp) if min_read_length_bp is not None else None
    )

    with h5py.File(file, "r") as h5:
        datasets: list[str] = [
            name for name, obj in h5.items() if isinstance(obj, h5py.Dataset)
        ]
        if "threshold" in h5:
            # we are looking at an .h5 file with the new, much better compressed format that does
            # not know the data type intrinsically for mod and val vectors, so we must check
            metadata_datasets = {"threshold", "threshold_by_motif_json"}
            readwise_datasets = [
                dataset for dataset in datasets if dataset not in metadata_datasets
            ]
            compressed_binary_datasets = ["mod_vector", "val_vector"]
            threshold_applied_to_h5 = h5["threshold"][()]
            binarized = not np.isnan(threshold_applied_to_h5)
        else:
            # backwards compatible with the old h5 file structure
            # If we remove backwards compatibility, beta test (Feb 2024) h5 extractions will not run
            readwise_datasets = datasets
            compressed_binary_datasets = []
            binarized = True  # in this case all this will do is make it so we don't apply a +1/512 correction to the mod_vector

        # Pre-load metadata so we can identify reads to pull from file
        read_chromosomes = np.array(h5["chromosome"], dtype=str)
        read_starts = np.array(h5["read_start"])
        read_ends = np.array(h5["read_end"])
        read_lengths = read_ends - read_starts
        read_motifs = np.array(h5["motif"], dtype=str)
        ref_strands = np.array(h5["strand"], dtype=str)

        # Identify reads to load, then load them
        if regions is not None:
            regions_dict = utils.regions_dict_from_input(
                regions=regions,
                window_size=window_size,
            )
            read_tuples_raw = []
            for chrom, region_list in regions_dict.items():
                for region_start, region_end, region_strand in region_list:
                    # Find the read indices that we want to load
                    # TODO: consider building this up and then loading all at the end, chunked
                    # TODO: consolidate logic into clear variables
                    relevant_read_indices = np.flatnonzero(
                        (read_ends > region_start)
                        & (read_starts < region_end)
                        & (read_starts <= region_start if span_full_window else True)
                        & (read_ends >= region_end if span_full_window else True)
                        & (
                            read_lengths >= min_read_length_bp
                            if min_read_length_bp is not None
                            else True
                        )
                        & np.isin(read_motifs, motifs)
                        & (read_chromosomes == chrom)
                        & (
                            (not single_strand)
                            | (region_strand not in ["+", "-"])
                            | (ref_strands == region_strand)
                        )
                    )
                    relevant_read_indices = _subset_indices(
                        relevant_read_indices,
                        subset_parameters=subset_parameters,
                        seed=random_state,  # deterministic subset when random_state given
                    )
                    read_tuples_raw += list(
                        zip(
                            *(
                                retrieve_h5_data(
                                    h5=h5,
                                    dataset=dataset,
                                    indices=relevant_read_indices,
                                    compressed=dataset in compressed_binary_datasets,
                                    dtype=np.uint8,
                                    decompressor=gzip.decompress,
                                    binarized=binarized,
                                )
                                for dataset in readwise_datasets
                            ),
                            [region_start for _ in relevant_read_indices],
                            [region_end for _ in relevant_read_indices],
                            [region_strand for _ in relevant_read_indices],
                            strict=False,
                        )
                    )
        else:
            regions_dict = None
            relevant_read_indices = np.flatnonzero(
                np.isin(read_motifs, motifs)
                & (
                    read_lengths >= min_read_length_bp
                    if min_read_length_bp is not None
                    else True
                )
            )
            relevant_read_indices = _subset_indices(
                relevant_read_indices,
                subset_parameters=subset_parameters,
                seed=random_state,  # deterministic subset when random_state given
            )
            read_tuples_raw = list(
                zip(
                    *(
                        retrieve_h5_data(
                            h5=h5,
                            dataset=dataset,
                            indices=relevant_read_indices,
                            compressed=dataset in compressed_binary_datasets,
                            dtype=np.uint8,
                            decompressor=gzip.decompress,
                            binarized=binarized,
                        )
                        for dataset in readwise_datasets
                    ),
                    [-1 for _ in relevant_read_indices],
                    [-1 for _ in relevant_read_indices],
                    ["." for _ in relevant_read_indices],
                    strict=False,
                )
            )
    #  We add region information (start, end, and strand; chromosome is already present!)
    # so that it is possible to sort by and process based on these
    readwise_datasets += ["region_start", "region_end", "region_strand", "read_length"]

    # This is sanitizing the dataset entries and adjusting prob values if needed
    if binarized:
        read_tuples_processed = [
            convert_bytes_to_strings(tup) for tup in read_tuples_raw
        ]
    else:
        read_tuples_processed = [
            adjust_mod_probs_in_tuples(
                convert_bytes_to_strings(tup),
                readwise_datasets.index("mod_vector"),
                readwise_datasets.index("val_vector"),
            )
            for tup in read_tuples_raw
        ]

    # Optional global random sample by unique read name, applied after filtering/loading.
    if random_sample_n_reads is not None and read_tuples_processed:
        read_name_idx_pre = readwise_datasets.index("read_name")
        unique_read_names = np.array(
            sorted({row[read_name_idx_pre] for row in read_tuples_processed}, key=str)
        )
        if random_sample_n_reads < unique_read_names.size:
            sampled_names = set(
                # deterministic global read-name sample when random_state given
                utils.random_sample(
                    unique_read_names, n=random_sample_n_reads, seed=random_state
                ).tolist()
            )
            read_tuples_processed = [
                row
                for row in read_tuples_processed
                if row[read_name_idx_pre] in sampled_names
            ]

    read_start_idx = readwise_datasets.index("read_start")
    read_end_idx = readwise_datasets.index("read_end")
    read_name_idx = readwise_datasets.index("read_name")
    motif_idx = readwise_datasets.index("motif")
    mod_vector_idx = readwise_datasets.index("mod_vector")
    val_vector_idx = readwise_datasets.index("val_vector")

    if calculate_mod_fractions:
        # Add the MOTIF_mod_fraction entries to the readwise_datasets list for future reference in sorting
        readwise_datasets += [f"{motif}_mod_fraction" for motif in motifs]
        # dict[read_name][motif_index]=modified fraction of motif in read, float
        motif_idx_by_name = {motif: i for i, motif in enumerate(motifs)}
        mod_fractions_by_read_name: dict[str, list[float]] = {}
        for read_tuple in read_tuples_processed:
            motif = read_tuple[motif_idx]
            motif_fraction_idx = motif_idx_by_name.get(motif)
            if motif_fraction_idx is None:
                continue
            mod_sum = np.sum(read_tuple[mod_vector_idx])
            val_sum = np.sum(read_tuple[val_vector_idx])
            mod_fraction = mod_sum / val_sum if val_sum > 0 else 0
            read_name = read_tuple[read_name_idx]
            if read_name not in mod_fractions_by_read_name:
                mod_fractions_by_read_name[read_name] = [0.0] * len(motifs)
            mod_fractions_by_read_name[read_name][motif_fraction_idx] = mod_fraction

        empty_mod_fractions = tuple(0.0 for _ in motifs)
        mod_fractions_by_read_name_tup = {
            read_name: tuple(mod_fractions)
            for read_name, mod_fractions in mod_fractions_by_read_name.items()
        }

        read_tuples_all = []
        for read_tuple in read_tuples_processed:
            read_length = read_tuple[read_end_idx] - read_tuple[read_start_idx]
            read_tuples_all.append(
                tuple(val for val in read_tuple)
                + (read_length,)
                + mod_fractions_by_read_name_tup.get(
                    read_tuple[read_name_idx], empty_mod_fractions
                )
            )
    else:
        read_tuples_all = []
        for read_tuple in read_tuples_processed:
            read_length = read_tuple[read_end_idx] - read_tuple[read_start_idx]
            read_tuples_all.append(tuple(val for val in read_tuple) + (read_length,))

    ## Sort the reads

    # Normalize sort_by to a list of tuples (field, order)
    # Support formats:
    #   - "field" or ["field1", "field2"] -> [("field1", "asc"), ("field2", "asc")]
    #   - [("field1", "desc"), "field2"] -> [("field1", "desc"), ("field2", "asc")]

    # Enforce that sort_by is a list
    if not isinstance(sort_by, list):
        sort_by = [sort_by]

    # Parse into (field, order) tuples
    sort_by_normalized = []
    for item in sort_by:
        if isinstance(item, tuple):
            field, order = item
            if order not in ["asc", "desc"]:
                raise ValueError(
                    f"Sort order must be 'asc' or 'desc', got '{order}' for field '{field}'"
                )
            sort_by_normalized.append((field, order))
        else:
            # Default to ascending order
            sort_by_normalized.append((item, "asc"))

    # If 'shuffle' appears anywhere in sort_by, we first shuffle the list
    if any(field == "shuffle" for field, _ in sort_by_normalized):
        utils.rng.shuffle(read_tuples_all)

    # Build sorting configuration
    try:
        sort_config = [
            (readwise_datasets.index(field), order == "desc")
            for field, order in sort_by_normalized
            if field != "shuffle"
        ]
    except ValueError as e:
        raise ValueError(
            f"Sorting error. {e}. Datasets include {readwise_datasets}. If you need mod fraction sorting make sure you are not setting calculate_read_fraction to False."
        ) from e

    if len(sort_config) > 0:
        # Use stable sort from right to left (reverse order of sort keys)
        sorted_read_tuples = read_tuples_all
        for idx, is_reverse in reversed(sort_config):
            sorted_read_tuples = sorted(
                sorted_read_tuples, key=lambda x: x[idx], reverse=is_reverse
            )
    else:
        sorted_read_tuples = read_tuples_all

    return sorted_read_tuples, readwise_datasets, regions_dict


def readwise_binary_modification_arrays(
    file: str | Path,
    motifs: list[str],
    regions: str | Path | list[str | Path],
    window_size: int | None = None,
    regions_5to3prime: bool = False,
    single_strand: bool = False,
    sort_by: str | list[str] = ["chromosome", "region_start", "read_start"],
    thresh: float | None = None,
    relative: bool = True,
    quiet: bool = True,  # currently unused; change to default False when pbars are implemented
    cores: int | None = None,  # currently unused
    subset_parameters: dict | None = None,
    random_sample_n_reads: int | None = None,
    min_read_length_bp: int | None = None,
) -> tuple[list[np.ndarray], np.ndarray[int], np.ndarray[str], dict | None]:
    """
    Primarily designed as a helper function for single-read plotting, but can be used by a user.

    Pulls a list of read data out of a file containing processed read vectors, formatted with
    seaborn plotting in mind. Currently we only support .h5 files.

    After running read_vectors_from_hdf5, this function takes the baton to convert the names of
    the sorted reads into integer indices, then goes through the reads and strips down the mod
    vectors to simply a list of modified positions (applying a threshold if one has not already
    been applied). Mod positions are by default expressed relative to the region from which
    the read was identified, allowing for nice plotting, but can also be expressed in absolute
    coordinates. If positions are relative, regions_5to3prime can be used to show all regions
    as upstream-to-downstream along their respective strands.

    TODO: Implement progress bars and parallelization as with pileup loaders

    Args:
        file: Path to an hdf5 (.h5) file containing modification data for single reads,
            stored in datasets read_name, chromosome, read_start,
            read_end, base modification motif, mod_vector, and val_vector.
        regions: Single or list of Path objects or strings. Path objects must point to .bed
            files, strings can be .bed paths or region string in the format chrX:XXX-XXX.
            All should all be regions for which your original .bam file had reads extracted,
            although by design this method will not raise an error if any region contains
            zero reads, as this may simply be a matter of low read depth.
        motifs: types of modification to extract data for. Motifs are specified as
            {DNA_sequence},{position_of_modification}. For example, a methylated adenine is specified
            as 'A,0' and CpG methylation is specified as 'CG,0'.
        window_size: An optional parameter for creating centered windows for the provided regions.
            If provided, all regions will be adjusted to be the same size and centered. If not provided,
            all regions should already be the same size, or there should be only one.
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        regions_5to3prime: True means negative strand regions get flipped, False means no flipping
        sort_by: Read properties by which to sort, either one string or a list of strings. Options
            include chromosome, region_start, region_end, read_start, read_end, and motif. More to
            be added in future.
        thresh: A modification calling threshold. If the .h5 is already modification-called, this does
            nothing. If the .h5 files is not modification-called, i.e. its modification data is saved
            as floating point array, thresh must be provided to have valid binary outputs.
        relative: If True, modification coordinates are specified relative to their respective regions
            in the genomes, centered at the center of the region. If False, absolute coordinates are provided.
            There is not currently a check for all reads being on the same chromosome if relative=False, but
            this could create unexpected behaviour for a the standard visualizations.
        quiet: silences progress bars (currently unused)
        cores: cores across which to parallelize processes (currently unused)
        subset_parameters: Parameters to pass to the utils.random_sample() method, to subset the
            reads to be returned. If not None, at least one of n or frac must be provided. The array
            parameter should not be provided here.
        random_sample_n_reads: Optional global random sample size by unique read name.
        min_read_length_bp: Optional minimum read length filter (bp) applied before loading vectors.

    Returns:
        Returns a tuple of three arrays, of length (N_READS * len(mod_names)), and a dict of regions.
        The arrays contain the following:
        * positions at which the specified modification was found in a read, after a binary call
        * unique integer ID for the read for each modification position. These integers are ordered
            based on the specified sorting.
        * modification represented by the positions, in the motif format
        The regions_dict contains the following:
        * keys: chromosomes/contigs
        * values: lists of tuples in the format (region_start,region_end)
        For example, if called on a dataset with a single read and two modification types, each array would have two entries. The unique IDs would be the same, as both entries would represent the same single read. The mods and positions would be different, as they would extact different mods.
    """
    file = Path(file)
    if file.suffix == ".h5" or file.suffix == ".hdf5":
        sorted_read_data_converted, datasets, regions_dict = read_vectors_from_hdf5(
            file=file,
            motifs=motifs,
            regions=regions,
            window_size=window_size,
            single_strand=single_strand,
            sort_by=sort_by,
            quiet=quiet,
            cores=cores,
            subset_parameters=subset_parameters,
            random_sample_n_reads=random_sample_n_reads,
            min_read_length_bp=min_read_length_bp,
        )
        read_name_index = datasets.index("read_name")
        mod_vector_index = datasets.index("mod_vector")
        motif_index = datasets.index("motif")
        region_start_index = datasets.index("region_start")
        region_end_index = datasets.index("region_end")
        read_start_index = datasets.index("read_start")
        region_strand_index = datasets.index("region_strand")

        if len(sorted_read_data_converted) == 0:
            return (
                np.array([], dtype=int),
                np.array([], dtype=int),
                np.array([], dtype=str),
                regions_dict,
            )

        # Check this .h5 file was created with a threshold, i.e. that the mod calls are binarized
        if thresh is None:
            if not (sorted_read_data_converted[0][mod_vector_index].dtype == np.bool_):
                raise ValueError(
                    "No threshold has been applied to this .h5 single read data. You must provide a threshold using the thresh parameter in order to extract binarized modification arrays."
                )
        else:
            thresh = utils.adjust_threshold(thresh)

        read_ints_list = []
        mod_coords_list = []
        motifs_list = []

        read_ids_by_region_key: dict[tuple[str, int, int, str], int] = {}
        read_ints = np.empty(len(sorted_read_data_converted), dtype=int)
        next_read_int = 0
        for idx, read_data in enumerate(sorted_read_data_converted):
            read_key = (
                read_data[read_name_index],
                read_data[region_start_index],
                read_data[region_end_index],
                read_data[region_strand_index],
            )
            read_int = read_ids_by_region_key.get(read_key)
            if read_int is None:
                read_int = next_read_int
                read_ids_by_region_key[read_key] = read_int
                next_read_int += 1
            read_ints[idx] = read_int

        for read_int, read_data in zip(
            read_ints, sorted_read_data_converted, strict=False
        ):
            if thresh is None:
                mod_pos_in_read = np.flatnonzero(read_data[mod_vector_index])
            else:
                mod_pos_in_read = np.flatnonzero(read_data[mod_vector_index] > thresh)

            if relative:
                if regions_5to3prime and read_data[region_strand_index] == "-":
                    # Here we want to show the regions each along their 5 prime to 3 prime direction
                    # This means that negative strand regions need to be flipped
                    mod_pos_record = -(
                        mod_pos_in_read
                        + read_data[read_start_index]
                        - (read_data[region_start_index] + read_data[region_end_index])
                        // 2
                    )
                else:
                    # This is the default case: just make the coordinates relative using
                    # the reference genome coordinate system. Normal, easy, chill, nice.
                    mod_pos_record = (
                        mod_pos_in_read
                        + read_data[read_start_index]
                        - (read_data[region_start_index] + read_data[region_end_index])
                        // 2
                    )
            else:
                # If we aren't using relative coordinates, then I think the 5prime to 3prime argument
                # can just be ignored, and I think it's nicer if that's silent - less clutter in the output
                # Basically if you are keeping different regions separate using other metadata (such as
                # just keeping their actual real genomic coordinates) it is superfluous to do the 5to3 flip.
                mod_pos_record = mod_pos_in_read + read_data[read_start_index]

            mod_coords_list += list(mod_pos_record)
            read_ints_list += [read_int] * len(mod_pos_record)
            motifs_list += [read_data[motif_index]] * len(mod_pos_record)

        return (
            np.array(mod_coords_list),
            np.array(read_ints_list),
            np.array(motifs_list),
            regions_dict,
        )

    else:
        raise ValueError(
            f"File {file} does not have a recognized extension for single read data."
        )


""" TEMPORARY STUB VARS """
STUB_HALFSIZE = 500
STUB_N_READS = 500


def reads_from_fake(
    file: Path, regions: Path, motifs: list[str]
) -> tuple[list[np.ndarray], np.ndarray[int], np.ndarray[str], dict]:
    """
    Helper function to support testing.

    TODO: What does the bed file represent in this method? This one is breaking my brain a bit.
    TODO: Variable names in this method stink.
    TODO: Currently assumes mod calling (thresholding probabilities) was already performed elsewhere

    Args:
        file: Path to file containing modification data for single reads
        bed_file: Path to bed file specifying regions (WHAT DO THESE REPRESENT???)
        mod_names: types of modification to extract data for

    Returns:
        Returns three parallel arrays, of length (N_READS * len(mod_names)), containing the following for each index:
        * array of positions at which the specified modification was found in a read
        * unique integer ID for the read
        * modification represented by the positions
        For example, if called on a dataset with a single read and two modification types, each array would have two entries. The unique IDs would be the same, as both entries would represent the same single read. The mods and positions would be different, as they would extact different mods.
    """
    reads = []
    read_names = []
    mods = []
    for mod_name in motifs:
        match mod_name:
            case "A,0":
                mod_reads = [
                    test_data.fake_read_mod_positions(STUB_HALFSIZE, "peak", 0.7)
                    for _ in range(STUB_N_READS)
                ]
            case "CG,0":
                mod_reads = [
                    test_data.fake_read_mod_positions(
                        STUB_HALFSIZE, "inverse_peak", 0.4
                    )
                    for _ in range(STUB_N_READS)
                ]
            case _:
                raise ValueError(f"No stub settings for requested mod {mod_name}")
        reads += mod_reads
        read_names.append(np.arange(len(mod_reads)))
        mods.append([mod_name] * len(mod_reads))

    read_names = np.concatenate(read_names)
    mods = np.concatenate(mods)
    return reads, read_names, mods, {}


# def convert_bytes(item):
#     """Convert bytes to string if item is bytes, otherwise return as is."""
#     if isinstance(item, bytes):
#         return item.decode()
#     return item


def convert_bytes_to_strings(tup):
    """
    Helper function for single read loading.
    Convert all bytes elements in a tuple to strings.
    """
    return tuple(item.decode() if isinstance(item, bytes) else item for item in tup)
    # tuple(convert_bytes(item) for item in tup)


def adjust_mod_probs_in_arrays(mod_array, val_array):
    """
    Helper function to correct for an idiosyncracy in modkit single-read parsing wherein 0-255
    "mod quality" values are parsed as floating-point values from 1/512 to 511/512.
    """
    mod_array[np.flatnonzero(val_array)] += 1 / 512
    return mod_array


def adjust_mod_probs_in_tuples(tup, mod_idx, val_idx):
    """
    Helper function to apply mod prob adjustments
    """
    return tuple(
        item if index != mod_idx else adjust_mod_probs_in_arrays(item, tup[val_idx])
        for index, item in enumerate(tup)
    )


def binary_to_np_array(compressed_bytes, dtype, decompressor, binarized, int8tofloat):
    """
    Helper function to decompress binary data to boolean or floating point arrays
    """
    if binarized:
        return np.frombuffer(decompressor(compressed_bytes), dtype=dtype).astype(bool)
    elif int8tofloat:
        return (
            (np.frombuffer(decompressor(compressed_bytes), dtype=dtype).astype(float))
            / 256
        ).astype(np.float16)
    else:
        # TODO: Verify whether this branch can be merged with the "binarized" branch in the long run
        return np.frombuffer(decompressor(compressed_bytes), dtype=dtype).astype(bool)


def retrieve_h5_data(h5, dataset, indices, compressed, dtype, decompressor, binarized):
    """
    Load the requested dataset from the h5 file at the relevant indices.

    For compressed vector data, decompress each dataset element to numpy array.
    """
    if compressed:
        # Determine if int8tofloat should be applied
        int8tofloat = "mod_vector" in dataset
        # Logic for compressed data
        loaded_uint8_list = h5[dataset][list(indices)]
        return [
            binary_to_np_array(
                loaded_uint8.tobytes(), dtype, decompressor, binarized, int8tofloat
            )
            for loaded_uint8 in loaded_uint8_list
        ]
    else:
        # Logic for non-compressed data
        return h5[dataset][list(indices)]
