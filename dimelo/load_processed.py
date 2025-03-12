import gzip
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import h5py
import numpy as np
import pysam
from tqdm.auto import tqdm

from . import test_data, utils


def process_region(region_string, function_handle, **kwargs):
    """
    process_region simply exists to convert position arguments into keyword arguments to make executor.map work

    Args:
        region_string: passed down with regions keyword
        function_handle: function to call with regions and other kwargs
        **kwargs: all keyword arguments passed to regions_to_list. These must be sufficient for whichever load_processed function
            if being referenced by function_handle
    Returns:
        function_handle return value
    """
    return function_handle(regions=region_string, **kwargs)


def regions_to_list(
    function_handle,
    regions,
    window_size: int | None = None,
    cores: int | None = None,
    **kwargs,
):
    """
    Run any standard load_processed pileup or extract loader loading each region from the region
    specifier into a new element of a list.

    Args:
        function_handle: the loader function you want to run.
        regions: the region specifier
        window_size: window around centers of regions, defaults to None
        cores: process count across which to parallelize. Each individual region will only ever get one core.
        **kwargs: all necessary keyword arguments to pass down to the loader

    Returns:
        List(function_handle return objects per region)
    """
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

    cores_to_run = utils.cores_to_run(cores)

    if cores_to_run > 1:
        with ProcessPoolExecutor(max_workers=cores_to_run) as executor:
            # Use functools.partial to pre-fill arguments
            process_partial = partial(
                process_region, function_handle=function_handle, cores=1, **kwargs
            )

            # Use executor.map without lambda
            results = list(
                tqdm(
                    executor.map(process_partial, region_strings),
                    total=len(region_strings),
                    desc=f"Processing regions in parallel across {cores_to_run}",
                )
            )
    else:
        # Single-threaded fallback
        results = [
            process_region(
                region_string=region, function_handle=function_handle, cores=1, **kwargs
            )
            for region in tqdm(region_strings, desc="Processing regions")
        ]

    return results


def pileup_counts_from_bedmethyl(
    bedmethyl_file: str | Path,
    motif: str,
    regions: str | Path | list[str | Path] | None = None,
    window_size: int | None = None,
    single_strand: bool = False,
    cores: int | None = None,  # currently unused
) -> tuple[int, int]:
    """
    Extract number of modified bases and total number of bases from the given bedmethyl file.
    Called by plotters or by the user.

    This function loops through all the provided regions and pulls those regions up in the input
    sorted and indexed bedmethyl file. For rows within those regions, checks that the motif
    is correct (i.e. sequence context, modified base, mod code, and optionally strand). All
    correct locations are included in the sum counts that get returned.

    If no regions are specified, returns the sum total for the motif of interest across the
    entire bedmethyl file.

    TODO: Consider renaming this method, e.g. counts_from_pileup

    Args:
        bedmethyl_file: Path to bedmethyl file
        regions: Path to bed file specifying regions
        motif: type of modification to extract data for
        window_size: (currently disabled) window around center of region, +-window_size
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        cores: cores across which to parallelize processes (currently unused)

    Returns:
        tuple containing counts of (modified_bases, total_bases)
    """

    source_tabix = pysam.TabixFile(str(bedmethyl_file))
    # Don't need vectors, just need counts; also not guaranteed that windows are the same length
    valid_base_count = 0
    modified_base_count = 0

    parsed_motif = utils.ParsedMotif(motif)

    if regions is not None:
        # Get counts from the specified regions
        regions_dict = utils.regions_dict_from_input(
            regions,
            window_size,
        )
        for chromosome, region_list in regions_dict.items():
            for start_coord, end_coord, strand in region_list:
                # TODO: change to try-except
                if chromosome in source_tabix.contigs:
                    for row in source_tabix.fetch(chromosome, start_coord, end_coord):
                        # TODO Consider using csv module
                        # TODO: probably this whole block should share logic with vectors_from_bedmethyl AND from export module functions
                        tabix_fields = row.split("\t")
                        pileup_basemod = tabix_fields[3]
                        pileup_strand = tabix_fields[5]
                        keep_basemod = False
                        if single_strand and pileup_strand != strand:
                            # This entry is on the wrong strand - skip it
                            continue
                        elif len(pileup_basemod.split(",")) == 3:
                            pileup_modname, pileup_motif, pileup_mod_coord = (
                                pileup_basemod.split(",")
                            )
                            if (
                                pileup_motif == parsed_motif.motif_seq
                                and int(pileup_mod_coord) == parsed_motif.modified_pos
                                and pileup_modname in parsed_motif.mod_codes
                            ):
                                keep_basemod = True
                        elif len(pileup_basemod.split(",")) == 1:
                            if pileup_basemod in parsed_motif.mod_codes:
                                keep_basemod = True
                        else:
                            raise ValueError(
                                f"Unexpected format in bedmethyl file: {row} contains {pileup_basemod} which cannot be parsed."
                            )
                        # TODO: consolidate the above into a function; just do adding outside
                        if keep_basemod:
                            pileup_info = tabix_fields[9].split(" ")
                            valid_base_count += int(pileup_info[0])
                            modified_base_count += int(pileup_info[2])
    else:
        # Get counts from the whole input file
        for row in source_tabix.fetch():
            tabix_fields = row.split("\t")
            pileup_basemod = tabix_fields[3]
            keep_basemod = False
            if len(pileup_basemod.split(",")) == 3:
                pileup_modname, pileup_motif, pileup_mod_coord = pileup_basemod.split(
                    ","
                )
                if (
                    pileup_motif == parsed_motif.motif_seq
                    and int(pileup_mod_coord) == parsed_motif.modified_pos
                    and pileup_modname in parsed_motif.mod_codes
                ):
                    keep_basemod = True
            elif len(pileup_basemod.split(",")) == 1:
                if pileup_basemod in parsed_motif.mod_codes:
                    keep_basemod = True
            else:
                raise ValueError(
                    f"Unexpected format in bedmethyl file: {row} contains {pileup_basemod} which cannot be parsed."
                )
            if keep_basemod:
                pileup_info = tabix_fields[9].split(" ")
                valid_base_count += int(pileup_info[0])
                modified_base_count += int(pileup_info[2])

    return (modified_base_count, valid_base_count)


def counts_from_fake(*args, **kwargs) -> tuple[int, int]:
    """
    Generates a fake set of enrichment counts. Ignores all arguments.

    Returns:
        tuple containing counts of (modified_bases, total_bases)
    """
    window_halfsize = 500
    return test_data.fake_peak_enrichment(halfsize=window_halfsize, peak_height=0.15)


def pileup_vectors_from_bedmethyl(
    bedmethyl_file: str | Path,
    motif: str,
    regions: str | Path | list[str | Path],
    window_size: int | None = None,
    single_strand: bool = False,
    regions_5to3prime: bool = False,
    cores: int | None = None,  # currently unused
) -> tuple[np.ndarray, np.ndarray]:
    """
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

    TODO: Consider renaming this method, e.g. vectors_from_pileup

    Args:
        bedmethyl_file: Path to bedmethyl file
        regions: Path to bed file specifying centered equal-length regions
        motif: type of modification to extract data for
        window_size: the extent in either direction for windows around the center of regions.
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        regions_5to3prime: True means negative strand regions get flipped, False means no flipping
        cores: cores across which to parallelize processes (currently unused)

    Returns:
        tuple containing (modified_base_counts, valid_base_counts)
    """

    source_tabix = pysam.TabixFile(str(bedmethyl_file))

    parsed_motif = utils.ParsedMotif(motif)

    regions_dict = utils.regions_dict_from_input(regions, window_size)

    # Peek at a region to figure out what size the vectors should be
    first_key = next(iter(regions_dict))
    first_tuple = regions_dict[first_key][0]
    region_len = first_tuple[1] - first_tuple[0]

    valid_base_counts = np.zeros(region_len, dtype=int)
    modified_base_counts = np.zeros(region_len, dtype=int)

    for chromosome, region_list in regions_dict.items():
        for start_coord, end_coord, strand in region_list:
            # TODO: This is not used anywhere; disabling for now
            # center_coord = (start_coord+end_coord)//2
            if chromosome in source_tabix.contigs:
                for row in source_tabix.fetch(chromosome, start_coord, end_coord):
                    # TODO: can we consolidate this with pileup_counts_from_bedmethyl?
                    # Just the checks?
                    # TODO: probably this whole block should share logic with counts_from_bedmethyl AND from export functions
                    tabix_fields = row.split("\t")
                    pileup_basemod = tabix_fields[3]
                    pileup_strand = tabix_fields[5]
                    keep_basemod = False
                    if single_strand and pileup_strand.strip() != strand:
                        # We are on the wrong strand, skip the rest of the steps for this row
                        continue
                    elif len(pileup_basemod.split(",")) == 3:
                        pileup_modname, pileup_motif, pileup_mod_coord = (
                            pileup_basemod.split(",")
                        )
                        if (
                            pileup_motif == parsed_motif.motif_seq
                            and int(pileup_mod_coord) == parsed_motif.modified_pos
                            and pileup_modname in parsed_motif.mod_codes
                        ):
                            keep_basemod = True
                    elif len(pileup_basemod.split(",")) == 1:
                        if pileup_basemod in parsed_motif.mod_codes:
                            keep_basemod = True
                    else:
                        raise ValueError(
                            f"Unexpected format in bedmethyl file: {row} contains {pileup_basemod} which cannot be parsed."
                        )
                    if keep_basemod:
                        pileup_info = tabix_fields[9].split(" ")
                        genomic_coord = int(tabix_fields[1])
                        if regions_5to3prime and strand == "-":
                            # We want to flip the coordinates for this region so that it is recorded along the 5 prime to 3 prime direction
                            # This will enable analyses where the orientation of protein binding / transcriptional dynamics / etc is relevant for our pileup signal
                            pileup_coord_relative = end_coord - genomic_coord - 1
                        else:
                            # Normal coordinates are the default. This will be used both for the '+' case and the '.' (no strand specified) case
                            pileup_coord_relative = genomic_coord - start_coord
                        if pileup_coord_relative > region_len:
                            print(
                                f"WARNING: You have specified a region {chromosome}:{start_coord}-{end_coord} that is longer than the first region; the end of the region will be skipped. To make a profile plot with differently-sized region, consider using the window_size parameter to make a profile across centered windows."
                            )
                            break
                        else:
                            valid_base_counts[pileup_coord_relative] += int(
                                pileup_info[0]
                            )
                            modified_base_counts[pileup_coord_relative] += int(
                                pileup_info[2]
                            )
    return modified_base_counts, valid_base_counts


def vector_from_fake(window_size: int, *args, **kwargs) -> np.ndarray:
    """
    Generates a fake peak trace. Ignores all arguments except window_size.

    Args:
        window_size: halfsize of the window; how far the window stretches on either side of the center point

    Returns:
        vector of fraction modified bases calculated for each position; float values between 0 and 1
    """
    return test_data.fake_peak_enrichment_profile(
        halfsize=window_size, peak_height=0.15
    )


def read_vectors_from_hdf5(
    file: str | Path,
    motifs: list[str],
    regions: str | Path | list[str | Path] | None = None,
    window_size: int | None = None,
    single_strand: bool = False,
    sort_by: str | list[str] = ["chromosome", "region_start", "read_start"],
    calculate_mod_fractions: bool = True,
    cores: int | None = None,  # currently unused
    subset_parameters: dict | None = None,
) -> tuple[list[tuple], list[str], dict | None]:
    """
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
        cores: cores across which to parallelize processes (currently unused)
        subset_parameters: Parameters to pass to the utils.random_sample() method, to subset the
            reads to be returned. If not None, at least one of n or frac must be provided. The array
            parameter should not be provided here.

    Returns:
        a list of tuples, each tuple containing all datasets corresponding to an individual read that
        was within the specified regions.
        a list of strings, naming the datasets returned.
        a regions_dict, containing lists of (region_start,region_end) coordinates by chromosome/contig.

    TODO: The way the subsetting is implemented is confusing, in that you need to pass all but one of
        the available parameters.
    """
    with h5py.File(file, "r") as h5:
        datasets: list[str] = [
            name for name, obj in h5.items() if isinstance(obj, h5py.Dataset)
        ]
        if "threshold" in h5:
            # we are looking at an .h5 file with the new, much better compressed format that does
            # not know the data type intrinsically for mod and val vectors, so we must check
            readwise_datasets = [
                dataset for dataset in datasets if dataset not in ["threshold"]
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
                        & np.isin(read_motifs, motifs)
                        & (read_chromosomes == chrom)
                        & (
                            (not single_strand)
                            | (region_strand not in ["+", "-"])
                            | (ref_strands == region_strand)
                        )
                    )
                    if subset_parameters is not None:
                        relevant_read_indices = np.sort(
                            utils.random_sample(
                                relevant_read_indices, **subset_parameters
                            )
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
                        )
                    )
        else:
            regions_dict = None
            relevant_read_indices = np.flatnonzero(np.isin(read_motifs, motifs))
            if subset_parameters is not None:
                relevant_read_indices = np.sort(
                    utils.random_sample(relevant_read_indices, **subset_parameters)
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
                )
            )
    #  We add region information (start, end, and strand; chromosome is already present!)
    # so that it is possible to sort by and process based on these
    readwise_datasets += ["region_start", "region_end", "region_strand"]

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

    if calculate_mod_fractions:
        # Add the MOTIF_mod_fraction entries to the readwise_datasets list for future reference in sorting
        readwise_datasets += [f"{motif}_mod_fraction" for motif in motifs]
        # dict[read_name][motif]=modified fraction of motif in read, float
        mod_fractions_by_read_name_by_motif: defaultdict[
            str, defaultdict[str, float]
        ] = defaultdict(lambda: defaultdict(lambda: 0.0))
        for motif in motifs:
            for read_tuple in read_tuples_processed:
                if read_tuple[readwise_datasets.index("motif")] == motif:
                    mod_sum = np.sum(read_tuple[readwise_datasets.index("mod_vector")])
                    val_sum = np.sum(read_tuple[readwise_datasets.index("val_vector")])
                    mod_fraction = mod_sum / val_sum if val_sum > 0 else 0
                    mod_fractions_by_read_name_by_motif[
                        read_tuple[readwise_datasets.index("read_name")]
                    ][motif] = mod_fraction

        read_tuples_all = []
        for read_tuple in read_tuples_processed:
            read_tuples_all.append(
                tuple(val for val in read_tuple)
                + tuple(
                    mod_frac
                    for mod_frac in mod_fractions_by_read_name_by_motif[
                        read_tuple[readwise_datasets.index("read_name")]
                    ].values()
                )
            )
    else:
        read_tuples_all = read_tuples_processed

    ## Sort the reads

    # Enforce that sort_by is a list
    if not isinstance(sort_by, list):
        sort_by = [sort_by]

    # If 'shuffle' appears anywhere in sort_by, we first shuffle the list
    if "shuffle" in sort_by:
        utils.rng.shuffle(read_tuples_all)

    try:
        sort_by_indices = [
            readwise_datasets.index(sort_item)
            for sort_item in sort_by
            if sort_item != "shuffle"
        ]
    except ValueError as e:
        raise ValueError(
            f"Sorting error. {e}. Datasets include {readwise_datasets}. If you need mod fraction sorting make sure you are not setting calculate_read_fraction to False."
        ) from e

    if len(sort_by_indices) > 0:
        sorted_read_tuples = sorted(
            read_tuples_all, key=lambda x: tuple(x[index] for index in sort_by_indices)
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
    cores: int | None = None,  # currently unused
    subset_parameters: dict | None = None,
) -> tuple[list[np.ndarray], np.ndarray[int], np.ndarray[str], dict | None]:
    """
    Pulls a list of read data out of a file containing processed read vectors, formatted with
    seaborn plotting in mind. Currently we only support .h5 files.

    After running read_vectors_from_hdf5, this function takes the baton to convert the names of
    the sorted reads into integer indices, then goes through the reads and strips down the mod
    vectors to simply a list of modified positions (applying a threshold if one has not already
    been applied). Mod positions are by default expressed relative to the region from which
    the read was identified, allowing for nice plotting, but can also be expressed in absolute
    coordinates. If positions are relative, regions_5to3prime can be used to show all regions
    as upstream-to-downstream along their respective strands.

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
        cores: cores across which to parallelize processes (currently unused)
        subset_parameters: Parameters to pass to the utils.random_sample() method, to subset the
            reads to be returned. If not None, at least one of n or frac must be provided. The array
            parameter should not be provided here.

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
        )
        read_name_index = datasets.index("read_name")
        mod_vector_index = datasets.index("mod_vector")
        motif_index = datasets.index("motif")
        region_start_index = datasets.index("region_start")
        region_end_index = datasets.index("region_end")
        read_start_index = datasets.index("read_start")
        region_strand_index = datasets.index("region_strand")

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

        read_names = np.array(
            [read_data[read_name_index] for read_data in sorted_read_data_converted]
        )
        # TODO: handle the case where a read shows up in more than one different region
        _, unique_first_indices = np.unique(read_names, return_index=True)
        unique_in_order = read_names[np.sort(unique_first_indices)]
        string_to_int = {
            read_name: index for index, read_name in enumerate(unique_in_order)
        }
        read_ints = np.array([string_to_int[read_name] for read_name in read_names])

        for read_int, read_data in zip(read_ints, sorted_read_data_converted):
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
    """Convert all bytes elements in a tuple to strings."""
    return tuple(item.decode() if isinstance(item, bytes) else item for item in tup)
    # tuple(convert_bytes(item) for item in tup)


def adjust_mod_probs_in_arrays(mod_array, val_array):
    mod_array[np.flatnonzero(val_array)] += 1 / 512
    return mod_array


def adjust_mod_probs_in_tuples(tup, mod_idx, val_idx):
    return tuple(
        item if index != mod_idx else adjust_mod_probs_in_arrays(item, tup[val_idx])
        for index, item in enumerate(tup)
    )


def binary_to_np_array(compressed_bytes, dtype, decompressor, binarized, int8tofloat):
    if binarized:
        return np.frombuffer(decompressor(compressed_bytes), dtype=dtype).astype(bool)
    elif int8tofloat:
        return (
            (np.frombuffer(decompressor(compressed_bytes), dtype=dtype).astype(float))
            / 256
        ).astype(np.float16)
    else:
        return np.frombuffer(decompressor(compressed_bytes), dtype=dtype).astype(int)


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
