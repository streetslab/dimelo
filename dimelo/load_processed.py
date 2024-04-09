from pathlib import Path
from collections import defaultdict
import random
import gzip

import pysam
import numpy as np
import h5py

from . import test_data
from . import utils


def pileup_counts_from_bedmethyl(
    bedmethyl_file: Path,
    motif: str,
    regions: str | Path | list[str | Path] = None,
) -> tuple[int, int]:
    """
    Extract number of modified bases and total number of bases from the given bedmethyl file

    TODO: How to name this method?

    Args:
        bedmethyl_file: Path to bedmethyl file
        regions: Path to bed file specifying regions
        motif: type of modification to extract data for

    Returns:
        tuple containing counts of (modified_bases, total_bases)
    """

    source_tabix = pysam.TabixFile(str(bedmethyl_file))
    # Don't need vectors, just need counts; also not guaranteed that windows are the same length
    valid_base_count = 0
    modified_base_count = 0

    mod_motif = motif.split(",")[0]
    mod_coord_in_motif = motif.split(",")[1]

    if regions is not None:
        # Get counts from the specified regions
        regions_dict = utils.regions_dict_from_input(regions)
        for chromosome, region_list in regions_dict.items():
            for start_coord, end_coord, strand in region_list:
                if chromosome in source_tabix.contigs:
                    for row in source_tabix.fetch(chromosome, start_coord, end_coord):
                        tabix_fields = row.split("\t")
                        pileup_basemod = tabix_fields[3]
                        keep_basemod = False
                        if len(pileup_basemod.split(",")) == 3:
                            pileup_modname, pileup_motif, pileup_mod_coord = (
                                pileup_basemod.split(",")
                            )
                            if (
                                pileup_motif == mod_motif
                                and pileup_mod_coord == mod_coord_in_motif
                                and pileup_modname
                                in utils.BASEMOD_NAMES_DICT[
                                    mod_motif[int(mod_coord_in_motif)]
                                ]
                            ):
                                keep_basemod = True
                        elif len(pileup_basemod.split(",")) == 1:
                            if (
                                pileup_basemod
                                in utils.BASEMOD_NAMES_DICT[
                                    mod_motif[int(mod_coord_in_motif)]
                                ]
                            ):
                                keep_basemod = True
                        else:
                            raise ValueError(
                                f"Unexpected format in bedmethyl file: {row} contains {pileup_basemod} which cannot be parsed."
                            )
                        if keep_basemod:
                            pileup_info = tabix_fields[9].split(" ")
                            valid_base_count += int(pileup_info[0])
                            modified_base_count += int(pileup_info[2])
    else:
        # Get counts from the whole input file
        for row in source_tabix.fetch():
            tabix_fields = row.split("\t")
            pileup_basemod = tabix_fields[3]
            if mod_motif in pileup_basemod and mod_coord_in_motif in pileup_basemod:
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
    window_size: int = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Generate trace for the specified modification aggregated across all regions in the given bed file. Called by profile plotters, can also be used by a user directly.

    TODO: Merge regions / region_str / window_size handling into a unified function somewhere
    TODO: How to name this method?
    TODO: I feel like stuff like this should be shared functionality
    TODO: I _THINK_ this should return a value for every position, with non-modified positions as zeros
    TODO: Currently, this redundantly generates a centered bed file and then never uses it; still doing centering in the actual parsing loop.
        Obvious solution is to just remove the second centering operation and reference the windowed file. But is there a smarter way to structure the code?

    Args:
        bedmethyl_file: Path to bedmethyl file
        regions: Path to bed file specifying centered equal-length regions
        motif: type of modification to extract data for
        window_size: the extent in either direction for windows around the center of regions.

    Returns:
        vector of fraction modifiied bases (e.g. mA/A) calculated for each position; float values between 0 and 1
    """

    source_tabix = pysam.TabixFile(str(bedmethyl_file))

    mod_motif = motif.split(",")[0]
    mod_coord_in_motif = motif.split(",")[1]

    regions_dict = utils.regions_dict_from_input(regions, window_size)

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
                    tabix_fields = row.split("\t")
                    pileup_basemod = tabix_fields[3]
                    keep_basemod = False
                    if len(pileup_basemod.split(",")) == 3:
                        pileup_modname, pileup_motif, pileup_mod_coord = (
                            pileup_basemod.split(",")
                        )
                        if (
                            pileup_motif == mod_motif
                            and pileup_mod_coord == mod_coord_in_motif
                            and pileup_modname
                            in utils.BASEMOD_NAMES_DICT[
                                mod_motif[int(mod_coord_in_motif)]
                            ]
                        ):
                            keep_basemod = True
                    elif len(pileup_basemod.split(",")) == 1:
                        if (
                            pileup_basemod
                            in utils.BASEMOD_NAMES_DICT[
                                mod_motif[int(mod_coord_in_motif)]
                            ]
                        ):
                            keep_basemod = True
                    else:
                        raise ValueError(
                            f"Unexpected format in bedmethyl file: {row} contains {pileup_basemod} which cannot be parsed."
                        )
                    if keep_basemod:
                        pileup_info = tabix_fields[9].split(" ")
                        pileup_coord_relative = int(tabix_fields[1]) - start_coord
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


def convert_bytes(item):
    """Convert bytes to string if item is bytes, otherwise return as is."""
    if isinstance(item, bytes):
        return item.decode()
    return item


def convert_tuple_elements(tup):
    """Convert all bytes elements in a tuple to strings."""
    return tuple(convert_bytes(item) for item in tup)


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


def process_data(h5, dataset, indices, compressed, dtype, decompressor, binarized):
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


def read_vectors_from_hdf5(
    file: str | Path,
    motifs: list[str],
    regions: str | Path | list[str | Path] = None,
    window_size: int = None,
    sort_by: str | list[str] = ["chromosome", "region_start", "read_start"],
    calculate_mod_fractions: bool = True,
) -> (list[tuple], list[str], dict):
    """
    Pulls a list of read data out of an .h5 file containing processed read vectors, formatted
    for read-by-read vector processing downstream use cases.

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
        window_size: An optional parameter for creating centered windows for the provided regions.
            If provided, all regions will be adjusted to be the same size and centered. If not provided,
            all regions should already be the same size, or there should be only one.
        sort_by: Read properties by which to sort, either one string or a list of strings. Options
            include chromosome, region_start, region_end, read_start, read_end, and motif. More to
            be added in future.

    Returns:
        a list of tuples, each tuple containing all datasets corresponding to an individual read that
        was within the specified regions.
        a list of strings, naming the datasets returned.
        a regions_dict, containing lists of (region_start,region_end) coordinates by chromosome/contig.
    """

    with h5py.File(file, "r") as h5:
        datasets = [name for name, obj in h5.items() if isinstance(obj, h5py.Dataset)]
        if "threshold" in h5:
            # we are looking at an .h5 file with the new, much better compressed format that does
            # not know the data type intrinsically for mod and val vectors, so we must check
            readwise_datasets = [
                dataset for dataset in datasets if dataset not in ["threshold"]
            ]
            compressed_binary_datasets = ["mod_vector", "val_vector"]
            threshold_applied_to_h5 = h5["threshold"][()]
            if np.isnan(threshold_applied_to_h5):
                binarized = False
            else:
                binarized = True
        else:
            # backwards compatible with the old h5 file structure
            readwise_datasets = datasets
            compressed_binary_datasets = []
            binarized = True  # in this case all this will do is make it so we don't apply a +1/512 correction to the mod_vector

        read_chromosomes = np.array(h5["chromosome"], dtype=str)
        read_starts = np.array(h5["read_start"])
        read_ends = np.array(h5["read_end"])
        read_motifs = np.array(h5["motif"], dtype=str)
        ref_strands = np.array(h5["strand"], dtype=str)

        if regions is not None:
            regions_dict = utils.regions_dict_from_input(
                regions=regions,
                window_size=window_size,
            )
            read_data_list = []
            for chrom, region_list in regions_dict.items():
                for region_start, region_end, region_strand in region_list:
                    relevant_read_indices = np.flatnonzero(
                        (read_ends > region_start)
                        & (read_starts < region_end)
                        & np.isin(read_motifs, motifs)
                        & (read_chromosomes == chrom)
                        & (
                            (region_strand not in ["+", "-"])
                            | (ref_strands == region_strand)
                        )
                    )
                    read_data_list += list(
                        zip(
                            *(
                                process_data(
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
                        )
                    )
        else:
            regions_dict = None
            relevant_read_indices = np.flatnonzero(np.isin(read_motifs, motifs))
            read_data_list = list(
                zip(
                    *(
                        process_data(
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
                )
            )
    if binarized:
        read_data_converted = [convert_tuple_elements(tup) for tup in read_data_list]
    else:
        read_data_converted = [
            adjust_mod_probs_in_tuples(
                convert_tuple_elements(tup),
                readwise_datasets.index("mod_vector"),
                readwise_datasets.index("val_vector"),
            )
            for tup in read_data_list
        ]

    readwise_datasets += ["region_start", "region_end"]
    # We add region information (start and end; chromosome is already present!) so that it is possible to sort by these
    if calculate_mod_fractions:
        # # Add MOTIF_mod_fraction for each unique motif in the read_data_list
        # unique_motifs = np.unique(read_motifs)
        # Add the MOTIF_mod_fraction entries to the readwise_datasets list for future reference in sorting
        readwise_datasets += [f"{motif}_mod_fraction" for motif in motifs]
        mod_fractions_by_read_name_by_motif = defaultdict(
            lambda: defaultdict(lambda: 0.0)
        )
        for motif in motifs:
            for read_data in read_data_converted:
                if read_data[readwise_datasets.index("motif")] == motif:
                    mod_sum = np.sum(read_data[readwise_datasets.index("mod_vector")])
                    val_sum = np.sum(read_data[readwise_datasets.index("val_vector")])
                    mod_fraction = mod_sum / val_sum if val_sum > 0 else 0
                    mod_fractions_by_read_name_by_motif[
                        read_data[readwise_datasets.index("read_name")]
                    ][motif] = mod_fraction

        read_data_all = []
        for read_data in read_data_converted:
            read_data_all.append(
                tuple(val for val in read_data)
                + tuple(
                    mod_frac
                    for mod_frac in mod_fractions_by_read_name_by_motif[
                        read_data[readwise_datasets.index("read_name")]
                    ].values()
                )
            )
    else:
        read_data_all = read_data_converted
    # Enforce that sort_by is a list
    if not isinstance(sort_by, list):
        sort_by = [sort_by]

    # If 'shuffle' appears anywhere in sort_by, we first shuffle the list
    if "shuffle" in sort_by:
        random.shuffle(read_data_all)

    try:
        sort_by_indices = [
            readwise_datasets.index(sort_item)
            for sort_item in sort_by
            if sort_item != "shuffle"
        ]
    except ValueError as e:
        raise ValueError(
            f"Sorting error. {e}. Datasets include {readwise_datasets}. If you need mod fraction sorting make sure you are not setting calculate_read_fraction to False."
        )

    if len(sort_by_indices) > 0:
        sorted_read_data = sorted(
            read_data_all, key=lambda x: tuple(x[index] for index in sort_by_indices)
        )
    else:
        sorted_read_data = read_data_all

    return sorted_read_data, readwise_datasets, regions_dict


def readwise_binary_modification_arrays(
    file: str | Path,
    motifs: list[str],
    regions: str | Path | list[str | Path],
    window_size: int = None,
    sort_by: str | list[str] = ["chromosome", "region_start", "read_start"],
    thresh: float = None,
    relative: bool = True,
) -> tuple[list[np.ndarray], np.ndarray[int], np.ndarray[str]]:
    """
    Pulls a list of read data out of a file containing processed read vectors, formatted with
    seaborn plotting in mind. Currently we only support .h5 files.

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
            sort_by=sort_by,
        )
        read_name_index = datasets.index("read_name")
        mod_vector_index = datasets.index("mod_vector")
        motif_index = datasets.index("motif")
        region_start_index = datasets.index("region_start")
        region_end_index = datasets.index("region_end")
        read_start_index = datasets.index("read_start")

        # Check that this .h5 file was created with a threshold, i.e. that the mod calls are binarized
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
        # string_to_int = {read_name: index for index, read_name in zip(first_indices,unique_read_names)}
        read_ints = np.array([string_to_int[read_name] for read_name in read_names])
        for read_int, read_data in zip(read_ints, sorted_read_data_converted):
            if thresh is None:
                mod_pos_in_read = np.flatnonzero(read_data[mod_vector_index])
            else:
                mod_pos_in_read = np.flatnonzero(read_data[mod_vector_index] > thresh)

            if relative:
                mod_pos_record = (
                    mod_pos_in_read
                    + read_data[read_start_index]
                    - (read_data[region_start_index] + read_data[region_end_index]) // 2
                )
            else:
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
) -> tuple[list[np.ndarray], np.ndarray[int], np.ndarray[str]]:
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
