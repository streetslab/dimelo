import argparse
import pickle
from pathlib import Path

from cases import test_matrix
from tqdm.auto import tqdm

from dimelo import load_processed, parse_bam
from dimelo.test import DiMeLoParsingTestCase, RelativePath, filter_kwargs_for_func

ref_genome_file = Path(RelativePath("./output/chm13.draft_v1.0.fasta"))
# Base input and output directories
test_data_dir = Path(RelativePath("./data"))
output_dir = test_data_dir / "test_targets"

output_dir.mkdir(exist_ok=True)


def generate_pileup(test_matrix, case_subset):
    print("Generating pileup files...")
    for case in case_subset if case_subset is not None else test_matrix.keys():
        kwargs, results = test_matrix[case]
        kwargs_pileup = filter_kwargs_for_func(parse_bam.pileup, kwargs)
        pileup_file, pileup_regions = parse_bam.pileup(
            **kwargs_pileup,
            ref_genome=ref_genome_file,
        )
        results["pileup"] = (
            RelativePath(pileup_file),
            RelativePath(pileup_regions),
        )


def generate_extract(test_matrix, case_subset):
    print("Generating extract files...")
    for case in case_subset if case_subset is not None else test_matrix.keys():
        kwargs, results = test_matrix[case]
        kwargs_extract = filter_kwargs_for_func(parse_bam.extract, kwargs)
        if "cores" in kwargs_extract:
            del kwargs_extract["cores"]
        extract_file, extract_regions = parse_bam.extract(
            **kwargs_extract,
            ref_genome=ref_genome_file,
            cores=1,
        )
        results["extract"] = (
            RelativePath(extract_file),
            RelativePath(extract_regions),
        )


def generate_pileup_counts_from_bedmethyl(test_matrix, case_subset):
    for case in tqdm(
        case_subset if case_subset is not None else test_matrix.keys(),
        desc="Generating pileup counts",
    ):
        kwargs, results = test_matrix[case]
        results["pileup_counts_from_bedmethyl"] = {}
        kwargs_func = filter_kwargs_for_func(
            load_processed.pileup_counts_from_bedmethyl, kwargs
        )
        for motif in kwargs["motifs"]:
            results["pileup_counts_from_bedmethyl"][motif] = (
                load_processed.pileup_counts_from_bedmethyl(
                    bedmethyl_file=results["pileup"][0],
                    **kwargs_func,
                    motif=motif,
                )
            )


def generate_pileup_vectors_from_bedmethyl(test_matrix, case_subset):
    for case in tqdm(
        case_subset if case_subset is not None else test_matrix.keys(),
        desc="Generating pileup vectors",
    ):
        kwargs, results = test_matrix[case]
        results["pileup_vectors_from_bedmethyl"] = {}
        kwargs_func = filter_kwargs_for_func(
            load_processed.pileup_vectors_from_bedmethyl, kwargs
        )
        for motif in kwargs["motifs"]:
            results["pileup_vectors_from_bedmethyl"][motif] = (
                load_processed.pileup_vectors_from_bedmethyl(
                    bedmethyl_file=results["pileup"][0],
                    **kwargs_func,
                    motif=motif,
                )
            )


def generate_read_vectors_from_hdf5(test_matrix, case_subset):
    for case in tqdm(
        case_subset if case_subset is not None else test_matrix.keys(),
        desc="Generating read vectors",
    ):
        kwargs, results = test_matrix[case]
        extract_file, regions_bed = results["extract"]
        if extract_file is not None and regions_bed is not None:
            kwargs_func = filter_kwargs_for_func(
                load_processed.read_vectors_from_hdf5, kwargs
            )
            read_data_list, datasets, _ = load_processed.read_vectors_from_hdf5(
                file=extract_file,
                **kwargs_func,
            )
            read_data_dict = {}
            # Pull out the data from the first read
            for idx, dataset in enumerate(datasets):
                for read_data in read_data_list:
                    read_data_dict[dataset] = read_data[idx]
                    break
            results["read_vectors_from_hdf5"] = read_data_dict


def main(test_matrix):
    """
    The main function runs applicable generators based on the test_matrix defined in cases.py and kwargs.

    Args:
        test_matrix: the dict containing test cases. This will be modified in-place to contain existing targets
            if generators won't be getting re-run, based on the command line parsed arguments described below.
    """
    # Set up input files, including ref genome download
    DiMeLoParsingTestCase.setup_class()

    parser = argparse.ArgumentParser(
        description="Generate target data from test cases."
    )

    valid_subsets = [
        "pileup",
        "extract",
        "pileup_counts_from_bedmethyl",
        "pileup_vectors_from_bedmethyl",
        "read_vectors_from_hdf5",
    ]

    parser.add_argument(
        "--target-subset",
        nargs="*",
        help=f"""Specify one or more subsets of test targets, separated by spaces (default: all).
                        The following are valid options: {valid_subsets}""",
    )

    parser.add_argument(
        "--case-subset",
        nargs="*",
        help=f"""Specify one or more subsets of cases, separated by spaces (default: all).
                        The following are valid options based on your current cases.py file: {test_matrix.keys()}""",
    )

    parser.add_argument(
        "--initial-target-pickle",
        type=str,
        default=RelativePath("data/test_targets/test_matrix.pickle"),
        help="A test target pickle to start with, when updating only a subset of cases.",
    )

    args = parser.parse_args()

    if args.target_subset is None:
        args.target_subset = valid_subsets
        print(
            f"Running {'all cases' if args.case_subset is None else args.case_subset} through all target generators from scratch based on test_cases.py"
        )
    else:
        if Path(args.initial_target_pickle).exists():
            with open(RelativePath(args.initial_target_pickle), "rb") as file:
                old_test_matrix = pickle.load(file)
                # loop through the old test matrix
                for key, (old_kwargs, old_results) in old_test_matrix.items():
                    # if the old test case is in the new test matrix, bring over the results
                    # either the test targets will be regenerated, and the results replaced, or they won't be regenerated, and we'll need the old results
                    if key in test_matrix:
                        new_kwargs = test_matrix[key][0]
                        # if case will be covered, bring in new kwargs
                        if args.case_subset is None or key in args.case_subset:
                            test_matrix[key] = (new_kwargs, old_results)
                        # if case will not be covered, keep old kwargs
                        else:
                            test_matrix[key] = (old_kwargs, old_results)
                print(
                    f"Running {'all cases' if args.case_subset is None else args.case_subset} through {args.target_subset} to supplement test targets from {args.initial_target_pickle}. Any new tests in cases.py will be added to the test matrix."
                )
        else:
            raise ValueError(
                f"Cannot run subset {args.target_subset} without a pre-existing complete set of generated targets, {args.initial_target_pickle} does not exist. Either specify an --initial-target-pickle that does exist or run without subsetting."
            )

    print("Generating targets for the following test matrix kwargs")
    for test_name, (kwargs, _) in test_matrix.items():
        print(test_name)
        print(kwargs)

    for subset in args.target_subset:
        function_name = f"generate_{subset}"
        globals()[function_name](test_matrix, args.case_subset)

    with open(RelativePath("./data/test_targets/test_matrix.pickle"), "wb") as file:
        pickle.dump(test_matrix, file)


if __name__ == "__main__":
    main(test_matrix)
