from pathlib import Path

from dimelo.test import RelativePath

"""
This module contains the test cases that will be run through most dimelo_test submodules. To add new test cases, simple create new dict entries
in test_matrix with descriptive names. Each dict entry must contain a tuple of two dicts, the first containing kwargs and the second containing
test result targets. dimelo_test module run the applicable kwargs through different dimelo modules to test functionality, comparing results against
those stored in the results entry. The output targets can be populated by running generate_targets.py.
"""

# Base input and output directories
test_data_dir = Path("./data")
output_dir = test_data_dir / "test_targets"

region = "chr1:114357437-114359753"  # this region is just one of the ones from the peaks bed file
# 'chr1:9167177-9169177' # nothing is really different about this one; I swapped at one point while debugging

# Paths to input files
ctcf_bam_file = test_data_dir / "ctcf_demo.sorted.bam"
ctcf_target_regions = RelativePath(test_data_dir / "ctcf_demo_peak.bed")
ctcf_off_target_regions = RelativePath(test_data_dir / "ctcf_demo_not_peak.bed")

ctcf_bam_file_updated = RelativePath("./output/ctcf_demo.updated.bam")
output_dir = RelativePath(output_dir)

test_matrix = {
    "megalodon_peaks_190": (
        # input kwargs
        {
            "input_file": ctcf_bam_file_updated,
            "output_name": "megalodon_peaks_190",
            "output_directory": output_dir,
            "regions": [ctcf_target_regions, ctcf_off_target_regions],
            "motifs": ["A,0", "CG,0"],
            "thresh": 190,
            "window_size": 5000,
            "sort_by": ["read_start", "read_name", "motif"],
            "smooth_window": 1,
            "title": "megalodon_peaks_190",
            "single_strand": False,
            "regions_5to3prime": False,
            "chunk_size": 1_000_000,
            "cores": 1,
        },
        # outputs dict function:values
        {},  # populated by generate_targets.py
    ),
    "megalodon_single_190": (
        # input kwargs
        {
            "input_file": ctcf_bam_file_updated,
            "output_name": "megalodon_single_190",
            "output_directory": output_dir,
            "regions": region,
            "motifs": ["A,0", "CG,0"],
            "thresh": 190,
            "window_size": None,
            "sort_by": ["read_start", "read_name", "motif"],
            "smooth_window": 10,
            "title": "megalodon_single_190",
            "single_strand": False,
            "regions_5to3prime": False,
            "chunk_size": 100,
            "cores": 2,
        },
        # outputs dict function:values
        {},  # populated by generate_targets.py
    ),
    "megalodon_single_and_peaks_190": (
        # input kwargs
        {
            "input_file": ctcf_bam_file_updated,
            "output_name": "megalodon_single_and_peaks_190",
            "output_directory": output_dir,
            "regions": [region, ctcf_target_regions, ctcf_off_target_regions],
            "motifs": ["A,0", "CG,0"],
            "thresh": 190,
            "window_size": 5000,
            "sort_by": ["read_start", "read_name", "motif"],
            "smooth_window": 100,
            "title": "megalodon_single_and_peaks_190",
            "single_strand": True,
            "regions_5to3prime": True,
            "chunk_size": 10_000,
            "cores": 3,
        },
        # outputs dict function:values
        {},  # populated by generate_targets.py
    ),
    "megalodon_peaks_nothresh": (
        # input kwargs
        {
            "input_file": ctcf_bam_file_updated,
            "output_name": "megalodon_peaks_nothresh",
            "output_directory": output_dir,
            "regions": [ctcf_target_regions, ctcf_off_target_regions],
            "motifs": ["A,0", "CG,0"],
            "thresh": None,
            "window_size": 5000,
            "sort_by": ["read_start", "read_name", "motif"],
            "smooth_window": 100,
            "title": "megalodon_peaks_nothresh",
            "single_strand": True,
            "regions_5to3prime": False,
            "chunk_size": 1000,
            "cores": 4,
        },
        # outputs dict function:values
        {},  # populated by generate_targets.py
    ),
    "megalodon_single_nothresh": (
        # input kwargs
        {
            "input_file": ctcf_bam_file_updated,
            "output_name": "megalodon_single_nothresh",
            "output_directory": output_dir,
            "regions": region,
            "motifs": ["A,0", "CG,0"],
            "thresh": None,
            "window_size": 5000,
            "sort_by": ["read_start", "read_name", "motif"],
            "smooth_window": 1,
            "title": "megalodon_single_nothresh",
            "single_strand": False,
            "regions_5to3prime": True,
            "chunk_size": 100,
            "cores": None,
        },
        # outputs dict function:values
        {},  # populated by generate_targets.py
    ),
}
