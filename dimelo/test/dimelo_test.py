import subprocess
import unittest
from pathlib import Path

import dimelo as dm
from dimelo.test import DiMeLoTestCase

# import pandas as pd


# from dimelo.test.helper.test_utils import (
#  create_methylation_objects,
# extract_methylation_data_subset,
# )

"""
Inputs

NOTE: Changing any of these may require recomputing the hard-coded database
hashes in the tests below.
"""
# TODO: Is this a reasonable way to specify input files? Where is this intended to be run from?
input_bams = [
    Path("dimelo/test/data/mod_mappings_subset.bam"),
    Path("dimelo/test/data/winnowmap_guppy_merge_subset.bam"),
]
input_sample_names = ["test1", "test2"]
input_bed = Path("dimelo/test/data/test.bed")

"""
Outputs
"""
# TODO: When implemented elsewhere, replace these explicit database paths with modular calls
output_dbs = [f.with_suffix(".db").name for f in input_bams]


def db_hash(db_path: Path) -> bytes:
    """Computes a hash of the sqlite db at the specified path."""
    # Use the sqlite3 .sha3sum command to compute hash
    sqlite_output = subprocess.run(
        ["sqlite3", db_path, ".sha3sum"], capture_output=True
    )
    return sqlite_output.stdout


class TestParseBam(DiMeLoTestCase):
    def test_parse_bam_bedFile(self):
        """Tests parsing a bam file into a database, specifying windows using a bed file.

        Notes:
            - cores set to 1 to ensure that hashes come out the same each time
            - thresholds set super low to ensure that meaningful rows are inserted for all modifications
        """
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]
        sample_name = input_sample_names[bam_idx]
        bed_file = input_bed

        # Outputs
        db_file = output_dbs[bam_idx]

        dm.parse_bam(
            fileName=str(bam_file),
            sampleName=sample_name,
            outDir=str(self.outDir),
            bedFile=str(bed_file),
            basemod="A+CG",
            center=True,
            windowSize=500,
            threshA=1,
            threshC=1,
            extractAllBases=False,
            cores=1,
        )
        database_path = self.outDir / db_file
        # Check whether database contents are the same as expected
        self.assertEqual(
            db_hash(database_path),
            b"8a8fe3984448ee6f215d1d71c01a4d8edde7691bcdb0bf28f2c01cd2\n",
        )

    def test_parse_bam_region(self):
        """Tests parsing a bam file into a database, specifying windows using a region string.

        Notes:
            - cores set to 1 to ensure that hashes come out the same each time
            - thresholds set super low to ensure that meaningful rows are inserted for all modifications
        """
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]
        sample_name = input_sample_names[bam_idx]

        # Outputs
        db_file = output_dbs[bam_idx]

        dm.parse_bam(
            fileName=str(bam_file),
            sampleName=sample_name,
            outDir=str(self.outDir),
            region="chr1:2907273-2909473",
            basemod="A+CG",
            threshA=1,
            threshC=1,
            cores=1,
        )
        database_path = self.outDir / db_file
        # Check whether database contents are the same as expected
        self.assertEqual(
            db_hash(database_path),
            b"baf72c009547e8a2e87402db04695d698e648c3856bcc9b8c0c6cf8a\n",
        )

    def test_parse_bam_bedFile_region_mutual_exclusion(self):
        """Verifies that bedFile and region arguments remain mutually exclusive."""
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]
        sample_name = input_sample_names[bam_idx]
        bed_file = input_bed

        with self.assertRaises(RuntimeError):
            dm.parse_bam(
                fileName=str(bam_file),
                sampleName=sample_name,
                outDir=str(self.outDir),
                bedFile=str(bed_file),
                region="chr1:2907273-2909473",
            )

    def test_parse_bam_region_center_incompatible(self):
        """Verifies that region and center arguments are incompatible."""
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]
        sample_name = input_sample_names[bam_idx]

        with self.assertRaises(RuntimeError):
            dm.parse_bam(
                fileName=str(bam_file),
                sampleName=sample_name,
                outDir=str(self.outDir),
                region="chr1:2907273-2909473",
                center=True,
            )


# TODO: More robust qc_report tests
class TestQCReport(DiMeLoTestCase):
    def test_qc_report_one_sample(self):
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]
        sample_name = input_sample_names[bam_idx]

        # Outputs
        db_file = output_dbs[bam_idx]
        qc_report_file = f"{sample_name}_qc_report.pdf"

        dm.qc_report(
            fileNames=str(bam_file),
            sampleNames=sample_name,
            outDir=str(self.outDir),
        )

        self.assertOutputFileExists(db_file)
        self.assertOutputFileExists(qc_report_file)

    def test_qc_report_multi_sample(self):
        dm.qc_report(
            fileNames=[str(f) for f in input_bams],
            sampleNames=input_sample_names,
            outDir=str(self.outDir),
        )

        for db_file in output_dbs:
            self.assertOutputFileExists(db_file)
        for sample_name in input_sample_names:
            self.assertOutputFileExists(f"{sample_name}_qc_report.pdf")


class TestPlotBrowser(DiMeLoTestCase):
    def test_plot_browser_html(self):
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]
        sample_name = input_sample_names[bam_idx]

        # Outputs
        db_file = output_dbs[bam_idx]

        dm.plot_browser(
            fileNames=str(bam_file),
            sampleNames=sample_name,
            region="chr1:2907273-2909473",
            basemod="A+CG",
            outDir=str(self.outDir),
            static=False,
        )

        self.assertOutputFileExists(db_file)

        # TODO: It's difficult to get the output file name. Find a better way to do this check.
        n_html_files = len(list(self.outDir.glob("*.html")))
        self.assertEqual(n_html_files, 1)

        for basemod in ["A", "C"]:
            for plot_type in ["fraction", "total"]:
                rolling_avg_file = (
                    f"{sample_name}_{basemod}_sm_rolling_avg_{plot_type}.pdf"
                )
                self.assertOutputFileExists(rolling_avg_file)

        # TODO: It's difficult to get the output file name. Find a better way to do this check.
        n_bed_files = len(list(self.outDir.glob("*.bed")))
        self.assertEqual(n_bed_files, 2)

    def test_plot_browser_pdf(self):
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]
        sample_name = input_sample_names[bam_idx]

        # Outputs
        db_file = output_dbs[bam_idx]

        dm.plot_browser(
            fileNames=str(bam_file),
            sampleNames=sample_name,
            region="chr1:2907273-2909473",
            basemod="A+CG",
            outDir=str(self.outDir),
            static=True,
        )

        self.assertOutputFileExists(db_file)

        # TODO: It's difficult to get the output file name. Find a better way to do this check.
        n_pdf_files = len(list(self.outDir.glob("*.pdf")))
        self.assertEqual(n_pdf_files, 5)

        for basemod in ["A", "C"]:
            for plot_type in ["fraction", "total"]:
                rolling_avg_file = (
                    f"{sample_name}_{basemod}_sm_rolling_avg_{plot_type}.pdf"
                )
                self.assertOutputFileExists(rolling_avg_file)

        # TODO: It's difficult to get the output file name. Find a better way to do this check.
        n_bed_files = len(list(self.outDir.glob("*.bed")))
        self.assertEqual(n_bed_files, 2)


class TestPlotEnrichment(DiMeLoTestCase):
    def test_plot_enrichment_2_bams(self):
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]

        # Outputs
        db_file = output_dbs[bam_idx]
        pdf_file = f"region_{input_bed.stem}_CG_enrichment_barplot.pdf"

        dm.plot_enrichment(
            fileNames=[str(bam_file), str(bam_file)],
            sampleNames=input_sample_names,
            bedFiles=str(input_bed),
            basemod="CG",
            outDir=str(self.outDir),
            threshC=129,
        )

        self.assertOutputFileExists(db_file)
        self.assertOutputFileExists(pdf_file)

    def test_plot_enrichment_2_beds(self):
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]

        # Outputs
        db_file = output_dbs[bam_idx]
        pdf_file = f"sample_{bam_file.stem}_CG_enrichment_barplot.pdf"

        dm.plot_enrichment(
            fileNames=str(bam_file),
            sampleNames=input_sample_names,
            bedFiles=[str(input_bed), str(input_bed)],
            basemod="CG",
            outDir=str(self.outDir),
            threshC=129,
        )

        self.assertOutputFileExists(db_file)
        self.assertOutputFileExists(pdf_file)

    def test_plot_enrichment_incompatible_bed_bam(self):
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]

        with self.assertRaises(RuntimeError):
            dm.plot_enrichment(
                fileNames=[str(bam_file), str(bam_file)],
                sampleNames=input_sample_names,
                bedFiles=[str(input_bed), str(input_bed)],
                basemod="CG",
                outDir=str(self.outDir),
            )

    def test_plot_enrichment_incompatible_basemod(self):
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]

        with self.assertRaises(RuntimeError):
            dm.plot_enrichment(
                fileNames=[str(bam_file), str(bam_file)],
                sampleNames=input_sample_names,
                bedFiles=str(input_bed),
                basemod="A+CG",
                outDir=str(self.outDir),
            )


class TestPlotEnrichmentProfile(DiMeLoTestCase):
    def test_plot_enrichment_profile_single(self):
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]
        sample_name = input_sample_names[bam_idx]

        # Outputs
        db_file = output_dbs[bam_idx]
        enrichment_plot = f"{sample_name}_A+CG_sm_rolling_avg.pdf"
        single_molecule_plot = f"{sample_name}_A+CG_sm_scatter.png"
        base_count_plots = [
            f"{sample_name}_{basemod}_base_count.png"
            for basemod in ["A", "CG"]
        ]

        dm.plot_enrichment_profile(
            fileNames=str(bam_file),
            sampleNames=sample_name,
            bedFiles=str(input_bed),
            basemod="A+CG",
            outDir=str(self.outDir),
            windowSize=500,
            dotsize=1,
        )

        self.assertOutputFileExists(db_file)
        self.assertOutputFileExists(enrichment_plot)
        self.assertOutputFileExists(single_molecule_plot)
        for f in base_count_plots:
            self.assertOutputFileExists(f)

    def test_plot_enrichment_profile_sample_overlay(self):
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]

        # Outputs
        db_file = output_dbs[bam_idx]
        overlay_plot = f"sample_{bam_file.stem}_A_sm_rolling_avg_overlay.pdf"

        dm.plot_enrichment_profile(
            fileNames=str(bam_file),
            sampleNames=input_sample_names,
            bedFiles=[str(input_bed), str(input_bed)],
            basemod="A",
            outDir=str(self.outDir),
            windowSize=500,
            dotsize=1,
        )

        self.assertOutputFileExists(db_file)
        self.assertOutputFileExists(overlay_plot)

    def test_plot_enrichment_profile_region_overlay(self):
        # Outputs
        overlay_plot = f"region_{input_bed.stem}_A_sm_rolling_avg_overlay.pdf"

        dm.plot_enrichment_profile(
            fileNames=[str(f) for f in input_bams],
            sampleNames=input_sample_names,
            bedFiles=str(input_bed),
            basemod="A",
            outDir=str(self.outDir),
            windowSize=500,
            dotsize=1,
        )

        for db_file in output_dbs:
            self.assertOutputFileExists(db_file)
        self.assertOutputFileExists(overlay_plot)

    def test_plot_enrichment_profile_overlay_incompatible_basemod(self):
        with self.assertRaises(RuntimeError):
            dm.plot_enrichment_profile(
                fileNames=[str(f) for f in input_bams],
                sampleNames=input_sample_names,
                bedFiles=str(input_bed),
                basemod="A+CG",
                outDir=str(self.outDir),
            )

    def test_plot_enrichment_profile_overlay_incompatible_bed_bam(self):
        with self.assertRaises(RuntimeError):
            dm.plot_enrichment_profile(
                fileNames=[str(f) for f in input_bams],
                sampleNames=input_sample_names,
                bedFiles=[str(input_bed), str(input_bed)],
                basemod="A",
                outDir=str(self.outDir),
            )


# class TestDiMeLo(DiMeLoTestCase):

# def test_parse_ont_bam_number(self):
#     """
#     This function calls helper functions:
#     - parse_ont_bam_by_window
#     - get_modified_reference_positions
#     - get_mod_reference_positions_by_mod
#     Therefore, tests all four total functions in parse_bam.py.
#     Tests under various conditions:
#     - Bam: megalodon or guppy+winnowmap
#     - Basemod: A, CG or A+CG
#     - Center: True or False
#     For reference, read names extracted are:
#     'f2f51b7c-7818-4771-9101-258a6ee8ecb4',
#     '35dbe94b-025a-4d0c-94e5-724bb2bdba0d',
#     '360915ea-cab1-4bf4-85ea-2eaa3e79e5a8'
#     This is testing correct number of bases are extracted
#     """
#     # all_tests = create_methylation_objects()

#     # test extracting correct number of bases
#     # placeholder for now
#     assert 1 == 1

# def test_parse_ont_bam_position(self):
#     """
#     This function calls helper functions:
#     - parse_ont_bam_by_window
#     - get_modified_reference_positions
#     - get_mod_reference_positions_by_mod
#     Therefore, tests all four total functions in parse_bam.py.
#     Tests under various conditions:
#     - Bam: megalodon or guppy+winnowmap
#     - Basemod: A, CG or A+CG
#     - Center: True or False
#     For reference, read names extracted are:
#     'f2f51b7c-7818-4771-9101-258a6ee8ecb4',
#     '35dbe94b-025a-4d0c-94e5-724bb2bdba0d',
#     '360915ea-cab1-4bf4-85ea-2eaa3e79e5a8'
#     Test with a forward read and a reverse read.
#     This is testing correct positions of bases are extracted for a test read
#     """
#     # placeholder for now
#     assert 1 == 1

# # test qualities parsed correctly
# def test_parse_ont_bam_probability(self):
#     """
#     This function calls helper functions:
#     - parse_ont_bam_by_window
#     - get_modified_reference_positions
#     - get_mod_reference_positions_by_mod
#     Therefore, tests all four total functions in parse_bam.py.
#     Tests under various conditions:
#     - Bam: megalodon or guppy+winnowmap
#     - Basemod: A, CG or A+CG
#     - Center: True or False
#     For reference, read names extracted are:
#     'f2f51b7c-7818-4771-9101-258a6ee8ecb4',
#     '35dbe94b-025a-4d0c-94e5-724bb2bdba0d',
#     '360915ea-cab1-4bf4-85ea-2eaa3e79e5a8'
#     Test with a forward read and a reverse read.
#     This is testing correct probabilities of methylation are extracted for a test read
#     """

#     # placeholder for now
#     assert 1 == 1

# TODO add md5sum check for plot
# def test_enrich_sm_roi(self):
#     dirPath = self.tmpFile()
#     enrich_sm_roi(
#         "dimelo/test/data/mod_mappings_subset.bam",
#         "test",
#         "dimelo/test/data/test.bed",
#         "A+CG",
#         "Users/annie/Desktop/test",
#         threshA=0, threshC=200,
#         windowSize=500,
#         dotsize=2,
#     )
#     assert os.path.exists(dirPath + "test_A+CG_sm_scatter.png")

# TODO add md5sum check for plot
# def test_browser_sm_roi(self):
#     dirPath = self.tmpFile()
#     browser_sm_roi(
#         ["dimelo/test/data/mod_mappings_subset.bam"],
#         ["test"],
#         "chr1:7504361-7506361",
#         "A+CG",
#         "/Users/annie/Desktop/dimelo_test",
#         threshA=0, threshC=200,
#         static=True,
#     )
#     assert os.path.exists(
#         dirPath + "methylation_browser_chr1_7504361_7506361.html"
#     )


if __name__ == "__main__":
    unittest.main()
