import subprocess
import unittest

import dimelo as dm
from dimelo.test import DiMeLoTestCase

# import pandas as pd


# from dimelo.test.helper.test_utils import (
#  create_methylation_objects,
# extract_methylation_data_subset,
# )


class TestParseBam(DiMeLoTestCase):
    def test_parse_bam_bedFile(self):
        """Tests parsing a bam file into a database, specifying windows using a bed file.

        Notes:
            - cores set to 1 to ensure that hashes come out the same each time
            - thresholds set super low to ensure that meaningful rows are inserted for all modifications
        """
        # TODO: Is this a reasonable way to specify input files? Where is this intended to be run from?
        dm.parse_bam(
            fileName="dimelo/test/data/mod_mappings_subset.bam",
            sampleName="test",
            outDir=str(self.outDir),
            bedFile="dimelo/test/data/test.bed",
            basemod="A+CG",
            center=True,
            windowSize=500,
            threshA=1,
            threshC=1,
            extractAllBases=False,
            cores=1,
        )
        # TODO: When implemented elsewhere, replace this explicit database path with a modular call
        database_path = self.outDir / "mod_mappings_subset.db"

        # Check whether database contents are the same as expected, using sqlite3 .sha3sum command
        try:
            db_hash_output = subprocess.run(
                ["sqlite3", database_path, ".sha3sum"], capture_output=True,
                check=True
            )
            # db_hash_output = subprocess.run(
            #     ["sqlite3", "-version"], capture_output=True,
            #     check=True
            # )
            # raise subprocess.CalledProcessError(1, 'test', db_hash_output.stdout, db_hash_output.stderr)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"CalledProcessError: {e.stdout}, {e.stderr}")
        self.assertEqual(
            db_hash_output.stdout,
            b"d9f626468e4221a30318639f19dc2cd8a3f00fa6e12a5280c1b58204\n",
        )

    def test_parse_bam_region(self):
        """Tests parsing a bam file into a database, specifying windows using a region string.

        Notes:
            - cores set to 1 to ensure that hashes come out the same each time
            - thresholds set super low to ensure that meaningful rows are inserted for all modifications
        """
        # TODO: Is this a reasonable way to specify input files? Where is this intended to be run from?
        dm.parse_bam(
            fileName="dimelo/test/data/mod_mappings_subset.bam",
            sampleName="test",
            outDir=str(self.outDir),
            region="chr1:2907273-2909473",
            basemod="A+CG",
            threshA=1,
            threshC=1,
            cores=1,
        )
        # TODO: When implemented elsewhere, replace this explicit database path with a modular call
        database_path = self.outDir / "mod_mappings_subset.db"
        # Check whether database contents are the same as expected, using sqlite3 .sha3sum command
        db_hash_output = subprocess.run(
            ["sqlite3", database_path, ".sha3sum"], capture_output=True
        )
        self.assertEqual(
            db_hash_output.stdout,
            b"ae1b8af9121734a96dbcba0df7c0369b14e995a9cf7163a0e9769af6\n",
        )

    def test_parse_bam_bedFile_region_mutual_exclusion(self):
        """Verifies that bedFile and region arguments remain mutually exclusive."""
        # TODO: Is this a reasonable way to specify input files? Where is this intended to be run from?
        with self.assertRaises(RuntimeError):
            dm.parse_bam(
                fileName="dimelo/test/data/mod_mappings_subset.bam",
                sampleName="test",
                outDir=str(self.outDir),
                bedFile="dimelo/test/data/test.bed",
                region="chr1:2907273-2909473",
            )

    def test_parse_bam_region_center_incompatible(self):
        """Verifies that region and center arguments are incompatible."""
        # TODO: Is this a reasonable way to specify input files? Where is this intended to be run from?
        with self.assertRaises(RuntimeError):
            dm.parse_bam(
                fileName="dimelo/test/data/mod_mappings_subset.bam",
                sampleName="test",
                outDir=str(self.outDir),
                region="chr1:2907273-2909473",
                center=True,
            )


# class TestPlotEnrichment(DiMeLoTestCase):
#     """
#     TODO:
#         - exactly one of fileNames and bedFiles must be of length one, and the other must be of length 2
#         - only one of the basemods is valid at a time
#     """
#     def test_plot_enrichment_2_bams(self):
#         dm.plot_enrichment(
#             fileNames=["dimelo/test/data/mod_mappings_subset.bam", "dimelo/test/data/mod_mappings_subset.bam"],
#             sampleNames=["test1", "test2"],
#             bedFiles="dimelo/test/data/test.bed",
#             basemod="CG",
#             outDir=str(self.outDir),
#             threshC=129
#         )


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
