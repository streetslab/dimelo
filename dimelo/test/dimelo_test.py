import unittest

from dimelo.test import DiMeLoTestCase

# import pandas as pd


# from dimelo.test.helper.test_utils import (
#  create_methylation_objects,
# extract_methylation_data_subset,
# )


class TestDiMeLo(DiMeLoTestCase):
    def test_parse_ont_bam_number(self):
        """
        This function calls helper functions:
        - parse_ont_bam_by_window
        - get_modified_reference_positions
        - get_mod_reference_positions_by_mod
        Therefore, tests all four total functions in parse_bam.py.
        Tests under various conditions:
        - Bam: megalodon or guppy+winnowmap
        - Basemod: A, CG or A+CG
        - Center: True or False
        For reference, read names extracted are:
        'f2f51b7c-7818-4771-9101-258a6ee8ecb4',
        '35dbe94b-025a-4d0c-94e5-724bb2bdba0d',
        '360915ea-cab1-4bf4-85ea-2eaa3e79e5a8'
        This is testing correct number of bases are extracted
        """
        # all_tests = create_methylation_objects()

        # test extracting correct number of bases
        # placeholder for now
        assert 1 == 1

    def test_parse_ont_bam_position(self):
        """
        This function calls helper functions:
        - parse_ont_bam_by_window
        - get_modified_reference_positions
        - get_mod_reference_positions_by_mod
        Therefore, tests all four total functions in parse_bam.py.
        Tests under various conditions:
        - Bam: megalodon or guppy+winnowmap
        - Basemod: A, CG or A+CG
        - Center: True or False
        For reference, read names extracted are:
        'f2f51b7c-7818-4771-9101-258a6ee8ecb4',
        '35dbe94b-025a-4d0c-94e5-724bb2bdba0d',
        '360915ea-cab1-4bf4-85ea-2eaa3e79e5a8'
        Test with a forward read and a reverse read.
        This is testing correct positions of bases are extracted for a test read
        """
        # placeholder for now
        assert 1 == 1

    # test qualities parsed correctly
    def test_parse_ont_bam_probability(self):
        """
        This function calls helper functions:
        - parse_ont_bam_by_window
        - get_modified_reference_positions
        - get_mod_reference_positions_by_mod
        Therefore, tests all four total functions in parse_bam.py.
        Tests under various conditions:
        - Bam: megalodon or guppy+winnowmap
        - Basemod: A, CG or A+CG
        - Center: True or False
        For reference, read names extracted are:
        'f2f51b7c-7818-4771-9101-258a6ee8ecb4',
        '35dbe94b-025a-4d0c-94e5-724bb2bdba0d',
        '360915ea-cab1-4bf4-85ea-2eaa3e79e5a8'
        Test with a forward read and a reverse read.
        This is testing correct probabilities of methylation are extracted for a test read
        """

        # placeholder for now
        assert 1 == 1

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
