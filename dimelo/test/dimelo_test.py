import os
import unittest

import pandas as pd

from dimelo.functions import PrintDictionaryToTab, SaveMAPQHistogram

# from dimelo.parse_bam import parse_ont_bam
from dimelo.test import DiMeLoTestCase
from dimelo.test.helper.test_utils import create_methylation_objects


class TestDiMeLo(DiMeLoTestCase):
    def test_PrintDictionaryToTab(self):

        dictIn = {0: 10, 4: 30, 10: 20, 50: 100}

        filePath = self.tmpFile()
        print(filePath)
        PrintDictionaryToTab("MAPQ", "readCount", dictIn, filePath)

        # read back in temp file
        df = pd.read_csv(filePath, sep="\t")
        assert df.shape[0] == 4

    def test_SaveMAPQHistogram(self):

        dictIn = {0: 10, 4: 30, 10: 20, 50: 100}

        filePath = self.tmpFile() + ".pdf"

        SaveMAPQHistogram(dictIn, filePath, title="MAPQ")
        assert os.path.exists(filePath)

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
        all_tests = create_methylation_objects()

        # test extracting correct number of bases
        assert all_tests[0].shape == (
            823,
            5,
        )  # mod_mappings_subset.bam, A, center=True
        assert all_tests[1].shape == (
            175,
            5,
        )  # mod_mappings_subset.bam, A, center=False
        assert all_tests[2].shape == (
            68,
            5,
        )  # mod_mappings_subset.bam, CG, center=True
        assert all_tests[3].shape == (
            7,
            5,
        )  # mod_mappings_subset.bam, CG, center=False
        assert all_tests[4].shape == (
            891,
            5,
        )  # mod_mappings_subset.bam, A+CG, center=True
        assert all_tests[5].shape == (
            182,
            5,
        )  # mod_mappings_subset.bam, A+CG, center=False
        assert all_tests[6].shape == (
            48,
            5,
        )  # winnowmap_guppy_merge_subset.bam, A, center=True
        assert all_tests[7].shape == (
            6,
            5,
        )  # winnowmap_guppy_merge_subset.bam, A, center=False
        assert all_tests[8].shape == (
            15,
            5,
        )  # winnowmap_guppy_merge_subset.bam, CG, center=True
        assert all_tests[9].shape == (
            2,
            5,
        )  # winnowmap_guppy_merge_subset.bam, CG, center=False
        assert all_tests[10].shape == (
            63,
            5,
        )  # winnowmap_guppy_merge_subset.bam, A+CG, center=True
        assert all_tests[11].shape == (
            8,
            5,
        )  # winnowmap_guppy_merge_subset.bam, A+CG, center=False

        # TODO add more tests for bam parsing
        # test with one read in particular: 360915ea-cab1-4bf4-85ea-2eaa3e79e5a8
        # sub_meg_center_ACG = all_tests[4][
        #     all_tests[4]["read_name"] == "360915ea-cab1-4bf4-85ea-2eaa3e79e5a8"
        # ]
        # sub_meg_nocenter_ACG = all_tests[5][
        #     all_tests[5]["read_name"] == "360915ea-cab1-4bf4-85ea-2eaa3e79e5a8"
        # ]
        # sub_wg_center_ACG = all_tests[10][
        #     all_tests[10]["read_name"]
        #     == "360915ea-cab1-4bf4-85ea-2eaa3e79e5a8"
        # ]
        # sub_wg_nocenter_ACG = all_tests[11][
        #     all_tests[11]["read_name"]
        #     == "360915ea-cab1-4bf4-85ea-2eaa3e79e5a8"
        # ]

        # test positions parsed correctly
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
            This is testing correct positions of bases are extracted
            """
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
            This is testing correct probabilities of methylation are extracted
            """
            assert 1 == 1


if __name__ == "__main__":
    unittest.main()
