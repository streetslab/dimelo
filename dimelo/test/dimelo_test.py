import os
import unittest

import pandas as pd

from dimelo.functions import PrintDictionaryToTab, SaveMAPQHistogram
from dimelo.test import DiMeLoTestCase
from dimelo.test.helper.test_utils import (
    create_methylation_objects,
    extract_methylation_data_subset,
)

# from dimelo.parse_bam import parse_ont_bam
# from dimelo.visualize import browser_sm_roi, enrich_sm_roi


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
            6,
        )  # mod_mappings_subset.bam, A, center=True
        assert all_tests[1].shape == (
            175,
            6,
        )  # mod_mappings_subset.bam, A, center=False
        assert all_tests[2].shape == (
            68,
            6,
        )  # mod_mappings_subset.bam, CG, center=True
        assert all_tests[3].shape == (
            7,
            6,
        )  # mod_mappings_subset.bam, CG, center=False
        assert all_tests[4].shape == (
            891,
            6,
        )  # mod_mappings_subset.bam, A+CG, center=True
        assert all_tests[5].shape == (
            182,
            6,
        )  # mod_mappings_subset.bam, A+CG, center=False
        assert all_tests[6].shape == (
            48,
            6,
        )  # winnowmap_guppy_merge_subset.bam, A, center=True
        assert all_tests[7].shape == (
            6,
            6,
        )  # winnowmap_guppy_merge_subset.bam, A, center=False
        assert all_tests[8].shape == (
            15,
            6,
        )  # winnowmap_guppy_merge_subset.bam, CG, center=True
        assert all_tests[9].shape == (
            2,
            6,
        )  # winnowmap_guppy_merge_subset.bam, CG, center=False
        assert all_tests[10].shape == (
            63,
            6,
        )  # winnowmap_guppy_merge_subset.bam, A+CG, center=True
        assert all_tests[11].shape == (
            8,
            6,
        )  # winnowmap_guppy_merge_subset.bam, A+CG, center=False

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
        key = pd.read_csv("dimelo/test/data/key.csv")
        sub_meg_center_ACG_KEY_rev = (
            key[["sub_meg_center_ACG_KEY_rev_pos"]].stack().tolist()
        )
        sub_meg_nocenter_ACG_KEY_rev = (
            key[["sub_meg_nocenter_ACG_KEY_rev_pos"]].stack().tolist()
        )
        sub_wg_center_ACG_KEY_rev = (
            key[["sub_wg_center_ACG_KEY_rev_pos"]].stack().tolist()
        )
        sub_wg_nocenter_ACG_KEY_rev = (
            key[["sub_wg_nocenter_ACG_KEY_rev_pos"]].stack().tolist()
        )

        sub_meg_center_ACG_KEY_fwd = (
            key[["sub_meg_center_ACG_KEY_fwd_pos"]].stack().tolist()
        )
        sub_meg_nocenter_ACG_KEY_fwd = (
            key[["sub_meg_nocenter_ACG_KEY_fwd_pos"]].stack().tolist()
        )
        sub_wg_center_ACG_KEY_fwd = (
            key[["sub_wg_center_ACG_KEY_fwd_pos"]].stack().tolist()
        )
        sub_wg_nocenter_ACG_KEY_fwd = (
            key[["sub_wg_nocenter_ACG_KEY_fwd_pos"]].stack().tolist()
        )

        all_tests = create_methylation_objects()

        read_name_rev = "360915ea-cab1-4bf4-85ea-2eaa3e79e5a8"  # reverse read
        read_name_fwd = "f2f51b7c-7818-4771-9101-258a6ee8ecb4"  # forward read

        read_names = [read_name_rev, read_name_fwd]

        for read_name in read_names:
            sub_meg_center_A = extract_methylation_data_subset(
                all_tests, 0, read_name
            )
            sub_meg_nocenter_A = extract_methylation_data_subset(
                all_tests, 1, read_name
            )
            sub_meg_center_CG = extract_methylation_data_subset(
                all_tests, 2, read_name
            )
            sub_meg_nocenter_CG = extract_methylation_data_subset(
                all_tests, 3, read_name
            )
            sub_meg_center_ACG = extract_methylation_data_subset(
                all_tests, 4, read_name
            )
            sub_meg_nocenter_ACG = extract_methylation_data_subset(
                all_tests, 5, read_name
            )
            sub_wg_center_A = extract_methylation_data_subset(
                all_tests, 6, read_name
            )
            sub_wg_nocenter_A = extract_methylation_data_subset(
                all_tests, 7, read_name
            )
            sub_wg_center_CG = extract_methylation_data_subset(
                all_tests, 8, read_name
            )
            sub_wg_nocenter_CG = extract_methylation_data_subset(
                all_tests, 9, read_name
            )
            sub_wg_center_ACG = extract_methylation_data_subset(
                all_tests, 10, read_name
            )
            sub_wg_nocenter_ACG = extract_methylation_data_subset(
                all_tests, 11, read_name
            )

            # check A matches the A in A+CG and CG matches the CG in A+CG
            # To check A, check A calls for 0 = 4, 1 = 5, 6 = 10, 7 = 11
            # To check CG, check CG calls for 2 = 4, 3 = 5, 8 = 10, 9 = 11
            assert (
                sub_meg_center_A["pos"].tolist()
                == sub_meg_center_ACG[sub_meg_center_ACG["mod"] == "A+Y"][
                    "pos"
                ].tolist()
            )
            assert (
                sub_meg_nocenter_A["pos"].tolist()
                == sub_meg_nocenter_ACG[sub_meg_nocenter_ACG["mod"] == "A+Y"][
                    "pos"
                ].tolist()
            )
            assert (
                sub_meg_center_CG["pos"].tolist()
                == sub_meg_center_ACG[sub_meg_center_ACG["mod"] == "C+Z"][
                    "pos"
                ].tolist()
            )
            assert (
                sub_meg_nocenter_CG["pos"].tolist()
                == sub_meg_nocenter_ACG[sub_meg_nocenter_ACG["mod"] == "C+Z"][
                    "pos"
                ].tolist()
            )
            assert (
                sub_wg_center_A["pos"].tolist()
                == sub_wg_center_ACG[sub_wg_center_ACG["mod"] == "A+a"][
                    "pos"
                ].tolist()
            )
            assert (
                sub_wg_nocenter_A["pos"].tolist()
                == sub_wg_nocenter_ACG[sub_wg_nocenter_ACG["mod"] == "A+a"][
                    "pos"
                ].tolist()
            )
            assert (
                sub_wg_center_CG["pos"].tolist()
                == sub_wg_center_ACG[sub_wg_center_ACG["mod"] == "C+m"][
                    "pos"
                ].tolist()
            )
            assert (
                sub_wg_nocenter_CG["pos"].tolist()
                == sub_wg_nocenter_ACG[sub_wg_nocenter_ACG["mod"] == "C+m"][
                    "pos"
                ].tolist()
            )

            # because extraction of mod alone matches in pair, can just check pos correct for pair A+CG
            if read_name == "360915ea-cab1-4bf4-85ea-2eaa3e79e5a8":
                assert (
                    sub_meg_center_ACG["pos"].tolist()
                    == sub_meg_center_ACG_KEY_rev
                )
                assert (
                    sub_meg_nocenter_ACG["pos"].tolist()
                    == sub_meg_nocenter_ACG_KEY_rev
                )
                assert (
                    sub_wg_center_ACG["pos"].tolist()
                    == sub_wg_center_ACG_KEY_rev
                )
                assert (
                    sub_wg_nocenter_ACG["pos"].tolist()
                    == sub_wg_nocenter_ACG_KEY_rev
                )

            if read_name == "f2f51b7c-7818-4771-9101-258a6ee8ecb4":
                assert (
                    sub_meg_center_ACG["pos"].tolist()
                    == sub_meg_center_ACG_KEY_fwd
                )
                assert (
                    sub_meg_nocenter_ACG["pos"].tolist()
                    == sub_meg_nocenter_ACG_KEY_fwd
                )
                assert (
                    sub_wg_center_ACG["pos"].tolist()
                    == sub_wg_center_ACG_KEY_fwd
                )
                assert (
                    sub_wg_nocenter_ACG["pos"].tolist()
                    == sub_wg_nocenter_ACG_KEY_fwd
                )

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

        key = pd.read_csv("dimelo/test/data/key.csv")
        sub_meg_center_ACG_KEY_rev = (
            key[["sub_meg_center_ACG_KEY_rev_prob"]].stack().tolist()
        )
        sub_meg_nocenter_ACG_KEY_rev = (
            key[["sub_meg_nocenter_ACG_KEY_rev_prob"]].stack().tolist()
        )
        sub_wg_center_ACG_KEY_rev = (
            key[["sub_wg_center_ACG_KEY_rev_prob"]].stack().tolist()
        )
        sub_wg_nocenter_ACG_KEY_rev = (
            key[["sub_wg_nocenter_ACG_KEY_rev_prob"]].stack().tolist()
        )

        sub_meg_center_ACG_KEY_fwd = (
            key[["sub_meg_center_ACG_KEY_fwd_prob"]].stack().tolist()
        )
        sub_meg_nocenter_ACG_KEY_fwd = (
            key[["sub_meg_nocenter_ACG_KEY_fwd_prob"]].stack().tolist()
        )
        sub_wg_center_ACG_KEY_fwd = (
            key[["sub_wg_center_ACG_KEY_fwd_prob"]].stack().tolist()
        )
        sub_wg_nocenter_ACG_KEY_fwd = (
            key[["sub_wg_nocenter_ACG_KEY_fwd_prob"]].stack().tolist()
        )

        all_tests = create_methylation_objects()

        read_name_rev = "360915ea-cab1-4bf4-85ea-2eaa3e79e5a8"  # reverse read
        read_name_fwd = "f2f51b7c-7818-4771-9101-258a6ee8ecb4"  # forward read

        read_names = [read_name_rev, read_name_fwd]

        for read_name in read_names:

            sub_meg_center_A = extract_methylation_data_subset(
                all_tests, 0, read_name
            )
            sub_meg_nocenter_A = extract_methylation_data_subset(
                all_tests, 1, read_name
            )
            sub_meg_center_CG = extract_methylation_data_subset(
                all_tests, 2, read_name
            )
            sub_meg_nocenter_CG = extract_methylation_data_subset(
                all_tests, 3, read_name
            )
            sub_meg_center_ACG = extract_methylation_data_subset(
                all_tests, 4, read_name
            )
            sub_meg_nocenter_ACG = extract_methylation_data_subset(
                all_tests, 5, read_name
            )
            sub_wg_center_A = extract_methylation_data_subset(
                all_tests, 6, read_name
            )
            sub_wg_nocenter_A = extract_methylation_data_subset(
                all_tests, 7, read_name
            )
            sub_wg_center_CG = extract_methylation_data_subset(
                all_tests, 8, read_name
            )
            sub_wg_nocenter_CG = extract_methylation_data_subset(
                all_tests, 9, read_name
            )
            sub_wg_center_ACG = extract_methylation_data_subset(
                all_tests, 10, read_name
            )
            sub_wg_nocenter_ACG = extract_methylation_data_subset(
                all_tests, 11, read_name
            )

            # check A matches the A in A+CG and CG matches the CG in A+CG
            # To check A, check A calls for 0 = 4, 1 = 5, 6 = 10, 7 = 11
            # To check CG, check CG calls for 2 = 4, 3 = 5, 8 = 10, 9 = 11
            assert (
                sub_meg_center_A["prob"].tolist()
                == sub_meg_center_ACG[sub_meg_center_ACG["mod"] == "A+Y"][
                    "prob"
                ].tolist()
            )
            assert (
                sub_meg_nocenter_A["prob"].tolist()
                == sub_meg_nocenter_ACG[sub_meg_nocenter_ACG["mod"] == "A+Y"][
                    "prob"
                ].tolist()
            )
            assert (
                sub_meg_center_CG["prob"].tolist()
                == sub_meg_center_ACG[sub_meg_center_ACG["mod"] == "C+Z"][
                    "prob"
                ].tolist()
            )
            assert (
                sub_meg_nocenter_CG["prob"].tolist()
                == sub_meg_nocenter_ACG[sub_meg_nocenter_ACG["mod"] == "C+Z"][
                    "prob"
                ].tolist()
            )
            assert (
                sub_wg_center_A["prob"].tolist()
                == sub_wg_center_ACG[sub_wg_center_ACG["mod"] == "A+a"][
                    "prob"
                ].tolist()
            )
            assert (
                sub_wg_nocenter_A["prob"].tolist()
                == sub_wg_nocenter_ACG[sub_wg_nocenter_ACG["mod"] == "A+a"][
                    "prob"
                ].tolist()
            )
            assert (
                sub_wg_center_CG["prob"].tolist()
                == sub_wg_center_ACG[sub_wg_center_ACG["mod"] == "C+m"][
                    "prob"
                ].tolist()
            )
            assert (
                sub_wg_nocenter_CG["prob"].tolist()
                == sub_wg_nocenter_ACG[sub_wg_nocenter_ACG["mod"] == "C+m"][
                    "prob"
                ].tolist()
            )

            # because extraction of mod alone matches in pair, can just check pos correct for pair A+CG
            if read_name == "360915ea-cab1-4bf4-85ea-2eaa3e79e5a8":
                assert (
                    sub_meg_center_ACG["prob"].tolist()
                    == sub_meg_center_ACG_KEY_rev
                )
                assert (
                    sub_meg_nocenter_ACG["prob"].tolist()
                    == sub_meg_nocenter_ACG_KEY_rev
                )
                assert (
                    sub_wg_center_ACG["prob"].tolist()
                    == sub_wg_center_ACG_KEY_rev
                )
                assert (
                    sub_wg_nocenter_ACG["prob"].tolist()
                    == sub_wg_nocenter_ACG_KEY_rev
                )

            if read_name == "f2f51b7c-7818-4771-9101-258a6ee8ecb4":
                assert (
                    sub_meg_center_ACG["prob"].tolist()
                    == sub_meg_center_ACG_KEY_fwd
                )
                assert (
                    sub_meg_nocenter_ACG["prob"].tolist()
                    == sub_meg_nocenter_ACG_KEY_fwd
                )
                assert (
                    sub_wg_center_ACG["prob"].tolist()
                    == sub_wg_center_ACG_KEY_fwd
                )
                assert (
                    sub_wg_nocenter_ACG["prob"].tolist()
                    == sub_wg_nocenter_ACG_KEY_fwd
                )

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
