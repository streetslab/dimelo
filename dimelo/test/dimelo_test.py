import os
import unittest

import pandas as pd

from dimelo.functions import PrintDictionaryToTab, SaveMAPQHistogram

# from dimelo.parse_bam import parse_ont_bam
from dimelo.test import DiMeLoTestCase
from dimelo.test.helper.test_utils import (
    create_methylation_objects,
    extract_methylation_data_subset,
)


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
        This is testing correct positions of bases are extracted for a test read
        """
        all_tests = create_methylation_objects()

        # test with one read in particular: 360915ea-cab1-4bf4-85ea-2eaa3e79e5a8
        read_name = "360915ea-cab1-4bf4-85ea-2eaa3e79e5a8"

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

        sub_meg_center_ACG_KEY = [
            199.0,
            204.0,
            209.0,
            212.0,
            213.0,
            218.0,
            224.0,
            225.0,
            227.0,
            233.0,
            246.0,
            249.0,
            258.0,
            259.0,
            260.0,
            264.0,
            266.0,
            267.0,
            268.0,
            269.0,
            273.0,
            274.0,
            289.0,
            294.0,
            295.0,
            297.0,
            298.0,
            299.0,
            300.0,
            301.0,
            302.0,
            303.0,
            304.0,
            305.0,
            306.0,
            307.0,
            308.0,
            309.0,
            310.0,
            311.0,
            312.0,
            313.0,
            314.0,
            315.0,
            317.0,
            318.0,
            323.0,
            330.0,
            331.0,
            332.0,
            335.0,
            337.0,
            338.0,
            340.0,
            341.0,
            350.0,
            355.0,
            360.0,
            363.0,
            368.0,
            369.0,
            370.0,
            372.0,
            375.0,
            379.0,
            386.0,
            392.0,
            399.0,
            400.0,
            406.0,
            408.0,
            409.0,
            411.0,
            414.0,
            418.0,
            424.0,
            428.0,
            431.0,
            435.0,
            440.0,
            441.0,
            451.0,
            459.0,
            467.0,
            470.0,
            471.0,
            472.0,
            473.0,
            474.0,
            476.0,
            478.0,
            479.0,
            480.0,
            481.0,
            482.0,
            485.0,
            492.0,
            496.0,
            497.0,
            498.0,
            504.0,
            506.0,
            507.0,
            516.0,
            519.0,
            521.0,
            523.0,
            527.0,
            530.0,
            533.0,
            535.0,
            537.0,
            538.0,
            541.0,
            544.0,
            551.0,
            552.0,
            557.0,
            565.0,
            568.0,
            573.0,
            574.0,
            581.0,
            582.0,
            591.0,
            605.0,
            606.0,
            608.0,
            609.0,
            611.0,
            616.0,
            618.0,
            619.0,
            624.0,
            625.0,
            631.0,
            632.0,
            640.0,
            641.0,
            642.0,
            643.0,
            656.0,
            658.0,
            659.0,
            664.0,
            667.0,
            670.0,
            671.0,
            674.0,
            676.0,
            677.0,
            681.0,
            685.0,
            693.0,
            698.0,
            699.0,
            709.0,
            711.0,
            713.0,
            715.0,
            716.0,
            727.0,
            728.0,
            732.0,
            735.0,
            737.0,
            738.0,
            739.0,
            742.0,
            743.0,
            744.0,
            748.0,
            766.0,
            768.0,
            773.0,
            774.0,
            782.0,
            789.0,
            790.0,
            792.0,
            797.0,
            802.0,
            810.0,
            812.0,
            818.0,
            825.0,
            827.0,
            838.0,
            840.0,
            841.0,
            849.0,
            858.0,
            865.0,
            866.0,
            875.0,
            878.0,
            881.0,
            882.0,
            883.0,
            884.0,
            888.0,
            890.0,
            891.0,
            901.0,
            904.0,
            907.0,
            909.0,
            910.0,
            914.0,
            928.0,
            931.0,
            936.0,
            939.0,
            943.0,
            945.0,
            950.0,
            952.0,
            955.0,
            962.0,
            970.0,
            971.0,
            972.0,
            973.0,
            974.0,
            977.0,
            981.0,
            987.0,
            994.0,
            997.0,
        ]
        sub_meg_nocenter_ACG_KEY = [
            7505860.0,
            7505865.0,
            7505870.0,
            7505873.0,
            7505874.0,
            7505879.0,
            7505885.0,
            7505886.0,
            7505888.0,
            7505894.0,
            7505907.0,
            7505910.0,
            7505919.0,
            7505920.0,
            7505921.0,
            7505925.0,
            7505927.0,
            7505928.0,
            7505929.0,
            7505930.0,
            7505934.0,
            7505935.0,
            7505950.0,
            7505955.0,
            7505956.0,
            7505958.0,
            7505959.0,
            7505960.0,
        ]
        sub_wg_center_ACG_KEY = [
            212.0,
            213.0,
            267.0,
            268.0,
            392.0,
            418.0,
            419.0,
            591.0,
            656.0,
            658.0,
            659.0,
            664.0,
            667.0,
            681.0,
            685.0,
            742.0,
            743.0,
            744.0,
            875.0,
            878.0,
            888.0,
            928.0,
            931.0,
            939.0,
            940.0,
            941.0,
            943.0,
            945.0,
            952.0,
            977.0,
            981.0,
            987.0,
            994.0,
        ]
        sub_wg_nocenter_ACG_KEY = [7505873.0, 7505874.0, 7505928.0, 7505929.0]

        assert sub_meg_center_ACG["pos"].tolist() == sub_meg_center_ACG_KEY
        assert sub_meg_nocenter_ACG["pos"].tolist() == sub_meg_nocenter_ACG_KEY
        assert sub_wg_center_ACG["pos"].tolist() == sub_wg_center_ACG_KEY
        assert sub_wg_nocenter_ACG["pos"].tolist() == sub_wg_nocenter_ACG_KEY

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
        This is testing correct probabilities of methylation are extracted for a test read
        """
        all_tests = create_methylation_objects()

        # test with one read in particular: 360915ea-cab1-4bf4-85ea-2eaa3e79e5a8
        read_name = "360915ea-cab1-4bf4-85ea-2eaa3e79e5a8"

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

        sub_meg_center_ACG_KEY = [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            18,
            22,
            26,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            3,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            3,
            0,
            0,
            96,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            236,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        sub_meg_nocenter_ACG_KEY = [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        sub_wg_center_ACG_KEY = [
            12,
            18,
            79,
            32,
            74,
            144,
            148,
            70,
            166,
            136,
            58,
            104,
            205,
            170,
            202,
            209,
            20,
            31,
            13,
            78,
            124,
            42,
            57,
            161,
            202,
            143,
            73,
            23,
            247,
            251,
            67,
            168,
            181,
        ]
        sub_wg_nocenter_ACG_KEY = [12, 18, 79, 32]

        assert sub_meg_center_ACG["prob"].tolist() == sub_meg_center_ACG_KEY
        assert (
            sub_meg_nocenter_ACG["prob"].tolist() == sub_meg_nocenter_ACG_KEY
        )
        assert sub_wg_center_ACG["prob"].tolist() == sub_wg_center_ACG_KEY
        assert sub_wg_nocenter_ACG["prob"].tolist() == sub_wg_nocenter_ACG_KEY

    def test_parse_ont_bam_fwd_rev(self):
        """
        check both forward and reverse parsing work
        """
        assert 1 == 1


if __name__ == "__main__":
    unittest.main()
