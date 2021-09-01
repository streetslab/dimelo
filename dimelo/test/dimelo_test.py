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
        Test with a forward read and a reverse read.
        This is testing correct positions of bases are extracted for a test read
        """
        sub_meg_center_ACG_KEY_rev = [
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
        sub_meg_nocenter_ACG_KEY_rev = [
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
        sub_wg_center_ACG_KEY_rev = [
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
        sub_wg_nocenter_ACG_KEY_rev = [
            7505873.0,
            7505874.0,
            7505928.0,
            7505929.0,
        ]

        sub_meg_center_ACG_KEY_fwd = [
            3.0,
            9.0,
            10.0,
            12.0,
            14.0,
            15.0,
            17.0,
            18.0,
            22.0,
            25.0,
            27.0,
            31.0,
            36.0,
            45.0,
            47.0,
            60.0,
            64.0,
            66.0,
            67.0,
            68.0,
            70.0,
            71.0,
            72.0,
            75.0,
            83.0,
            97.0,
            102.0,
            104.0,
            106.0,
            119.0,
            124.0,
            125.0,
            138.0,
            142.0,
            144.0,
            148.0,
            153.0,
            160.0,
            164.0,
            178.0,
            182.0,
            184.0,
            186.0,
            211.0,
            228.0,
            233.0,
            235.0,
            247.0,
            252.0,
            267.0,
            268.0,
            273.0,
            274.0,
            281.0,
            282.0,
            285.0,
            286.0,
            289.0,
            294.0,
            297.0,
            304.0,
            315.0,
            323.0,
            327.0,
            329.0,
            336.0,
            344.0,
            352.0,
            354.0,
            363.0,
            365.0,
            384.0,
            389.0,
            393.0,
            396.0,
            399.0,
            403.0,
            418.0,
            420.0,
            423.0,
            428.0,
            429.0,
            434.0,
            441.0,
            442.0,
            444.0,
            450.0,
            451.0,
            453.0,
            454.0,
            462.0,
            467.0,
            468.0,
            474.0,
            479.0,
            482.0,
            486.0,
            488.0,
            489.0,
            491.0,
            492.0,
            495.0,
            497.0,
            498.0,
            508.0,
            519.0,
            525.0,
            526.0,
            528.0,
            540.0,
            543.0,
            548.0,
            549.0,
            553.0,
            555.0,
            565.0,
            572.0,
            574.0,
            579.0,
            580.0,
            585.0,
            586.0,
            588.0,
            593.0,
            594.0,
            597.0,
            601.0,
            603.0,
            605.0,
            608.0,
            617.0,
            618.0,
            620.0,
            621.0,
            624.0,
            626.0,
            633.0,
            634.0,
            646.0,
            663.0,
            666.0,
            677.0,
            700.0,
            705.0,
            713.0,
            714.0,
            716.0,
            722.0,
            724.0,
            726.0,
            727.0,
            728.0,
            744.0,
            746.0,
            748.0,
            757.0,
            759.0,
            775.0,
            777.0,
            782.0,
            787.0,
            792.0,
            794.0,
            797.0,
            803.0,
            806.0,
            811.0,
            817.0,
            819.0,
            822.0,
            825.0,
            841.0,
            843.0,
            845.0,
            848.0,
            849.0,
            850.0,
            863.0,
            865.0,
            866.0,
            867.0,
            871.0,
            878.0,
            884.0,
            899.0,
            916.0,
            918.0,
            919.0,
            920.0,
            923.0,
            925.0,
            926.0,
            927.0,
            930.0,
            932.0,
            933.0,
            934.0,
            935.0,
            936.0,
            940.0,
            945.0,
            947.0,
            952.0,
            954.0,
            960.0,
            963.0,
            967.0,
            968.0,
            973.0,
            974.0,
            979.0,
            987.0,
            988.0,
            990.0,
            991.0,
            998.0,
            1000.0,
        ]
        sub_meg_nocenter_ACG_KEY_fwd = [
            2907376.0,
            2907382.0,
            2907383.0,
            2907385.0,
            2907387.0,
            2907388.0,
            2907390.0,
            2907391.0,
            2907395.0,
            2907398.0,
            2907400.0,
            2907404.0,
            2907409.0,
            2907418.0,
            2907420.0,
            2907433.0,
            2907437.0,
            2907439.0,
            2907440.0,
            2907441.0,
            2907443.0,
            2907444.0,
            2907445.0,
            2907448.0,
            2907456.0,
            2907470.0,
        ]
        sub_wg_center_ACG_KEY_fwd = [
            10.0,
            83.0,
            106.0,
            125.0,
            403.0,
            442.0,
            451.0,
            479.0,
            553.0,
            565.0,
            782.0,
            810.0,
            835.0,
            841.0,
        ]
        sub_wg_nocenter_ACG_KEY_fwd = [2907383.0, 2907456.0]

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

        sub_meg_center_ACG_KEY_rev = [
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
        sub_meg_nocenter_ACG_KEY_rev = [
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
        sub_wg_center_ACG_KEY_rev = [
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
        sub_wg_nocenter_ACG_KEY_rev = [12, 18, 79, 32]

        sub_meg_center_ACG_KEY_fwd = [
            0,
            0,
            13,
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
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            99,
            0,
            5,
            0,
            196,
            0,
            0,
            16,
            0,
            0,
            0,
            0,
            0,
            0,
            86,
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
            70,
            0,
            0,
            0,
            0,
            22,
            0,
            0,
            29,
            0,
            0,
            16,
            0,
            0,
            0,
            0,
            0,
            0,
            215,
            131,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            88,
            0,
            2,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            104,
            0,
            230,
            0,
            167,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            26,
            0,
            11,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            80,
            0,
            37,
            0,
            0,
            1,
            27,
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
            211,
            186,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            101,
            0,
            47,
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
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            3,
            47,
            61,
            0,
            0,
            0,
            0,
            7,
            12,
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
        sub_meg_nocenter_ACG_KEY_fwd = [
            0,
            0,
            13,
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
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            99,
            0,
        ]
        sub_wg_center_ACG_KEY_fwd = [
            15,
            93,
            191,
            14,
            115,
            30,
            12,
            214,
            18,
            232,
            208,
            12,
            27,
            37,
        ]
        sub_wg_nocenter_ACG_KEY_fwd = [15, 93]

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


if __name__ == "__main__":
    unittest.main()
