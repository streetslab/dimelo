import os
import unittest

import pandas as pd

from dimelo.functions import PrintDictionaryToTab, SaveMAPQHistogram
from dimelo.test import DiMeLoTestCase


class TestMAPQ(DiMeLoTestCase):
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


if __name__ == "__main__":
    unittest.main()
