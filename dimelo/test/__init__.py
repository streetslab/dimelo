import tempfile
import unittest
from pathlib import Path
import shutil


class DiMeLoTestCase(unittest.TestCase):
    """
    TODO:
        - Should these be setUpClass/tearDownClass or setUp/tearDown? Is it okay that only one temporary directory is created each time?
    """

    @classmethod
    def setUpClass(cls):
        # cls._outDir = tempfile.TemporaryDirectory()
        # cls.outDir = Path(cls._outDir.name)
        cls.outDir = Path('dimelo/dimelo_test')

    @classmethod
    def tearDownClass(cls):
        # cls._outDir.cleanup()
        shutil.rmtree(cls.outDir)

    # def tmpFile(self):
    #     tempFile = tempfile.NamedTemporaryFile(delete=True)
    #     tempFile.close()
    #     return tempFile.name
