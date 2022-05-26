import tempfile
import unittest
from pathlib import Path


class DiMeLoTestCase(unittest.TestCase):
    """
    TODO:
        - Should these be setUpClass/tearDownClass or setUp/tearDown? Is it okay that only one temporary directory is created each time?
    """

    @classmethod
    def setUpClass(cls):
        cls._outDir = tempfile.TemporaryDirectory()
        cls.outDir = Path(cls._outDir.name)

    @classmethod
    def tearDownClass(cls):
        cls._outDir.cleanup()

    def assertOutputFileExists(self, file_name: Path):
        """Fails test if the given file name is not found in the output directory"""
        file_path = self.outDir / file_name
        self.assertTrue(file_path.exists(), msg=f"{file_path} does not exist")

    # def tmpFile(self):
    #     tempFile = tempfile.NamedTemporaryFile(delete=True)
    #     tempFile.close()
    #     return tempFile.name
