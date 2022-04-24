import tempfile
import unittest
from pathlib import Path


class DiMeLoTestCase(unittest.TestCase):
    def setUp(self):
        self._outDir = tempfile.TemporaryDirectory()
        self.outDir = Path(self._outDir.name)

    def tearDown(self):
        self._outDir.cleanup()
    
    # def tmpFile(self):
    #     tempFile = tempfile.NamedTemporaryFile(delete=True)
    #     tempFile.close()
    #     return tempFile.name
