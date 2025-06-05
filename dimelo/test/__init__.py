import gzip
import subprocess
import tempfile
import urllib
from inspect import signature
from pathlib import Path

import pysam

ref_genome_url = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.0.fasta.gz"

script_location = Path(__file__).resolve().parent
output_dir = script_location / "output"
ref_genome_gz = output_dir / "chm13.draft_v1.0.fasta.gz"
ref_genome_fasta = output_dir / "chm13.draft_v1.0.fasta"

data_dir = script_location / "data"
bams_to_update = [data_dir / "ctcf_demo.sorted.bam"]

"""
The RelativePath class is a pickle-able Path class that will play nice with parse_bam.sanitize_path_inputs (using the __fspath__ function)
and then can also be included in dicts and other objects for pickling. This way we can store input and output test files using relative
paths but execute code using e.g. pytest, which operates in a different directory, and have the paths work by converting them to absolute
paths. Clumsy maybe? Could be a better approach out there? Who knows. Not I.
"""


class RelativePath:
    def __init__(self, path):
        path = Path(path)
        if path.is_absolute():
            # Convert the absolute path to a relative path
            try:
                self.relative_path = path.relative_to(script_location)
            except ValueError as e:
                raise ValueError(
                    "The provided path is not in the dimelo/test directory."
                ) from e
        else:
            self.relative_path = path
        self._update_absolute_path()

    def _update_absolute_path(self):
        # Dynamically determine the absolute path based on relative path and the location of this file
        self.base_path = Path(__file__).parent
        self.absolute_path = (self.base_path / self.relative_path).resolve()

    def __fspath__(self):
        # Allows the object to be used by pathlib.Path and anything that accepts path-like objects
        return str(self.absolute_path)

    def __str__(self):
        # Allows the object to be used by methods that want strings
        return str(self.absolute_path)

    def __getstate__(self):
        # Returns the state to be pickled; i.e. just the relative path
        return self.relative_path

    def __setstate__(self, state):
        # Restore state from the unpickled state; recreate the absolute path from the relative path
        self.relative_path = state
        self._update_absolute_path()


def filter_kwargs_for_func(func, kwargs, extra_args=[]):
    func_sig = signature(func)
    allowed_args = set(list(func_sig.parameters) + extra_args)
    filtered_kwargs = {k: v for k, v in kwargs.items() if k in allowed_args}
    return filtered_kwargs


def download_reference(force_redownload=False):
    """
    download the reference genome to enable downstream operations
    """
    if not output_dir.exists():
        output_dir.mkdir()
    if ref_genome_fasta.exists() and not force_redownload:
        print("Reference genome already downloaded.")
        return ref_genome_fasta
    else:
        urllib.request.urlretrieve(ref_genome_url, ref_genome_gz)

        with (
            gzip.open(ref_genome_gz, "rb") as gzip_file,
            open(ref_genome_fasta, "wb") as output_file,
        ):
            for chunk in gzip_file:
                output_file.write(chunk)

        ref_genome_gz.unlink()
        print("Reference genome downloaded and decompressed.")
        return ref_genome_fasta


def retag_bam(bam_path, force_retag=False):
    """
    use modkit to retag an input .bam file (to MM/ML vs mm/ml, with ?/. specified)
    by default only retags if the .updated file doesn't already exist
    """
    bam_updated_path = (
        bam_path.parent.parent
        / "output"
        / (Path(bam_path.stem).stem + ".updated" + bam_path.suffix)
    )
    bam_updated_index_path = (
        bam_path.parent.parent
        / "output"
        / (Path(bam_path.stem).stem + ".updated" + bam_path.suffix + ".bai")
    )
    if (
        bam_updated_path.exists()
        and bam_updated_index_path.exists()
        and not force_retag
    ):
        print("Input bam already retagged.")
        return bam_updated_path
    else:
        subprocess.run(
            [
                "modkit",
                "update-tags",
                str(bam_path),
                str(bam_updated_path),
                "--mode",
                "ambiguous",
            ]
        )
        pysam.index(str(bam_updated_path))
        return bam_updated_path


class DiMeLoParsingTestCase:
    """
    This is the base class for any DiMeLo tests that need to parse data or create output files.
    """

    @classmethod
    def setup_class(cls):
        # TODO: The existence of two variables with basically the same name is very confusing. Please clarify?
        cls._outDir = tempfile.TemporaryDirectory()
        cls.outDir = Path(cls._outDir.name)
        cls.reference_genome = download_reference()
        _ = [retag_bam(bam) for bam in bams_to_update]

    @classmethod
    def teardown_class(cls):
        cls._outDir.cleanup()

    def assertOutputFileExists(self, file_name: Path):
        """Fails test if the given file name is not found in the output directory"""
        """
        TODO: There are a couple of things wrong here:
        * This is never called anywhere; was this a holdover from my old infrastructure? Can it be deleted?
        * mypy error: dimelo/test/__init__.py:157: error: "DiMeLoParsingTestCase" has no attribute "outDir"  [attr-defined]
        """
        file_path = self.outDir / file_name
        assert file_path.exists(), f"{file_path} does not exist"
