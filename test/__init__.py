import tempfile
from pathlib import Path
import subprocess
import gzip
import urllib

import pysam

ref_genome_url = 'https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.0.fasta.gz'

script_location = Path(__file__).resolve().parent
output_dir = script_location / 'output'
ref_genome_gz = output_dir / 'chm13.draft_v1.0.fasta.gz'
ref_genome_fasta = output_dir / 'chm13.draft_v1.0.fasta'

data_dir = script_location / 'data'
input_bam = data_dir / 'ctcf_demo.sorted.bam'

def download_reference(force_redownload = False):
    """
    download the reference genome to enable downstream operations
    """
    if not output_dir.exists():
        output_dir.mkdir()
    if ref_genome_fasta.exists() and not force_redownload:
        print("Reference genome already downloaded.")
        return ref_genome_fasta
    else:
        urllib.request.urlretrieve(ref_genome_url,ref_genome_gz)

        with gzip.open(ref_genome_gz,'rb') as gzip_file:
            with open(ref_genome_fasta,'wb') as output_file:
                for chunk in gzip_file:
                    output_file.write(chunk)

        ref_genome_gz.unlink()
        print("Reference genome downloaded and decompressed.")
        return ref_genome_fasta

def retag_bam(bam_path,force_retag = False):
    """
    use modkit to retag an input .bam file (to MM/ML vs mm/ml, with ?/. specified)
    by default only retags if the .updated file doesn't already exist
    """
    bam_updated_path = bam_path.parent.parent / 'output' / (bam_path.stem + '.updated' + bam_path.suffix)
    bam_updated_index_path = bam_path.parent.parent / 'output' / (bam_path.stem + '.updated' + bam_path.suffix + '.bai')
    if bam_updated_path.exists() and bam_updated_index_path.exists() and not force_retag:
        print('Input bam already retagged.')
        return bam_updated_path
    else:
        subprocess.run(['modkit','update-tags',str(bam_path),str(bam_updated_path),'--mode','ambiguous'])
        pysam.index(str(bam_updated_path))
        return bam_updated_path

class DiMeLoParsingTestCase:
    """
    This is the base class for any DiMeLo tests that need to parse data or create output files.
    """
    def setup_class(cls):
        cls._outDir = tempfile.TemporaryDirectory()
        cls.outDir = Path(cls._outDir.name)
        cls.reference_genome = download_reference()
        cls.megalodon_input = retag_bam(input_bam)

    def teardown_class(cls):
        cls._outDir.cleanup()

    def assertOutputFileExists(self, file_name: Path):
        """Fails test if the given file name is not found in the output directory"""
        file_path = self.outDir / file_name
        self.assertTrue(file_path.exists(), msg=f"{file_path} does not exist")