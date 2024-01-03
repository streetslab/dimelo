import subprocess
def parse_bam_modkit(
    input_file,
    output_file,
):
    result = subprocess.run(['ls', '-l'], capture_output=True, text=True)
    return result