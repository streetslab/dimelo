import os

# I believe that pty does not currently work on Windows, although this may change in future releases: https://bugs.python.org/issue41663
# However, it may be that pywinpty, which is installable from pip, would work fine. That just needs to be tested with a Windows machine
# My current thinking is to wait on this until Nanopore puts Windows executables on Anaconda: https://anaconda.org/nanoporetech/modkit
import pty
import re
import select
import subprocess
import sys
from pathlib import Path
from typing import Optional, cast

from tqdm.auto import tqdm

# This should be updated in tandem with the environment.yml nanoporetech::modkit version
EXPECTED_MODKIT_VERSION = "0.2.4"

def ensure_modkit_available(quiet: bool = False) -> None:
    """
    Lazily check that modkit is on PATH and matches expected version.
    Called by parse functions to avoid import-time failures during analysis-only workflows.
    """

    # Add conda env bin folder to path if it is not already present
    current_interpreter = sys.executable
    env_bin_path = os.path.dirname(current_interpreter)
    if env_bin_path not in os.environ["PATH"]:
        if not quiet:
            print(
                f"PATH does not include the conda environment /bin folder. Adding {env_bin_path}."
            )
        os.environ["PATH"] = f"{env_bin_path}:{os.environ['PATH']}"

    try:
        result = subprocess.run(
            ["modkit", "--version"], stdout=subprocess.PIPE, text=True, check=True
        )
        modkit_version = result.stdout.split()
        if len(modkit_version) > 1 and modkit_version[1] != EXPECTED_MODKIT_VERSION:
            if not quiet:
                print(
                    f"modkit found with unexpected version {modkit_version[1]}. "
                    f"Expected {EXPECTED_MODKIT_VERSION}."
                )
    except Exception:
        raise FileNotFoundError(
            'Executable not found for modkit. Install dimelo using "conda env create -f environment.yml" '
            'or install modkit manually using "conda install nanoporetech::modkit==0.2.4". '
            "Without modkit you cannot run parse_bam functions."
        )


def run_with_progress_bars(
    command_list: list[str],
    input_file: Path,
    ref_genome: Path,
    motifs: list[str],
    load_fasta_regex: str,
    find_motifs_regex: str,
    contigs_progress_regex: str,
    single_contig_regex: str,
    buffer_size: int = 50,
    progress_granularity: int = 10,
    done_str: str = "Done",
    err_str: str = "Error",
    expect_done: bool = False,
    quiet: bool = False,
) -> str:
    r"""
    This function runs modkit with subprocess / pseudoterminal and grabs the progress outputs to populate progress bars

    Args:
        command_list: a list of commands to pass to subprocess: [modkit, pileup, ...] or [modkit, extract, ...]
        load_fasta_regex: a regular expression that captures the contig being loaded in the step where modkit
            reads fasta sequence. Should specify all the output context
            so that groups aren't captured unless there is whitespace on either end i.e. the whole output
            e.g. r'\s+\[.*?\]\s+(\d+)\s+Reading' for pileup in 0.2.4
        input_file: the bam file you are processing
        ref_genome: the reference genome to which your bam is aligned
        motifs: the list of motifs you are looking for
        find_motifs_regex: a regular expression that captures contigs-so-far and total-contigs-to-process
            in the step where modkit is finding motifs throughout the genome. Should specify all the output context
            so that groups aren't captured unless there is whitespace on either end i.e. the whole output
            has been loaded into the buffer
            e.g. r'\s+(\d+)/(\d+)\s+finding\s+([A-Za-z0-9,]+)\s+motifs' for pileup in 0.2.4
        contigs_progress_regex: a regular expression that captures currently-processing-contig and total-contigs
            -to-process in the step where modkit is running through the bam file. Should specify all the output context
            so that groups aren't captured unless there is whitespace on either end i.e. the whole output
            e.g. r'\s+(\d+)/(\d+)\s+contigs' for pileup in 0.2.4
        single_contig_regex: a regular expression that captures reads-processed, reads-total, and contig-name for
            a contig that is being processing from the bam file. Should specify all the output context
            so that groups aren't captured unless there is whitespace on either end i.e. the whole output
            e.g. r'\s+(\d+)/(\d+)\s+processing\s+([\w]+)[^\w]' for pileup in 0.2.4
        buffer_size: the length of the string that the modkit stderr output gets saved into. This size will not
            be respected if you hit Done or Error; in that case the rest of the output will be captured and returned
            or raised.
        progress_granularity: this tells the function how often to check the output buffer string for the various regex.
            Less frequent checking is good because it means fewer spurious updates and less overhead. However you
            need this to be sufficiently less than buffer_size that you can always capture the entirety of your
            relevant information.
        done_str: a string telling the function what to look for to know that modkit is done processing. Everything
            after this will get returned
        err_str: a string telling the function what to look for to know that modkit has encountered an error. Everything
            after this will be raised as a ValueError
        expect_done: specifies whether the command is expected to show a clear "Done" at the end of the output
        quiet: sending True will suppress all progress bars and stdout outputs.

    Returns:
        The command line stderr output string after the point where we detect modkit is done parsing
    """

    # modkit 0.2.4 does not like gzipped or bgzipped fasta files
    with open(ref_genome, "rb") as f:
        if f.read(2) == b"\x1f\x8b":
            raise ValueError(
                f"{ref_genome.name} is gzipped, which will cause modkit to fail.\ngunzip {ref_genome.name} and try again."
            )

    # Set up progress bar variables to display progress updates when not in quiet mode
    format_pre = "{bar}| {desc} {percentage:3.0f}% | {elapsed}"
    format_contigs = "{bar}| {desc} {percentage:3.0f}% | {elapsed}<{remaining}"
    format_chr = "{bar}| {desc} {percentage:3.0f}%"
    pbar_pre: Optional[tqdm] = None
    pbar_contigs: Optional[tqdm] = None
    pbar_chr: Optional[tqdm] = None

    # TODO: Is this the correct type annotation? I think it is, based on approx. line 280
    finding_progress_dict: dict[str, tuple[int, int]] = {}
    in_contig_progress = (0, 1)
    total_contigs = 0

    # Set up output buffer variables to capture modkit output
    buffer_bytes = bytearray()
    tail_buffer = ""

    # Set up flags for modkit error / modkit done from text outputs
    err_flag = False
    done_flag = False

    # Create a pseudo-terminal in which to run the modkit subprocess
    master_fd, slave_fd = pty.openpty()

    # Start modkit subprocess with the slave end as stdio
    process = subprocess.Popen(
        command_list,
        stdin=slave_fd,
        stdout=slave_fd,
        stderr=subprocess.STDOUT,
        close_fds=True,
    )
    os.close(slave_fd)

    readout_count = 0
    progress_bars_initialized = False
    region_parsing_started = False

    # Grab output bytes for as long as they're coming
    while True:
        # Wait for the process to be ready to provide bytes
        ready, _, _ = select.select([master_fd], [], [], 0.1)
        if ready:
            try:
                # Read a single byte
                data = os.read(master_fd, 1)
                if not data:
                    break  # No more data

                if quiet:
                    # If we are in quiet mode, nothing gets grabbed
                    continue
                else:
                    # Create the progress bars when first entering this code block
                    if not progress_bars_initialized:
                        pbar_pre = tqdm(
                            total=100,
                            desc=f"Step 1: Identify motif locations in {ref_genome.name}",
                            bar_format=format_pre,
                        )
                        pbar_contigs = tqdm(
                            total=100,
                            desc=f"Step 2: Parse regions in {input_file.name}",
                            bar_format=format_contigs,
                        )
                        pbar_chr = tqdm(
                            total=100,
                            desc="",
                            bar_format=format_chr,
                        )
                        progress_bars_initialized = True

                    buffer_bytes += data  # Accumulate bytes in the buffer

                    try:
                        # Try to decode the accumulated bytes
                        # This will throw a UnicodeDecodeError if not complete, which is ok! Then we just continue on
                        text = buffer_bytes.decode("utf-8")
                        readout_count += 1
                        buffer_bytes.clear()  # Clear the buffer after successful decoding
                        # If we have hit an error or modkit is done, just accumulate the rest of the output and then deal with it:
                        # no need to check the progress tracking stuff in that case
                        if err_flag or done_flag:
                            tail_buffer += text
                        # If we haven't hit an error or a done state, first check for that
                        else:
                            tail_buffer = (tail_buffer + text)[-buffer_size:]
                            if err_str in tail_buffer:
                                index = tail_buffer.find(err_str)
                                tail_buffer = tail_buffer[index:]
                                err_flag = True
                            elif done_str in tail_buffer:
                                index = tail_buffer.find(done_str)
                                tail_buffer = tail_buffer[index - 2 :]
                                done_flag = True
                            # If the process is ongoing, then go through the possible cases and create/adjust pbars accordingly
                            # We only sometimes want to update progress because otherwise the constant updates slow us down
                            elif (
                                readout_count % progress_granularity == 0
                                and progress_bars_initialized
                            ):
                                region_parsing_started, in_contig_progress = (
                                    update_progress_bars(
                                        pbar_pre=pbar_pre,
                                        pbar_contigs=pbar_contigs,
                                        pbar_chr=pbar_chr,
                                        tail_buffer=tail_buffer,
                                        contigs_progress_regex=contigs_progress_regex,
                                        single_contig_regex=single_contig_regex,
                                        find_motifs_regex=find_motifs_regex,
                                        load_fasta_regex=load_fasta_regex,
                                        region_parsing_started=region_parsing_started,
                                        in_contig_progress=in_contig_progress,
                                        finding_progress_dict=finding_progress_dict,
                                        ref_genome=ref_genome,
                                        input_file=input_file,
                                        motifs=motifs,
                                    )
                                )

                    except UnicodeDecodeError:
                        # If decoding fails, continue accumulating bytes
                        continue
                    except Exception as e:
                        raise e
            except OSError:
                break

    # After the data stops coming, we wait until the process is done so we can grab the return code
    process.wait()
    return_code = process.returncode
    # Modkit gives return code 0 if it terminates successfully; any other return code should be raised
    # This catches system kills caused by memory and disk space
    if return_code != 0 or err_flag:
        if err_flag:
            print(tail_buffer)
        if return_code == -9:
            print(
                "It looks like the process was killed with SIGKILL. This is often due to the system running out of memory.",
                "\nConsider setting 'cores=1' to reduce resource usage or reducing the data size being processed.",
                "\nNote that pileup and extract both create large intermediate plain-text files before reducing down to a compressed and indexed format.",
                "\nRunning out of disk space can indirectly cause memory issues."
                "\nYou might also want to check if there are any resource limits set for your user or in the environment where the code is running.",
            )
        raise subprocess.CalledProcessError(
            return_code, command_list, output=tail_buffer
        )

    # If modkit gives return code 0, it can still have had an unusual state.
    # If the expected done flag was seen, or if the command doesn't have a done
    # string (i.e. extract), update the progress bars to reflect the final status
    # If the progress bars are not initialized, then the code must have been run in
    # quiet mode
    elif done_flag or not expect_done:
        if progress_bars_initialized:
            pbar_pre = cast(tqdm, pbar_pre)
            pbar_contigs = cast(tqdm, pbar_contigs)
            pbar_chr = cast(tqdm, pbar_chr)
            pbar_pre.close()
            pbar_contigs.n = 100
            pbar_contigs.set_description(
                f"Step 2 complete. {total_contigs} contigs processed from {input_file.name}"
            )
            pbar_contigs.refresh()
            pbar_contigs.close()
            pbar_chr.n = 100
            ansi_escape_pattern = re.compile(r"(\[2K>)")
            pbar_chr.set_description(
                command_list[0]
                + " "
                + command_list[1]
                + " return code "
                + str(return_code)
                + " | "
                + ansi_escape_pattern.sub("", tail_buffer).strip()
            )
            pbar_chr.refresh()
            pbar_chr.close()
        return tail_buffer
    # Indicate unusual state
    else:
        if progress_bars_initialized:
            pbar_pre = cast(tqdm, pbar_pre)
            pbar_contigs = cast(tqdm, pbar_contigs)
            pbar_chr = cast(tqdm, pbar_chr)
            pbar_pre.close()
            pbar_contigs.set_description("Unexpected modkit outputs")
            pbar_contigs.refresh()
            pbar_contigs.close()
            pbar_chr.set_description(
                command_list[0]
                + " "
                + command_list[1]
                + " return code "
                + str(return_code)
            )
            pbar_chr.refresh()
            pbar_chr.close()
        print(
            'WARNING: the modkit command may not have completed normally. Consider re-running with "log=True" if you do not get the expected outputs.'
        )
        return tail_buffer


def update_progress_bars(
    pbar_pre,
    pbar_contigs,
    pbar_chr,
    tail_buffer,
    contigs_progress_regex,
    single_contig_regex,
    find_motifs_regex,
    load_fasta_regex,
    region_parsing_started,
    in_contig_progress,
    finding_progress_dict,
    ref_genome,
    input_file,
    motifs,
):
    # We check these in the reverse order from that in which they occur, which I guess will save a tiny
    # amount of processing time because we don't check for previous steps when on later steps
    # Once we are in the contig progress stage, step 1 is done by definition
    if contigs_progress_matches := re.search(contigs_progress_regex, tail_buffer):
        # If we get here we can be sure the pbars are initialized
        pbar_pre = cast(tqdm, pbar_pre)
        pbar_contigs = cast(tqdm, pbar_contigs)
        if not region_parsing_started:
            # These are now no longer indicating future steps, but rather counting the actual
            # time for step 2
            pbar_contigs.reset()
            pbar_chr.reset()
            region_parsing_started = True
            # Now that region parsing has started, we can close out the preprocessing pbar
            pbar_pre.n = 100
            pbar_pre.set_description(
                f"Step 1 complete. Located {motifs} in {ref_genome.name}"
            )
            pbar_pre.refresh()
            pbar_pre.close()
        # This progress bar tracks how many contigs/chromosomes have been processed
        current_contig = int(contigs_progress_matches.group(1))
        total_contigs = int(contigs_progress_matches.group(2))
        pbar_contigs.n = (
            (
                100
                * (
                    current_contig
                    + (
                        in_contig_progress[0] / in_contig_progress[1]
                        if in_contig_progress[1] > 0
                        else 0
                    )
                )
            )
            / total_contigs
            if total_contigs > 0
            else 0
        )
        pbar_contigs.set_description(
            f"Step 2: parsing {current_contig}/{total_contigs} from {input_file.name}"
        )
        pbar_contigs.refresh()
    elif region_parsing_started and (
        single_contig_matches := re.search(single_contig_regex, tail_buffer)
    ):
        # If we get here we can be sure the pbars are initialized
        pbar_chr = cast(tqdm, pbar_chr)
        # This progress bar tracks reads processed within a chromosomes
        chromosome = single_contig_matches.group(3)
        reads_done = int(single_contig_matches.group(1))
        reads_total = int(single_contig_matches.group(2))
        pbar_chr.n = 100 * reads_done / reads_total if reads_total > 0 else 0
        in_contig_progress = (reads_done, reads_total)
        pbar_chr.set_description(
            f"Step 2: {chromosome} {reads_done}/{reads_total} chunks processed"
        )
        pbar_chr.refresh()

    elif find_motifs_matches := re.search(find_motifs_regex, tail_buffer):
        # If we get here we can be sure the pbars are initialized
        pbar_pre = cast(tqdm, pbar_pre)
        finding_progress_dict[find_motifs_matches.group(3)] = (
            int(find_motifs_matches.group(1)),
            int(find_motifs_matches.group(2)),
        )
        num_sum, denom_sum = 0, 0
        for (
            num,
            denom,
        ) in finding_progress_dict.values():
            num_sum += num
            denom_sum += denom
        if denom_sum > 0:
            pbar_pre.n = 100 * num_sum / denom_sum
        else:
            pbar_pre.n = 0
        pbar_pre.set_description(f"Step 1b: finding motif(s) {motifs}")
        pbar_pre.refresh()
    elif load_fasta_match := re.search(load_fasta_regex, tail_buffer):
        # If we get here we can be sure the pbars are initialized
        pbar_pre = cast(tqdm, pbar_pre)
        pbar_pre.n = 100 * int(load_fasta_match.group(1)) / 24
        pbar_pre.set_description(f"Step 1a: reading {ref_genome.name}")
        pbar_pre.refresh()
    return region_parsing_started, in_contig_progress
