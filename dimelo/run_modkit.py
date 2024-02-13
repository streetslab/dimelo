import subprocess
# I believe that pty does not currently work on Windows, although this may change in future releases: https://bugs.python.org/issue41663
# However, it may be that pywinpty, which is installable from pip, would work fine. That just needs to be tested with a Windows machine
# My current thinking is to wait on this until Nanopore puts Windows executables on Anaconda: https://anaconda.org/nanoporetech/modkit
import pty 
import re
import os
from pathlib import Path
import select

from tqdm.auto import tqdm

from . import utils

def run_with_progress_bars(
    command_list: list[str],
    input_file: Path,
    ref_genome: Path,
    motifs: list[str],
    load_fasta_regex: str,
    find_motifs_regex: str,
    contigs_progress_regex: str,
    single_contig_regex: str,
    buffer_size: int=50,
    progress_granularity: int=10,
    done_str: str='Done',
    err_str: str='Error',
    expect_done: bool=False,
    quiet: bool = False,
) -> str:
    """
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

    # Set up progress bar variables
    pbar_pre,pbar_contigs,pbar_chr = None,None,None
    format_pre = '{bar}| {desc} {percentage:3.0f}% | {elapsed}'
    format_contigs = '{bar}| {desc} {percentage:3.0f}% | {elapsed}<{remaining}'
    format_chr = '{bar}| {desc} {percentage:3.0f}%'
    finding_progress_dict = {}
    in_contig_progress = (0,1)   

    # Set up output buffer variables
    buffer_bytes = bytearray()
    tail_buffer = ''

    # Set up flags
    err_flag = False
    done_flag = False

    # Create a pseudo-terminal
    master_fd, slave_fd = pty.openpty()

    # Start subprocess with the slave end as stdio
    process = subprocess.Popen(command_list,stdin=slave_fd,stdout=slave_fd,stderr=subprocess.STDOUT,close_fds=True)
    os.close(slave_fd)

    if quiet:
        # We need to grab the outputs for the process to terminate but that's it, don't have to do anything with them in this case
        while True:
            ready, _, _ = select.select([master_fd], [], [], 0.1)
            if ready:
                try:
                    data = os.read(master_fd,1)
                    if not data:
                        break # No more data

                    # buffer_bytes += data  # Accumulate bytes in the buffer
                    # 
                    # try:
                    #     # Try to decode the accumulated bytes
                    #     text = buffer_bytes.decode('utf-8')
                    #     buffer_bytes.clear()  # Clear the buffer after successful decoding
                        # # If we have hit an error or modkit is done, just accumulate the rest of the output and then deal with it
                        # if err_flag or done_flag:
                        #     tail_buffer+=text
                        # # If we haven't hit an error or a done state, first check for that
                        # else :
                        #     tail_buffer = (tail_buffer + text)[-buffer_size:]
                        #     if err_str in tail_buffer:
                        #         index = tail_buffer.find(err_str)
                        #         tail_buffer = tail_buffer[index:]
                        #         err_flag = True
                        #     elif done_str in tail_buffer:
                        #         index = tail_buffer.find(done_str)
                        #         tail_buffer = tail_buffer[index-2:]
                        #         done_flag = True
                    # except:
                    #     continue
                except OSError:
                    break  # Handle errors - or just don't!
        process.wait()
        return ''
    # Grab output bytes as they come
    readout_count = 0
    while True:
        ready, _, _ = select.select([master_fd], [], [], 0.1)
        if ready:
            try:
                # Read a single byte
                data = os.read(master_fd, 1)
                if not data:
                    break  # No more data

                buffer_bytes += data  # Accumulate bytes in the buffer
                
                try:
                    # Try to decode the accumulated bytes
                    text = buffer_bytes.decode('utf-8')
                    readout_count += 1
                    buffer_bytes.clear()  # Clear the buffer after successful decoding
                    # If we have hit an error or modkit is done, just accumulate the rest of the output and then deal with it
                    if err_flag or done_flag:
                        tail_buffer+=text
                    # If we haven't hit an error or a done state, first check for that
                    else :
                        tail_buffer = (tail_buffer + text)[-buffer_size:]
                        if err_str in tail_buffer:
                            index = tail_buffer.find(err_str)
                            tail_buffer = tail_buffer[index:]
                            err_flag = True
                        elif done_str in tail_buffer:
                            index = tail_buffer.find(done_str)
                            tail_buffer = tail_buffer[index-2:]
                            done_flag = True
                        # If the process is ongoing, then go through the possible cases and create/adjust pbars accordingly
                        elif readout_count%progress_granularity==0:
                            # We check these in the reverse order from that in which they occur, which I guess will save a tiny
                            # amount of processing time because we don't check for previous steps when on later steps
                            if contigs_progress_matches := re.search(contigs_progress_regex,tail_buffer):
                                # Create and close pbars
                                if pbar_chr is None or pbar_contigs is None:
                                    if pbar_pre is not None:
                                        pbar_pre.n=100
                                        pbar_pre.set_description(f'Preprocessing complete for motifs {motifs} in {ref_genome.name}')
                                        pbar_pre.refresh()
                                        pbar_pre.close()
                                    pbar_contigs = tqdm(
                                        total=100,
                                        desc=f'Processing {input_file.name}',
                                        bar_format = format_contigs,
                                    )
                                    pbar_chr = tqdm(
                                        total=100,
                                        desc='',
                                        bar_format = format_chr,
                                    )
                                # This progress bar tracks how many contigs/chromosomes have been processed
                                current_contig = int(contigs_progress_matches.group(1))
                                total_contigs = int(contigs_progress_matches.group(2))
                                pbar_contigs.n = (100*(current_contig+
                                                       (in_contig_progress[0]/in_contig_progress[1] if in_contig_progress[1]>0 else 0)
                                                       ))/total_contigs
                                pbar_contigs.set_description(f'Processing contig {current_contig}/{total_contigs} from {input_file.name}')
                                pbar_contigs.refresh()

                            elif single_contig_matches := re.search(single_contig_regex,tail_buffer):
                                # Create and close pbars
                                if pbar_chr is None or pbar_contigs is None:
                                    if pbar_pre is not None:
                                        pbar_pre.n=100
                                        pbar_pre.set_description('Preprocessing complete for motifs {motifs} in {ref_genome.name}')
                                        pbar_pre.refresh()
                                        pbar_pre.close()
                                    pbar_contigs = tqdm(
                                        total=100,
                                        desc=f'Processing {input_file.name}',
                                        bar_format = format_contigs,
                                    )
                                    pbar_chr = tqdm(
                                        total=100,
                                        desc='',
                                        bar_format = format_chr,
                                    )
                                # This progress bar tracks reads processed within a chromosomes
                                chromosome = single_contig_matches.group(3)
                                reads_done = int(single_contig_matches.group(1))
                                reads_total = int(single_contig_matches.group(2))
                                pbar_chr.n = 100*reads_done/reads_total if reads_total>0 else 0
                                in_contig_progress = (reads_done,reads_total)
                                pbar_chr.set_description(f'{chromosome}: {reads_done}/{reads_total} reads processed')
                                pbar_chr.refresh()

                            elif find_motifs_matches := re.search(find_motifs_regex,tail_buffer):
                                if pbar_pre is not None:
                                    # pbar_pre = tqdm(
                                    #     total=100,
                                    #     desc='Preprocessing',
                                    #     bar_format=format_pre,
                                    # )
                                    finding_progress_dict[find_motifs_matches.group(3)] = (int(find_motifs_matches.group(1)),int(find_motifs_matches.group(2)))
                                    num_sum,denom_sum = 0,0
                                    for num,denom in finding_progress_dict.values():
                                        num_sum+=num
                                        denom_sum+=denom
                                    pbar_pre.n = 100*num_sum/denom_sum
                                    pbar_pre.set_description(f'Preprocessing: finding motif(s) {motifs}')
                                    pbar_pre.refresh()
                            elif load_fasta_match := re.search(load_fasta_regex,tail_buffer):
                                if pbar_pre is None:
                                    pbar_pre = tqdm(
                                        total=100,
                                        desc='Preprocessing',
                                        bar_format=format_pre,
                                    )
                                pbar_pre.n = 100*int(load_fasta_match.group(1))/24
                                pbar_pre.set_description(f'Preprocessing: reading {ref_genome.name}')
                                pbar_pre.refresh()
                except UnicodeDecodeError:
                    # If decoding fails, continue accumulating bytes
                    continue
                except Exception as e:
                    print(f'WARNING: unexpected error in progress tracking:\n{e}')

            except OSError:
                break  # Handle errors - or just don't!
    process.wait()
    if err_flag:
        raise ValueError('modkit raised the following error:\n'+tail_buffer)
    elif done_flag or not expect_done:
        try:
            pbar_contigs.n=100
            pbar_contigs.set_description(f'All regions complete in {input_file.name}')
            pbar_contigs.refresh()
            pbar_contigs.close()
            pbar_chr.n=100
            ansi_escape_pattern = re.compile(r'(\[2K>)')
            pbar_chr.set_description(ansi_escape_pattern.sub('',tail_buffer).strip())
            pbar_chr.refresh()
            pbar_chr.close()
        except Exception as e:
            print(f'WARNING: unexpected error in progress tracking:\n{e}')
        return tail_buffer
    else:
        try:
            pbar_contigs.close()
            pbar_chr.close()
        except Exception as e:
            print(f'WARNING: unexpected error in progress tracking:\n{e}')
        print('WARNING: the modkit command may not have completed normally. Consider re-running with "log=True" if you do not get the expected outputs.')
        return tail_buffer