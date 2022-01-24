"Functions to parse_bam for QC"

import multiprocessing
import sqlite3
import time
from math import log

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam

# import seaborn as sns
from joblib import Parallel, delayed

from dimelo.parse_bam import get_modified_reference_positions, make_db
from dimelo.utils import execute_sql_command

# def create_sql_table(database_name, table_name, cols, d_types):
#     if os.path.exists(database_name):
#         os.remove(database_name)
#     if os.path.exists(database_name + '-journal'):
#         os.remove(database_name + '-journal')
#
#     clear_db(database_name)
#
#     conn = sqlite3.connect(database_name)
#     c = conn.cursor()
#
#     s = ''
#     for i in range(len(cols)):
#         if i == 0:
#             s = s + cols[i] + ' ' + d_types[i] + ' ' + 'primary key, '
#         elif i == len(cols) - 1:
#             s = s + cols[i] + ' ' + d_types[i]
#         else:
#             s = s + cols[i] + ' ' + d_types[i] + ',' + ' '
#
#     fs = '(' + s + ')'
#     c.execute('''create table if not exists ''' + table_name + ''' ''' + fs + ''';''')
#     conn.commit()
#
#
# def execute_sql_command(command: str, database_name: str, values) -> none:
#     """
#     function to execute a sql command from python.
#
#     parameters
#     ----------
#     command: str
#         sql command (use strings with three quotes on each side
#         so that it can be a multiline string
#     database_name: str
#         file name of the database (e.g, "my.db")
#     values: list or tuple
#         data per row
#
#     returns
#     -------
#     no return, executes the command
#     """
# """
#     # will create if not present
#     conn = sqlite3.connect(database_name, timeout=30)
#     c = conn.cursor()
#
#     if len(values) == 0:
#         c.execute(command)
#     elif type(values) == list:
#         c.executemany(command, values)
#     else:
#         c.execute(command, values)
#     # saves the changes
#     conn.commit()


# def get_modified_reference_positions(read, basemod, threshA, threshC):
#     """Extract mA and mC pos & prob information for the read
#     Args:
#             :param read: single read from bam file
#             :param basemod: which basemods, currently supported options are 'A', 'CG', 'A+CG'
#             :param window: window from bed file
#             :param center: report positions with respect to reference center (+/- window size) if
#             True or in original reference space if False
#             :param threshA: threshold above which to call an A base methylated
#             :param threshC: threshold above which to call a C base methylated
#             :param windowSize: window size around center point of feature of interest to plot (+/-);
#             only mods within this window are stored; only applicable for center=True
#     """
#     if (read.has_tag("Mm")) & (";" in read.get_tag("Mm")):
#         mod1 = read.get_tag("Mm").split(";")[0].split(",", 1)[0]
#         mod2 = read.get_tag("Mm").split(";")[1].split(",", 1)[0]
#         mod1_list = read.get_tag("Mm").split(";")[0].split(",", 1)
#         mod2_list = read.get_tag("Mm").split(";")[1].split(",", 1)
#         base = basemod[0]  # this will be A, C, or A
#         if basemod == "A+CG":
#             base2 = basemod[2]  # this will be C for A+C case
#         else:  # in the case of a single mod will just be checking that single base
#             base2 = base
#         if len(mod1_list) > 1 and (base in mod1 or base2 in mod1):
#             mod1_return = get_mod_reference_positions_by_mod(
#                 read,
#                 mod1,
#                 0,
#                 threshA,
#                 threshC,
#
#             )
#         else:
#             mod1_return = (None, [None], [None])
#         if len(mod2_list) > 1 and (base in mod2 or base2 in mod2):
#             mod2_return = get_mod_reference_positions_by_mod(
#                 read,
#                 mod2,
#                 1,
#                 threshA,
#                 threshC,
#             )
#             return mod1_return, mod2_return
#         else:
#             return mod1_return, (None, [None], [None])
#     else:
#         return (None, [None], [None]), (None, [None], [None])
#
#
# def get_mod_reference_positions_by_mod(read,
#                                        basemod,
#                                        index,
#                                        threshA,
#                                        threshC,
#                                        ):
#     """Get positions and probabilities of modified bases for a single read
#     Args:
#             :param read: one read in bam file
#             :param mod: which basemod, reported as base+x/y/m
#             :param window: window from bed file
#             :param center: report positions with respect to reference center (+/- window size)
#             if True or in original reference space if False
#             :param threshA: threshold above which to call an A base methylated
#             :param threshC: threshold above which to call a C base methylated
#             :param windowSize: window size around center point of feature of interest to plot (+/-);
#             only mods within this window are stored; only applicable for center=True
#     """
#     base, mod = basemod.split("+")
#     deltas = [
#         int(i) for i in read.get_tag("Mm").split(";")[index].split(",")[1:]
#     ]
#     num_base = len(read.get_tag("Mm").split(";")[index].split(",")) - 1
#     Ml = read.get_tag("Ml")
#     if index == 0:
#         probabilities = np.array(Ml[0:num_base], dtype=int)
#     if index == 1:
#         probabilities = np.array(Ml[0 - num_base:], dtype=int)
#     base_index = np.array(
#         [
#             i
#             for i, letter in enumerate(read.get_forward_sequence())
#             if letter == base
#         ]
#     )
#     # determine locations of the modified bases, where index_adj is the adjustment of the base_index
#     # based on the cumulative sum of the deltas
#     locations = np.cumsum(deltas)
#     # loop through locations and increment index_adj by the difference between the next location and current one + 1
#     # if the difference is zero, therefore, the index adjustment just gets incremented
#     by one because no base should be skipped
#     index_adj = []
#     index_adj.append(locations[0])
#     i = 0
#     for i in range(len(locations) - 1):
#         diff = locations[i + 1] - locations[i]
#         index_adj.append(index_adj[i] + diff + 1)
#     # get the indices of the modified bases
#     modified_bases = base_index[index_adj]
#     refpos = np.array(read.get_reference_positions(full_length=True))
#     if read.is_reverse:
#         refpos = np.flipud(refpos)
#         probabilities = probabilities[::-1]
#
#     # extract CpG sites only rather than all mC
#     keep = []  # ind 3
#     prob_keep = []
#     all_bases_index = []
#     i = 0
#     seq = read.get_forward_sequence()
#     # deal with None for refpos from soft clipped / unaligned bases
#     if "C" in basemod:
#         for b in base_index:
#             if (
#                     b < len(seq) - 1
#             ):  # if modified C is not the last base in the read
#                 if (refpos[b] is not None) & (refpos[b + 1] is not None):
#                     if seq[b + 1] == "G":
#                         if (
#                                 abs(refpos[b + 1] - refpos[b]) == 1
#                         ):  # ensure there isn't a gap
#                             all_bases_index.append(
#                                 b
#                             )  # add to all_bases_index whether or not modified
#                             if b in modified_bases:
#                                 if probabilities[i] >= threshC:
#                                     keep.append(b)  # ind 3
#                                     prob_keep.append(i)
#             # increment for each instance of modified base
#             if b in modified_bases:
#                 i = i + 1
#     else:  # for m6A no need to look at neighboring base; do need to remove refpos that are None
#         for b in base_index:
#             if refpos[b] is not None:
#                 all_bases_index.append(
#                     b
#                 )  # add to all_bases_index whether or not modified
#                 if b in modified_bases:
#                     if probabilities[i] >= threshA:
#                         keep.append(b)  # ind 3
#                         prob_keep.append(i)
#             # increment for each instance of modified base
#             if b in modified_bases:
#                 i = i + 1
#
#     return basemod, len(all_bases_index), probabilities[prob_keep]


def batch_read_generator(file_bamIn, batch_size, outDir):
    counter = 0
    r_list = []

    for read in file_bamIn.fetch(until_eof=True):
        isA = False
        isC = False
        r = [
            read.query_name,
            read.reference_name,
            read.reference_start,
            read.reference_end,
            read.query_length,
            "-" if read.is_reverse else "+",
            read.mapping_quality,
            ave_qual(read.query_qualities),
            ave_qual(read.query_alignment_qualities),
        ]
        [
            (mod, numA, probs),
            (mod2, numC, probs2),
        ] = get_modified_reference_positions(
            read,
            "A+CG",
            None,
            False,
            129,
            129,
            None,
            None,
            None,
            outDir,
            False,
            True,
        )

        if type(numA) == int:
            isA = True
            problist = [int(x) / 255 for x in probs]
            probagg_05 = prob_bin(problist)

        if type(numC) == int:
            isC = True
            prob2list = [int(x) / 255 for x in probs2]
            prob2agg_05 = prob_bin(prob2list)

        if isC and not isA:
            r.extend((np.nan, np.nan, numC, prob2agg_05))
        elif isA and not isC:
            r.extend((numA, probagg_05, np.nan, np.nan))
        elif isA and isC:
            r.extend((numA, probagg_05, numC, prob2agg_05))

        r = tuple(r)
        if counter < batch_size:
            r_list.append(r)
            counter += 1
        else:
            yield r_list
            counter = 0
            r_list = [r]
    yield r_list


def prob_bin(bin):
    # probability a base in the window (or across reads or across bases within a read) is methylated by:
    # calculating probability that no base in the window (or across reads) is methylated and then taking the complement
    # treat p=1 as 254/255 for prevent log(0)
    # print(bin)
    probs = [
        np.log(1 - p) for p in bin if ((p < 1) and (p >= 0.5))
    ]  # only consider probabilities > 0.5 and handle 1 on next line
    probs1 = [np.log(1 - 254 / 255) for p in bin if p == 1]
    probsAll = probs + probs1
    prob = 1 - np.exp(sum(probsAll))
    return prob


def errs_tab(n):
    """Generate list of error rates for qualities less than equal than n."""
    return [10 ** (q / -10) for q in range(n + 1)]


def ave_qual(quals, qround=False, tab=errs_tab(129)):
    """Calculate average basecall quality of a read.
    Receive the integer quality scores of a read and return the average quality for that read
    First convert Phred scores to probabilities,
    calculate average error probability
    convert average back to Phred scale
    """
    if quals:
        mq = -10 * log(sum([tab[q] for q in quals]) / len(quals), 10)
        if qround:
            return round(mq)
        else:
            return mq
    else:
        return None


def parse_bam_read(bamIn, outDir):
    file_bamIn = pysam.AlignmentFile(bamIn, "rb")

    # DB_NAME = bamIn.split('/')[-1].split('.')[0] + ".db"
    # cols = ['name', 'chr', 'start', 'end', 'length', 'strand', 'mapq', 'ave_baseq', 'ave_alignq', 'basemod',
    #        'numA', 'methAprob05', 'numC', 'methCprob05']
    # dt = ['TEXT', 'TEXT', 'INT', 'INT', 'INT', 'TEXT', 'INT', 'INT', 'INT',
    # 'TEXT', 'FLOAT', 'FLOAT', 'FLOAT', 'FLOAT']
    # table_name = 'READS'
    # create_sql_table(DB_NAME, table_name, cols, dt)

    DB_NAME, tables = make_db(bamIn, "", outDir, True, True)
    template_command = (
        """INSERT INTO """
        + tables[0]
        + """ VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?);"""
    )
    num_cores = multiprocessing.cpu_count() - 1

    Parallel(n_jobs=num_cores, verbose=10)(
        delayed(execute_sql_command)(template_command, DB_NAME, i)
        for i in batch_read_generator(file_bamIn, 100, outDir)
    )
    return DB_NAME, tables[0]


def get_runtime(f, inp1, inp2):
    start = time.time()
    re_val = f(inp1, inp2)
    time.sleep(1)
    end = time.time()
    return f"Runtime of the program is {end - start}", re_val


def qc_report(filebamIn):
    runtime = get_runtime(parse_bam_read, filebamIn, "out")
    print(runtime)
    # runtime = get_runtime(parse_bam_read, filebamIn, 'out')
    # print(runtime)

    DB_NAME, TABLE_NAME = parse_bam_read(filebamIn, "out")

    # DB_NAME = "out/winnowmap_guppy_merge_subset.db"
    TABLE_NAME = "reads"
    plot_feature_df = pd.read_sql(
        "SELECT * from " + TABLE_NAME, con=sqlite3.connect(DB_NAME)
    )
    print(plot_feature_df.columns)
    fig, ax = plt.subplots(figsize=(6, 4))
    x = plot_feature_df["length"]
    # sns.displot(data=x, color ='#559CAD', kind='hist', bins = 200)
    # plt.show()
    # plt.hist(x, bins=200, color = '#559CAD', density=True) #density = True
    colors = [
        "#BB4430",
        "#FFBC0A",
        "#053C5E",
        "#A9E5BB",
        "#610345",
        "#2D1E2F",
        "#559CAD",
        "#5E747F",
        "#F343F4",
    ]
    plt.hist(x, bins=200, color=colors[6])  # density = True
    plt.axvline(
        x.median(),
        color=colors[0],
        linestyle="dashed",
        linewidth=1.3,
        label="median: " + str(round(x.median())) + " bp",
    )
    plt.axvline(
        x.mean(),
        color=colors[2],
        linestyle="dashed",
        linewidth=1.3,
        label="mean: " + str(round(x.mean())) + " bp",
    )
    # sns.despine()
    plt.title(DB_NAME.split("/")[1][:-3] + " read length")
    plt.xlabel("Read Length (bp)")
    plt.ylabel("Count")
    plt.legend()
    plt.savefig(DB_NAME[:-3] + "_read_length_count.pdf")
    plt.clf()

    plt.hist(x, bins=200, color=colors[6], density=True)  #
    plt.axvline(
        x.median(),
        color=colors[0],
        linestyle="dashed",
        linewidth=1.3,
        label="median: " + str(round(x.median())) + " bp",
    )
    plt.axvline(
        x.mean(),
        color=colors[2],
        linestyle="dashed",
        linewidth=1.3,
        label="mean: " + str(round(x.mean())) + " bp",
    )
    # sns.despine()
    plt.title(DB_NAME.split("/")[1][:-3] + " read length")
    plt.xlabel("Read Length (bp)")
    plt.ylabel("Frequency")
    plt.legend()
    plt.savefig(DB_NAME[:-3] + "_read_length_freq.pdf")


def main():
    print("main")
