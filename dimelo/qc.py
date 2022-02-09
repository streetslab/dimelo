"Functions to parse_bam for QC"

import numpy as np
import pandas as pd
import pysam
import multiprocessing
import time
import sqlite3
from joblib import Parallel, delayed
from math import log
import seaborn as sns
import matplotlib.pyplot as plt
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


def batch_read_generator(file_bamIn, batch_size):
    counter = 0
    r_list = []

    for read in file_bamIn.fetch(until_eof=True):
        #isA = False
        #isC = False
        r = [read.query_name, read.reference_name, read.reference_start, read.reference_end, read.query_length,
             "-" if read.is_reverse else "+", read.mapping_quality, ave_qual(read.query_qualities),
             ave_qual(read.query_alignment_qualities)]
        #[(mod, numA, probs), (mod2, numC, probs2)] = get_modified_reference_positions(read, "A+CG", None, False,
         #                                                                             129, 129, None, None, None,
          #                                                                            outDir, False, True)

       # if type(numA) == int:
       #      isA = True
       #      problist = [int(x) / 255 for x in probs]
       #      probagg_05 = prob_bin(problist)
       #
       #  if type(numC) == int:
       #      isC = True
       #      prob2list = [int(x) / 255 for x in probs2]
       #      prob2agg_05 = prob_bin(prob2list)
       #
       #  if isC and not isA:
       #      r.extend((np.nan, np.nan, numC, prob2agg_05))
       #  elif isA and not isC:
       #      r.extend((numA, probagg_05, np.nan, np.nan))
       #  elif isA and isC:
       #      r.extend((numA, probagg_05, numC, prob2agg_05))

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
    probs = [np.log(1 - p) for p in bin if
             ((p < 1) and (p >= 0.5))]  # only consider probabilities > 0.5 and handle 1 on next line
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


def parse_bam_read(bamIn, outDir, cores):
    file_bamIn = pysam.AlignmentFile(bamIn, "rb")

    # DB_NAME = bamIn.split('/')[-1].split('.')[0] + ".db"
    # cols = ['name', 'chr', 'start', 'end', 'length', 'strand', 'mapq', 'ave_baseq', 'ave_alignq', 'basemod',
    #        'numA', 'methAprob05', 'numC', 'methCprob05']
    # dt = ['TEXT', 'TEXT', 'INT', 'INT', 'INT', 'TEXT', 'INT', 'INT', 'INT',
    # 'TEXT', 'FLOAT', 'FLOAT', 'FLOAT', 'FLOAT']
    # table_name = 'READS'
    # create_sql_table(DB_NAME, table_name, cols, dt)

    DB_NAME, tables = make_db(bamIn, '', outDir, True, True)
    template_command = '''INSERT INTO ''' + tables[0] + ''' VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?);'''
    #num_cores = multiprocessing.cpu_count() - 1 #need to keep all cores
    cores_avail = multiprocessing.cpu_count()
    if cores is None:
        num_cores = cores_avail
    else:
        # if more than available cores is specified, process with available cores
        if cores > cores_avail:
            num_cores = cores_avail
        else:
            num_cores = cores

    Parallel(n_jobs=num_cores, verbose=10)(
        delayed(execute_sql_command)(
            template_command,
            DB_NAME,
            i) for i in batch_read_generator(file_bamIn, 100))
    return DB_NAME, tables[0]


def get_runtime(f, inp1, inp2):
    start = time.time()
    re_val = f(inp1, inp2)
    time.sleep(1)
    end = time.time()
    return f"Runtime of the program is {end - start}", re_val


def qc_plot(x, sampleName, plotType, colors, num, axes):
    #plt.subplot(3, 2, num)
    an_array = np.array(x)
    q1 = np.quantile(an_array, 0.25)
    q3 = np.quantile(an_array, 0.75)
    iq = q3 - q1
    outlier = q3 + 3 * iq
    not_outlier = an_array <= outlier
    no_outliers = an_array[not_outlier]


    ptype = ''
    unit = ''
    xlabel = ''
    hasN50 = False
    if plotType == 'L':
        ptype = " Read Length"
        xlabel = "Read Length (bp)"
        unit = " bp"
        hasN50 = True
        n50 = calculate_N50(x)
    elif plotType == 'M':
        ptype = " Mapping Quality"
        xlabel = "Mapping Quality"
    elif plotType == 'B':
        ptype = " Basecall Quality"
        xlabel = "Average Basecall Quality"
    elif plotType == 'A':
        ptype = " Alignment Quality"
        xlabel = "Average Alignment Quality"


    plt.hist(no_outliers, bins=200, color=colors[6], density=True)  #
    plt.axvline(x.median(), color=colors[0], linestyle='dashed', linewidth=1.3,
                label='median: ' + str(round(x.median())) + unit)
    plt.axvline(x.mean(), color=colors[2], linestyle='dashed', linewidth=1.3,
                label='mean: ' + str(round(x.mean())) + unit)
    if hasN50:
        plt.axvline(n50, color=colors[1], linestyle='dashed', linewidth=1,
                label='N50: ' + str(round(n50)) + unit)
    plt.title(ptype)
    plt.xlabel(xlabel)
    plt.ylabel("Frequency")
    plt.xlim(0, )
    plt.plot([], [], ' ', label="max: " + str(round(max(x))) + unit)
    plt.legend()

    values = [round(min(x)), round(q1), round(x.median()),
              round(q3), round(max(x)), round(x.mean())]
    return plt, values
    #plt.savefig(outDir + "/" + sampleName + "_rlength_freq_no_outliers.pdf", bbox_inches='tight')

def calculate_N50(x):
    array_rl = np.array(x)
    N = np.sum(array_rl)
    array_rl[::-1].sort()
    rl_cumsum = np.cumsum(array_rl)
    n50 = array_rl[np.argmax(rl_cumsum > N / 2)]
    return n50


def qc_report(filebamInList, sampleNameList, outDir, testMode = False, colors = ["#BB4430", "#FFBC0A",
                                                                         "#053C5E", "#A9E5BB",
                                                                         "#610345", "#2D1E2F",
                                                                         "#559CAD", "#5E747F", "#F343F4"]):
    #runtime = get_runtime(parse_bam_read, filebamIn, 'out')
    #print(runtime)

    if type(filebamInList) != list:
        filebamInList = [filebamInList]
        sampleNameList = [sampleNameList]

    for index in range(len(filebamInList)):
        filebamIn = filebamInList[index]
        sampleName = sampleNameList[index]
        if testMode:
            DB_NAME = outDir + "/" + filebamIn.split("/")[-1][:-4] + ".db"
            #DB_NAME = "out/winnowmap_guppy_merge_subset.db"
            # DB_NAME = "out/mod_mappings_subset.db"
            TABLE_NAME = "reads"
        else:
            DB_NAME, TABLE_NAME = parse_bam_read(filebamIn, 'out')

        if sampleName is None:
            sampleName = DB_NAME.split('/')[1][:-3]


        plot_feature_df = pd.read_sql("SELECT * from " + TABLE_NAME, con=sqlite3.connect(DB_NAME))

        #fig, ax = plt.subplots(figsize=(10, 8))

        fig = plt.figure(figsize=(12,10))
        grid = plt.GridSpec(3, 2, figure=fig)

        ax_5 = plt.subplot2grid(shape=(3, 2), loc=(0, 0), colspan=2)
        ax_5.axis("off")
        pltTable = ax_5.table(cellText=report_table,
                              rowLabels=rows,
                              colLabels=columns,
                              loc='center')


        #gridsize = (3,2)
        #
        # fig = plt.figure()
        # fig.set_figheight(6)
        # fig.set_figwidth(6)
        # ax_1 = plt.subplot2grid(shape=(3, 2), loc=(0, 0), colspan=1)
        # ax_2 = plt.subplot2grid(shape=(3, 2), loc=(0, 1), colspan=1)
        # ax_3 = plt.subplot2grid(shape=(3, 2), loc=(1, 0), colspan=1)
        # #ax_5.axis("off")
        # ax_4 = plt.subplot2grid(shape=(3, 2), loc=(1, 1), colspan=1)
        # #ax_4.axis("off")
        # ax_5 = plt.subplot2grid(shape=(3, 2), loc=(2, 0), colspan=2)
        # ax_5.axis("off")

        #
        # report_table = np.array([[0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0]]).T
        # columns = ['Read Length', 'Mapping Quality', 'Basecall Quality', 'Alignment Quality']
        # rows = ['Min', '25%', 'Median', '75%', 'Max', 'Mean'] #get rid of N50, total bases (mean read length * total)
        # pltTable = ax_5.table(cellText=report_table,
        #                       rowLabels=rows,
        #                       colLabels=columns,
        #                       loc='center')
        #pltTable.scale(1, 1.5)


        #plt.show()

        #################
        # ax_5 = fig.add_subplot(111)
        # ax_5.axis("tight")
        # ax_5.axis("off")
        # ax_5.margins(3)
        #
        # Read Length
        #plt.figure()
        x = plot_feature_df['length']
        ax_1 = fig.add_subplot(grid[0,0])
        pltRL, valRL = qc_plot(x, sampleName, 'L', colors, 1, ax_1)
        pltRL.savefig(outDir + "/" + sampleName + "_rlength_freq_no_outliers.pdf", bbox_inches='tight')
        #
        # Mapping Quality
        #plt.figure()
        x = plot_feature_df['mapq']
        ax_2 = fig.add_subplot(grid[0,1])
        pltMQ, valMQ = qc_plot(x, sampleName, 'M', colors, 2, ax_2)
        pltMQ.savefig(outDir + "/" + sampleName + "_mapq_freq_no_outliers.pdf", bbox_inches='tight')

        # Basecall Quality
        #plt.figure()
        x = plot_feature_df['ave_baseq']
        ax_3 = fig.add_subplot(grid[1,0])
        pltBQ, valBQ = qc_plot(x, sampleName, 'B', colors, 3, ax_3)
        pltBQ.savefig(outDir + "/" + sampleName + "_baseq_freq_no_outliers.pdf", bbox_inches='tight')

        # Alignment Quality
        #plt.figure()
        x = plot_feature_df['ave_alignq']
        ax_4 = fig.add_subplot(grid[1,1])
        pltAQ, valAQ = qc_plot(x, sampleName, 'A', colors, 4, ax_4)
        pltAQ.savefig(outDir + "/" + sampleName + "_alignq_freq_no_outliers.pdf", bbox_inches='tight')

        report_table = np.array([valRL, valMQ, valBQ, valAQ]).T
        columns = ['Read Length', 'Mapping Quality', 'Basecall Quality', 'Alignment Quality']
        rows = ['Min', '25%', 'Median', '75%', 'Max', 'Mean'] #get rid of N50, total bases (mean read length * total)
        #have x num reads, x num bases description separate from the table
        print("mean length: ", valRL[5])
        print("num reads: ", len(x))
        print("num bases: ", round(valRL[5] * len(valRL)))
        print(report_table)

        #plt.table(report_table, loc='bottom')


        ##ax_5 = fig.add_subplot(grid[2,:])
        #ax_5 = plt.subplot2grid(gridsize, (2, 0), colspan=2, rowspan=2)
        # #ax_5.visible = False
        # ax_5.axis("tight")
        # #ax_5.axis("off")
        ax_5 = plt.subplot2grid(shape=(3, 2), loc=(2, 0), colspan=2)
        ax_5.axis("off")
        pltTable = ax_5.table(cellText=report_table,
                                rowLabels=rows,
                                colLabels=columns,
                                loc='center')
        #pltTable.scale(1, 1.5)
        #plt.subplots_adjust(wspace=0.4,
        #                    hspace=0.4)
        fig.tight_layout(w_pad = 2, h_pad = 4)
        #plt.subplots_adjust(left=0.3, bottom=0.4)
        #plt.title("", size=30)
        #pltTable.text(12, 3.4, '', size=30)
##############
        summary_data = "mean length: " + str(valRL[5]) + " bp"
        summary_data = summary_data + "; num reads: " + str(len(x))
        summary_data = summary_data + "; " + "num bases: " + str(round(valRL[5] * len(valRL)) + " bp")
        fig.suptitle(sampleName + " QC Summary Report", y = 1.05)
        #plt.title("mean length: " + str(valRL[5]), y = 0.8)
        plt.title(summary_data, y=0.8)
        #plt.title("TITLE: " + str(valRL[5]), fontsize = 12, y=1.3)
        plt.savefig(outDir + "/" + sampleName + "_qc_report.pdf", bbox_inches='tight')

    # an_array = np.array(x)
    # q1 = np.quantile(an_array, 0.25)
    # q3 = np.quantile(an_array, 0.75)
    # iq = q3 - q1
    # outlier = q3 + 3 * iq
    # not_outlier = an_array <= outlier
    # no_outliers = an_array[not_outlier]
    #
    #
    # n50 = calculate_N50(x)
    #
    # plt.hist(no_outliers, bins=200, color=colors[6], density = True)  #
    # plt.axvline(x.median(), color=colors[0], linestyle='dashed', linewidth=1.3,
    #             label='median: ' + str(round(x.median())) + ' bp')
    # plt.axvline(x.mean(), color=colors[2], linestyle='dashed', linewidth=1.3,
    #             label='mean: ' + str(round(x.mean())) + ' bp')
    # plt.axvline(n50, color=colors[1], linestyle='dashed', linewidth=1,
    #             label='N50: ' + str(round(n50)) + ' bp')
    # plt.title(sampleName + " read length")
    # plt.xlabel("Read Length (bp)")
    # plt.ylabel("Frequency")
    # plt.xlim(0, )
    # plt.plot([], [], ' ', label="max: " + str(round(max(x))) + ' bp')
    # plt.legend()
    # plt.savefig(outDir + "/" + sampleName + "_rlength_freq_no_outliers.pdf", bbox_inches='tight')
    #
    # # Basecall Quality
    # print(plot_feature_df.columns)
    # x = plot_feature_df['ave_baseq']
    #
    # an_array = np.array(x)
    # q1 = np.quantile(an_array, 0.25)
    # q3 = np.quantile(an_array, 0.75)
    # iq = q3 - q1
    # outlier = q3 + 3 * iq
    # not_outlier = an_array <= outlier
    # no_outliers = an_array[not_outlier]
    #
    # def calculate_N50(list_of_lengths):
    #     array_rl = np.array(list_of_lengths)
    #     N = np.sum(array_rl)
    #     array_rl[::-1].sort()
    #     rl_cumsum = np.cumsum(array_rl)
    #     n50 = array_rl[np.argmax(rl_cumsum > N / 2)]
    #     return n50
    #
    # n50 = calculate_N50(x)
    #
    # plt.figure()
    # plt.hist(no_outliers, bins=200, color=colors[6], density=True)  #
    # plt.axvline(x.median(), color=colors[0], linestyle='dashed', linewidth=1.3,
    #             label='median: ' + str(round(x.median())))
    # plt.axvline(x.mean(), color=colors[2], linestyle='dashed', linewidth=1.3,
    #             label='mean: ' + str(round(x.mean())))
    # plt.axvline(n50, color=colors[1], linestyle='dashed', linewidth=1,
    #             label='N50: ' + str(round(n50)))
    # plt.title(sampleName + " average basecall quality")
    # plt.xlabel("Average Basecall Quality")
    # plt.ylabel("Frequency")
    # #plt.xlim(0, )
    # plt.plot([], [], ' ', label="max: " + str(round(max(x))))
    # plt.legend()
    # plt.savefig(outDir + "/" + sampleName + "_baseq_freq_no_outliers.pdf", bbox_inches='tight')
    #
    # # Mapping Quality
    # x = plot_feature_df['mapq']
    #
    # an_array = np.array(x)
    # q1 = np.quantile(an_array, 0.25)
    # q3 = np.quantile(an_array, 0.75)
    # iq = q3 - q1
    # outlier = q3 + 3 * iq
    # not_outlier = an_array <= outlier
    # no_outliers = an_array[not_outlier]
    #
    # def calculate_N50(list_of_lengths):
    #     array_rl = np.array(list_of_lengths)
    #     N = np.sum(array_rl)
    #     array_rl[::-1].sort()
    #     rl_cumsum = np.cumsum(array_rl)
    #     n50 = array_rl[np.argmax(rl_cumsum > N / 2)]
    #     return n50
    #
    # n50 = calculate_N50(x)
    #
    # plt.figure()
    # plt.hist(no_outliers, bins=200, color=colors[6], density=True)  #
    # plt.axvline(x.median(), color=colors[0], linestyle='dashed', linewidth=1.3,
    #             label='median: ' + str(round(x.median())))
    # plt.axvline(x.mean(), color=colors[2], linestyle='dashed', linewidth=1.3,
    #             label='mean: ' + str(round(x.mean())))
    # plt.axvline(n50, color=colors[1], linestyle='dashed', linewidth=1,
    #             label='N50: ' + str(round(n50)))
    # plt.title(sampleName + " average mapping quality")
    # plt.xlabel("Average Mapping Quality")
    # plt.ylabel("Frequency")
    # plt.xlim(0, )
    # plt.plot([], [], ' ', label="max: " + str(round(max(x))))
    # plt.legend()
    # plt.savefig(outDir + "/" + sampleName + "_mapq_freq_no_outliers.pdf", bbox_inches='tight')
    #
    # # Alignment Quality
    # print(plot_feature_df['ave_alignq'])
    # x = plot_feature_df['ave_alignq']
    #
    # an_array = np.array(x)
    # q1 = np.quantile(an_array, 0.25)
    # q3 = np.quantile(an_array, 0.75)
    # iq = q3 - q1
    # outlier = q3 + 3 * iq
    # not_outlier = an_array <= outlier
    # no_outliers = an_array[not_outlier]
    #
    # def calculate_N50(list_of_lengths):
    #     array_rl = np.array(list_of_lengths)
    #     N = np.sum(array_rl)
    #     array_rl[::-1].sort()
    #     rl_cumsum = np.cumsum(array_rl)
    #     n50 = array_rl[np.argmax(rl_cumsum > N / 2)]
    #     return n50
    #
    # n50 = calculate_N50(x)
    #
    # plt.figure()
    # plt.hist(no_outliers, bins=200, color=colors[6], density=True)  #
    # plt.axvline(x.median(), color=colors[0], linestyle='dashed', linewidth=1.3,
    #             label='median: ' + str(round(x.median())))
    # plt.axvline(x.mean(), color=colors[2], linestyle='dashed', linewidth=1.3,
    #             label='mean: ' + str(round(x.mean())))
    # plt.axvline(n50, color=colors[1], linestyle='dashed', linewidth=1,
    #             label='N50: ' + str(round(n50)))
    # plt.title(sampleName + " average alignment quality")
    # plt.xlabel("Average Alignment Quality")
    # plt.ylabel("Frequency")
    # plt.xlim(0, )
    # plt.plot([], [], ' ', label="max: " + str(round(max(x))))
    # plt.legend()
    # plt.savefig(outDir + "/" + sampleName + "_alignq_freq_no_outliers.pdf", bbox_inches='tight')

def main():
    print("main")