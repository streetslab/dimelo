r"""
=================
plot_browser module
=================
.. currentmodule:: dimelo.plot_browser
.. autosummary::
    plot_browser

plot_browser plots single molecules with colored base modifications in region of interest


Portions of code adapted from methplotlib:
Copyright (c) 2018 Wouter De Coster
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

"""

# code adapted from methplotlib
# https://doi.org/10.1093/bioinformatics/btaa093

import argparse
import multiprocessing
import os
import sqlite3
import sys

import matplotlib.pyplot as plt
import pandas as pd
import plotly
import plotly.graph_objs as go
import pyranges as pr
import seaborn as sns

from dimelo.parse_bam import parse_bam

# import plotly.io as pio

COLOR_A = "#053C5E"
COLOR_C = "#BB4430"
DEFAULT_THRESH_A = 129
DEFAULT_THRESH_C = 129
DEFAULT_SMOOTH = 1000
DEFAULT_MIN_PERIODS = 100
DEFAULT_DOTSIZE = 4


class DataTraces(object):
    def __init__(self, traces, names):
        self.traces = traces
        self.names = names
        self.index = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.index == len(self.traces):
            raise StopIteration
        else:
            self.index += 1
            return self.traces[self.index - 1]


class Region(object):
    def __init__(self, region, fasta=None):
        if ":" in region:
            try:
                self.chromosome, interval = region.replace(",", "").split(":")
                self.begin, self.end = [int(i) for i in interval.split("-")]
            except ValueError:
                sys.exit(
                    "\n\nERROR: Region (-w/--region) inproperly formatted, "
                    "examples of accepted formats are:\n"
                    "'chr5:150200605-150423790'\n\n"
                )
            self.size = self.end - self.begin
            self.string = f"{self.chromosome}_{self.begin}_{self.end}"


def plot_browser(
    fileNames,
    sampleNames,
    region,
    basemod,
    outDir,
    threshA=DEFAULT_THRESH_A,
    threshC=DEFAULT_THRESH_C,
    bedFileFeatures=None,
    smooth=DEFAULT_SMOOTH,
    min_periods=DEFAULT_MIN_PERIODS,
    colorA=COLOR_A,
    colorC=COLOR_C,
    dotsize=DEFAULT_DOTSIZE,
    static=False,
    cores=None,
):
    """
    fileNames
        list of names of bam files with Mm and Ml tags; indexed; or single file name as string
    sampleNames
        list of names of samples for output plot name labelling; or single sample name as string; valid names contain [``a-zA-Z0-9_``].
    region
        formatted as for example: "chr1:1-100000"
    basemod
        One of the following:

        * ``'A'`` - extract mA only
        * ``'CG'`` - extract mCpG only
        * ``'A+CG'`` - extract mA and mCpG
    outDir
        directory to output plot
    threshA
        threshold for calling mA; default 129
    threshC
        threshold for calling mCG; default 129
    bedFileFeatures
        bedFile specifying regions to display in browser (optional); default None
    smooth
        window over which to smooth aggregate curve; default of 1000 bp
    min_periods
        minimum number of bases to consider for smoothing: default of 100 bp
    colorA
        color in hex for mA; default #053C5E
    colorC
        color in hex for mCG; default #BB4430
    dotsize
        size of points; default 4
    static
        One of the following:

        * ``'True'`` - pdf output
        * ``'False'`` - interactive html output; default is False
    cores
        number of cores over which to parallelize; default is all available

    **Example**

    >>> dm.plot_browser("dimelo/test/data/mod_mappings_subset.bam", "test", "chr1:2907273-2909473", "A+CG", "dimelo/dimelo_test", static=False)
    >>> dm.plot_browser(["dimelo/test/data/mod_mappings_subset.bam", "dimelo/test/data/mod_mappings_subset.bam"], ["test1", "test2"], "chr1:2907273-2909473", "A+CG", "dimelo/dimelo_test", static=False)

    **Return**

        * PDF or HTML file with single molecules displayed over region of interest. Modified bases are colored according to colorA and colorC.
        * PDFs of aggregate coverage and fraction of bases modified over region of interest.
        * A summary bed file is also produced to support visualizing aggregate data with any genome browser tool. The columns of this bed file are chr, start, end, methylated_bases, total_bases. For example, to take a summary output bed and create a file with fraction of modified bases with a window size of 100 bp for visualization with the WashU browser, you could run the below commands in terminal:

            * ``bedtools makewindows -g ref_genome.chromsizes.txt -w 100 > ref_genome_windows.100.bp.bed``
            * ``bedtools map -a ref_genome_windows.100.bp.bed -b outDir/fileName_sampleName_chr_start_end_A.bed -c 4,5 -o sum,sum -null 0 | awk -v "OFS=\\t" '{if($5>0){print $1,$2,$3,$4/$5}else{print $1,$2,$3,$5}}' > outDir/fileName_sampleName_chr_start_end_A.100.bed``
            * ``bgzip outDir/fileName_sampleName_chr_start_end_A.100.bed``
            * ``tabix -f -p bed outDir/fileName_sampleName_chr_start_end_A.100.bed.gz``

    **Example Plots**

    :ref:`sphx_glr_auto_examples_browser_example.py`
    """

    if not os.path.isdir(outDir):
        os.makedirs(outDir)

    cores_avail = multiprocessing.cpu_count()
    if cores is None:
        num_cores = cores_avail
    else:
        # if more than available cores is specified, process with available cores
        if cores > cores_avail:
            num_cores = cores_avail
        else:
            num_cores = cores

    # if  single bam file rather than list is entered, convert to list
    if type(fileNames) != list:
        fileNames = [fileNames]
    # if single sample name rather than list is entered, convert to list
    if type(sampleNames) != list:
        sampleNames = [sampleNames]

    all_data = []
    aggregate_counts = []
    for f, n in zip(fileNames, sampleNames):
        # extract all bases and both mods to get full extent of read in terms of any mod bases
        parse_bam(
            f,
            n,
            outDir,
            basemod="A+CG",
            region=region,  # pass string representation
            threshA=threshA,
            threshC=threshC,
            extractAllBases=True,
            cores=num_cores,
        )
        d = pd.read_sql(
            "SELECT * from methylationByBase_" + n,
            sqlite3.connect(
                outDir + "/" + f.split("/")[-1].replace(".bam", "") + ".db"
            ),
        )
        all_data.append(d)
        aggregate_counts.append(
            pd.read_sql(
                "SELECT * from methylationAggregate_" + n,
                sqlite3.connect(
                    outDir + "/" + f.split("/")[-1].replace(".bam", "") + ".db"
                ),
            )
        )

        # print number of reads for each sample
        print(
            "processing "
            + str(len(d["read_name"].unique()))
            + " reads for "
            + n
            + " for bam: "
            + f
        )

    meth_browser(
        all_data=all_data,
        aggregate_counts=aggregate_counts,
        basemod=basemod,
        region=Region(region),
        sampleNames=sampleNames,
        outDir=outDir,
        bed=bedFileFeatures,
        smooth=smooth,
        min_periods=min_periods,
        dotsize=dotsize,
        static=static,
        threshA=threshA,
        threshC=threshC,
        colorA=colorA,
        colorC=colorC,
    )

    # print output files to std out
    if static is False:
        ext = "html"
    else:
        ext = "pdf"

    db_paths = []
    f_paths = []
    t_paths = []
    b_paths = []

    for f, s in zip(fileNames, sampleNames):
        db = outDir + "/" + f.split("/")[-1].replace(".bam", "") + ".db"
        db_paths.append(db)
        f_base = f.split("/")[-1].replace(".bam", "")
        if "A" in basemod:
            f_path = (
                outDir + "/" + s + "_" + "A" + "_sm_rolling_avg_fraction.pdf"
            )
            t_path = outDir + "/" + s + "_" + "A" + "_sm_rolling_avg_total.pdf"
            b_path = f"{outDir}/{f_base}_{s}_{Region(region).string}_A.bed"
            f_paths.append(f_path)
            t_paths.append(t_path)
            b_paths.append(b_path)
        if "C" in basemod:
            f_path = (
                outDir + "/" + s + "_" + "CG" + "_sm_rolling_avg_fraction.pdf"
            )
            t_path = (
                outDir + "/" + s + "_" + "CG" + "_sm_rolling_avg_total.pdf"
            )
            b_path = f"{outDir}/{f_base}_{s}_{Region(region).string}_CG.bed"
            f_paths.append(f_path)
            t_paths.append(t_path)
            b_paths.append(b_path)

    w = Region(region)

    browser_path = f"{outDir}/methylation_browser_{w.string}.{ext}"
    str_out = f"Outputs\n_______\nDB file: {db_paths}\nbrowser plot: {browser_path}\nrolling average fraction bases methylated plot: {f_paths}\nrolling average total bases plot: {t_paths}\nsummary bed file: {b_paths}"
    print(str_out)


def create_subplots(num_methrows, names=None, annotation=True):
    """
    Prepare the panels (rows * 1 column) for the subplots.
    One row for each dataset, taking 90%/len(datasets) for heights
    if annotation is True (bed) then add a row with height 10%
    """
    return plotly.subplots.make_subplots(
        rows=num_methrows + annotation,
        cols=1,
        shared_xaxes=True,
        specs=[[{}] for i in range(num_methrows + annotation)],
        print_grid=False,
        subplot_titles=names,
        vertical_spacing=0.1 if num_methrows < 10 else 0.01,
        row_heights=[0.9 / num_methrows] * num_methrows + [0.1] * annotation,
    )


def create_output(fig, outfile, region, static, outDir):
    """
    write output pdf or html
    """
    if static:
        outfile = outDir + "/" + f"methylation_browser_{region.string}.pdf"
        fig.write_image(outfile, width=1000, height=400)
    if not static:
        outfile = outDir + "/" + f"methylation_browser_{region.string}.html"
        with open(outfile, "w+") as output:
            output.write(
                plotly.offline.plot(
                    fig,
                    output_type="div",
                    show_link=False,
                    include_plotlyjs="cdn",
                )
            )


def methylation(
    all_data,
    sampleNames,
    basemod,
    colorA=COLOR_A,
    colorC=COLOR_C,
    dotsize=4,
    threshA=DEFAULT_THRESH_A,
    threshC=DEFAULT_THRESH_C,
):
    """
    Plot methylation traces
    """
    traces = []
    names = []
    for m, n in zip(all_data, sampleNames):
        traces.append(
            make_per_read_meth_traces_phred(
                all_data=all_data,
                table=m,
                basemod=basemod,
                colorA=colorA,
                colorC=colorC,
                dotsize=dotsize,
                threshA=threshA,
                threshC=threshC,
            )
        )
        names.append(n)
    return DataTraces(traces=traces, names=names)


def make_per_read_meth_traces_phred(
    all_data,
    table,
    basemod,
    colorA,
    colorC,
    max_cov=1000,
    dotsize=4,
    threshA=DEFAULT_THRESH_A,
    threshC=DEFAULT_THRESH_C,
):
    """Make traces for each read"""
    minmax_table = find_min_and_max_pos_per_read(table)
    df_heights = assign_y_height_per_read(minmax_table, max_coverage=max_cov)
    table = pd.merge(table, df_heights, left_on="read_name", right_on="read")
    traces = []
    hidden = 0
    for read in table["read_name"].unique():
        try:
            traces.append(
                make_per_read_line_trace(
                    read_range=minmax_table.loc[read],
                    y_pos=df_heights.loc[read, "height"],
                )
            )
        except KeyError:
            hidden += 1
            continue
    if hidden:
        sys.stderr.write(
            f"Warning: hiding {hidden} reads because coverage above {max_cov}x.\n"
        )
    read_table_mC = table[table["mod"].str.contains("C")]
    read_table_mA = table[table["mod"].str.contains("A")]
    cmapA = ["white", colorA]
    cmapC = ["white", colorC]
    if "C" in basemod:
        traces.append(
            make_per_position_phred_scatter(
                all_data=all_data,
                read_table=read_table_mC[read_table_mC["prob"] >= threshC],
                mod="mC",
                thresh=threshC,
                dotsize=dotsize,
                colorscale=cmapC,
                offset=0.05,
            )
        )
    if "A" in basemod:
        traces.append(
            make_per_position_phred_scatter(
                all_data=all_data,
                read_table=read_table_mA[read_table_mA["prob"] >= threshA],
                mod="mA",
                thresh=threshA,
                dotsize=dotsize,
                colorscale=cmapA,
                offset=0.15,
            )
        )
    return traces


def make_per_position_phred_scatter(
    all_data, read_table, mod, thresh, dotsize=4, colorscale="Reds", offset=0
):
    """Make scatter plot per modified base per read"""
    # get min and max probabilities across all for legend and color consistency for comparisons
    if "C" in mod:
        m = "C"
    if "A" in mod:
        m = "A"
    min_overall = 255
    max_overall = 0
    for d in all_data:
        min_temp = d[(d["mod"].str.contains(m)) & (d["prob"] >= thresh)][
            "prob"
        ].min()
        max_temp = d[(d["mod"].str.contains(m)) & (d["prob"] >= thresh)][
            "prob"
        ].max()
        if min_temp < min_overall:
            min_overall = min_temp
        if max_temp > max_overall:
            max_overall = max_temp
    return go.Scatter(
        x=read_table["pos"],
        y=read_table["height"],
        mode="markers",
        showlegend=False,
        text=round(read_table["prob"] / 255, 2),
        hoverinfo="text",
        marker=dict(
            size=dotsize,
            color=read_table["prob"],
            colorscale=colorscale,
            colorbar=dict(
                title=mod + " probability",
                titleside="right",
                tickmode="array",
                tickvals=[min_overall, max_overall],
                ticktext=[
                    str(round(min_overall / 255, 2)),
                    str(round(max_overall / 255, 2)),
                ],
                ticks="outside",
                thickness=15,
                x=offset + 1,
            ),
        ),
    )


def find_min_and_max_pos_per_read(table):
    """Return a table with for every read the minimum and maximum position"""
    mm_table = (
        table.loc[:, ["read_name", "pos"]]
        .groupby("read_name")
        .min()
        .join(
            table.loc[:, ["read_name", "pos"]].groupby("read_name").max(),
            lsuffix="min",
            rsuffix="max",
        )
    )
    return mm_table


def assign_y_height_per_read(df, max_coverage=1000):
    """Assign height of the read in the per read traces
    Gets a dataframe of read_name, posmin and posmax.
    Sorting by position.
    Determines optimal height (y coordinate) for this read
    Returns a dictionary mapping read_name to y_coord
    """
    dfs = df.sort_values(by=["posmin", "posmax"], ascending=[True, False])
    heights = [[] for i in range(max_coverage)]
    y_pos = dict()
    for read in dfs.itertuples():
        for y, layer in enumerate(heights, start=1):
            if len(layer) == 0:
                layer.append(read.posmax)
                y_pos[read.Index] = y
                break
            if read.posmin > layer[-1]:
                layer.append(read.posmax)
                y_pos[read.Index] = y
                break
    return pd.DataFrame(
        {"read": list(y_pos.keys()), "height": list(y_pos.values())}
    ).set_index("read")


def make_per_read_line_trace(read_range, y_pos):
    """
    Make a grey line trace for a single read
    """
    return go.Scatter(
        x=[read_range["posmin"], read_range["posmax"]],
        y=[y_pos, y_pos],
        mode="lines",
        line=dict(width=1, color="lightgrey"),
        showlegend=False,
    )


def meth_browser(
    all_data,
    aggregate_counts,
    basemod,
    region,
    sampleNames,
    outDir,
    smooth,
    min_periods,
    bed=False,
    outfile=None,
    dotsize=4,
    static=False,
    threshA=DEFAULT_THRESH_A,
    threshC=DEFAULT_THRESH_C,
    colorA=COLOR_A,
    colorC=COLOR_C,
):
    """
    meth_data is a list of methylationByBase tables as dataframes
    all_dict is a list of methylationAggregate tables as dataframes
    annotation is optional and is a bed file
     then show one line per sample and one for the annotation, with methrows = number of datasets
    the trace to be used for annotation is thus always num_methrows + 1
    """
    meth_traces = methylation(
        all_data,
        sampleNames,
        basemod,
        colorA=colorA,
        colorC=colorC,
        dotsize=dotsize,
        threshA=threshA,
        threshC=threshC,
    )

    num_methrows = len(all_data)
    annot_row = num_methrows + 1
    annot_axis = f"yaxis{annot_row}"
    fig = create_subplots(
        num_methrows, names=meth_traces.names, annotation=bool(bed)
    )
    for y, sample_traces in enumerate(meth_traces, start=1):
        for meth_trace in sample_traces:
            fig.add_trace(trace=meth_trace, row=y, col=1)
        fig["layout"][f"yaxis{y}"].update(title="Reads")
    if bed:
        for annot_trace in bed_annotation(bed, region):
            fig.add_trace(trace=annot_trace, row=annot_row, col=1)
        y_max = -2
    if bed:
        fig["layout"][annot_axis].update(
            range=[-2, y_max + 1],
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks="",
            showticklabels=False,
        )
    fig["layout"]["xaxis"].update(
        tickformat="g",
        separatethousands=True,
        range=[region.begin, region.end],
    )
    fig["layout"].update(
        barmode="overlay",
        title=region.chromosome,
        hovermode="closest",
        plot_bgcolor="rgba(0,0,0,0)",
    )
    if num_methrows > 10:
        for i in fig["layout"]["annotations"]:
            i["font"]["size"] = 10
    create_output(fig, outfile, region, static, outDir)

    i = 0
    for d in aggregate_counts:
        plot_aggregate(
            sampleNames[i],
            d,
            smooth,
            min_periods,
            region,
            basemod,
            outDir,
            colorA,
            colorC,
        )
        i = i + 1


def bed_annotation(bed, region):
    return [
        go.Scatter(
            x=[begin, end],
            y=[-2, -2],
            mode="lines",
            line=dict(width=16, color="grey"),
            text=name,
            hoverinfo="text",
            showlegend=False,
        )
        for (begin, end, name) in parse_bed(bed, region)
    ]


def parse_bed(bed, region):
    gr = pr.read_bed(bed)[region.chromosome, region.begin : region.end]
    df = gr.unstrand().df
    df = df.drop(columns=["Chromosome", "Score", "Strand"], errors="ignore")
    if "Name" not in df.columns:
        df["Name"] = "noname"
    df_short = df[df.columns[0:3]]
    return df_short.itertuples(index=False, name=None)


def plot_aggregate(
    sampleName,
    aggregate_counts,
    smooth,
    min_periods,
    region,
    basemod,
    outDir,
    colorA,
    colorC,
):
    """
    plot rolling aggregate of frac methylated
    plot rolling aggregate of total bases
    """

    aggregate_counts["frac"] = (
        aggregate_counts["methylated_bases"] / aggregate_counts["total_bases"]
    )

    # plot aggregate of fraction and of total count coverage
    if "A" in basemod:
        aggregate_A = aggregate_counts[
            aggregate_counts["mod"].str.contains("A")
        ].copy()
        # need to sort first!
        aggregate_A.sort_values(["pos"], inplace=True)
        aggregate_A_rolling = aggregate_A.rolling(
            window=smooth, min_periods=min_periods, center=True, on="pos"
        ).mean()
        plot_aggregate_frac(
            aggregate_A_rolling, sampleName, "A", colorA, outDir
        )
        plot_aggregate_total(
            aggregate_A_rolling, sampleName, "A", colorA, outDir
        )
    if "C" in basemod:
        aggregate_C = aggregate_counts[
            aggregate_counts["mod"].str.contains("C")
        ].copy()
        # need to sort first!
        aggregate_C.sort_values(["pos"], inplace=True)
        aggregate_C_rolling = aggregate_C.rolling(
            window=smooth, min_periods=min_periods, center=True, on="pos"
        ).mean()
        plot_aggregate_frac(
            aggregate_C_rolling, sampleName, "C", colorC, outDir
        )
        plot_aggregate_total(
            aggregate_C_rolling, sampleName, "C", colorC, outDir
        )


def plot_aggregate_frac(aggregate_rolling, sampleName, mod, color, outDir):
    fig = plt.figure()
    sns.lineplot(
        x=aggregate_rolling["pos"],
        y=aggregate_rolling["frac"],
        color=color,
    )
    if "A" in mod:
        mod_name = "A"
    if "C" in mod:
        mod_name = "CG"

    plt.title(mod_name)
    plt.ylabel("m" + mod_name + "/" + mod_name)
    fig.savefig(
        outDir
        + "/"
        + sampleName
        + "_"
        + mod_name
        + "_sm_rolling_avg_fraction.pdf"
    )
    plt.close()


def plot_aggregate_total(aggregate_rolling, sampleName, mod, color, outDir):
    fig = plt.figure()
    sns.lineplot(
        x=aggregate_rolling["pos"],
        y=aggregate_rolling["total_bases"],
        color=color,
    )
    if "A" in mod:
        mod_name = "A"
    if "C" in mod:
        mod_name = "CG"

    plt.title(mod_name)
    plt.ylabel("total " + mod_name)
    fig.savefig(
        outDir
        + "/"
        + sampleName
        + "_"
        + mod_name
        + "_sm_rolling_avg_total.pdf"
    )
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="DiMeLo plot browser")

    # Required arguments
    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument(
        "-f", "--fileNames", required=True, nargs="+", help="bam file name(s)"
    )
    required_args.add_argument(
        "-s",
        "--sampleNames",
        required=True,
        nargs="+",
        help="sample name(s) for output file labelling",
    )
    required_args.add_argument(
        "-r",
        "--region",
        required=True,
        type=str,
        help='single region over which to extract base mods, e.g. "chr1:1-100000"',
    )
    required_args.add_argument(
        "-m",
        "--basemod",
        required=True,
        type=str,
        choices=["A", "CG", "A+CG"],
        help="which base modification to extract",
    )
    required_args.add_argument(
        "-o", "--outDir", required=True, help="directory to output plot"
    )

    # Smoothing options
    smoothing_args = parser.add_argument_group("smoothing options")
    smoothing_args.add_argument(
        "-t",
        "--smooth",
        type=int,
        default=DEFAULT_SMOOTH,
        help="window over which to smooth aggregate curve",
    )
    smoothing_args.add_argument(
        "-n",
        "--min_periods",
        type=int,
        default=DEFAULT_MIN_PERIODS,
        help="minimum number of bases to consider for smoothing",
    )

    # Plotting arguments
    plotting_args = parser.add_argument_group("plotting options")
    plotting_args.add_argument(
        "--colorA",
        type=str,
        default=COLOR_A,
        help='color in hex (e.g. "#BB4430") for mA',
    )
    plotting_args.add_argument(
        "--colorC",
        type=str,
        default=COLOR_C,
        help='color in hex (e.g. "#BB4430") for mCG',
    )
    plotting_args.add_argument(
        "-d",
        "--dotsize",
        type=float,
        default=DEFAULT_DOTSIZE,
        help="size of points",
    )

    # Optional arguments
    parser.add_argument(
        "-A",
        "--threshA",
        type=int,
        default=DEFAULT_THRESH_A,
        help="threshold above which to call an A base methylated",
    )
    parser.add_argument(
        "-C",
        "--threshC",
        type=int,
        default=DEFAULT_THRESH_C,
        help="threshold above which to call a C base methylated",
    )
    parser.add_argument(
        "-b",
        "--bedFileFeatures",
        help="bed file specifying annotation to display in browser",
    )
    parser.add_argument(
        "--static",
        action="store_true",
        help="output as PDF instead of interactive HTML",
    )
    parser.add_argument(
        "-p",
        "--cores",
        type=int,
        help="number of cores over which to parallelize",
    )

    args = parser.parse_args()
    plot_browser(**vars(args))
