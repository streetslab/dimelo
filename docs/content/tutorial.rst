Tutorial
====================================

Follow the steps here to run the dimelo package with a test dataset.

Files needed
-------------
Files can be found in directory: ``dimelo/test/data``

1. ctcf_demo.sorted.bam
2. ctcf_demo.sorted.bam.bai
3. ctcf_demo_peak.bed
4. ctcf_demo_not_peak.bed

Steps
-------------
1. :doc:`../content/installation`
2. Run below commands and ensure your output matches the expected plots below. This tutorial walks through running functions from python, but the dimelo package can also be used from the command line (see :doc:`../auto_examples/index`).

>>> import dimelo as dm
>>> bam = "dimelo/test/data/ctcf_demo.sorted.bam"
>>> sampleName = "CTCF_demo"
>>> outDir = "out"
>>> bedPeak = "dimelo/test/data/ctcf_demo_peak.bed"
>>> bedNotPeak = "dimelo/test/data/ctcf_demo_not_peak.bed"

>>> dm.qc_report(bam, sampleName, outDir)

.. image:: demo/CTCF_demo_qc_report.png
    :align: center

>>> dm.plot_enrichment(bam, ['CTCF_peak', 'CTCF_not_peak'], [bedPeak, bedNotPeak], "A", outDir, threshA=190)

.. image:: demo/sample_ctcf_demo.sorted_A_enrichment_barplot.png
    :align: center

>>> dm.plot_enrichment_profile(bam, sampleName, bedPeak,"A+CG",outDir,threshA=190,threshC=190)

Aggregate rolling average fraction of methylated bases

.. image:: demo/CTCF_demo_A+CG_sm_rolling_avg.png
    :align: center

Single molecules with binary mA and mCpG colored

.. image:: demo/CTCF_demo_A+CG_sm_scatter.png
    :align: center

A base count

.. image:: demo/CTCF_demo_A_base_count.png
    :align: center

CG base count

.. image:: demo/CTCF_demo_CG_base_count.png
    :align: center

>>> dm.plot_browser(bam, sampleName, "chr11:2086423-2091187", "A+CG", outDir, threshA=153, threshC=153, static=True, smooth=100, min_periods=10)

.. image:: demo/methylation_browser_chr11_2086423_2091187.png
    :align: center

.. image:: demo/CTCF_demo_A_sm_rolling_avg_fraction.png
    :align: center

.. image:: demo/CTCF_demo_A_sm_rolling_avg_total.png
    :align: center

.. image:: demo/CTCF_demo_CG_sm_rolling_avg_fraction.png
    :align: center

.. image:: demo/CTCF_demo_CG_sm_rolling_avg_total.png
    :align: center