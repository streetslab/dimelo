"""
QC Report Example
===========================

This creates a qc report from the test bams

"""
# %%
# Creating QC Report
# ------------------------
# Usually as a first step after receiving bam files, we want to do a quality check
# and assess our data. This code generates a PDF report of important QC statistics
# for long read methylation modification bam files

import dimelo as dm

# first we specify the locations of our bam files
in_bam = "../../dimelo/test/data/winnowmap_guppy_merge_subset.bam"
sample_name = "test"
out_dir = "../../dimelo/dimelo_test"

# next we run the "qc_report" function
dm.qc_report(in_bam, sample_name, out_dir)
# now our output directory that we specified, will have a file called
# "test_qc_report.pdf"