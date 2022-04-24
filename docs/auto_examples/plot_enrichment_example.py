"""
Enrichment Plot Example
===========================

This plots the fraction of bases modified within regions of interest defined by a bed file
"""
# %%
# Plotting Enrichment
# ------------------------
# add description here

import dimelo as dm

input_bams = [
    "../../dimelo/test/data/mod_mappings_subset.bam",
    "../../dimelo/test/data/mod_mappings_subset.bam",
]
sample_names = ["test1", "test2"]
regions = "../../dimelo/test/data/test.bed"
outDir = "../../dimelo/dimelo_test"
dm.plot_enrichment(
    input_bams, sample_names, regions, "CG", outDir, threshC=129
)

# This will return a barplot with overall fraction of bases
# modified within regions of interest specified by bedFile(s)
