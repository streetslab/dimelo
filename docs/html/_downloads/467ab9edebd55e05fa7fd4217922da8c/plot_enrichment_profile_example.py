"""
Enrichment Profile Plot Example
===========================

This plots single molecules centered at regions of interest defined in bed file and produces aggregate profile
"""
# %%
# Plotting Enrichment Profile
# ------------------------
# add description here

import dimelo as dm

input_bam = "../../dimelo/test/data/mod_mappings_subset.bam"
sample_names = ["test1", "test2"]
regions = [
    "../../dimelo/test/data/test.bed",
    "../../dimelo/test/data/test.bed",
]
outDir = "../../dimelo/dimelo_test"
dm.plot_enrichment_profile(
    input_bam, sample_names, regions, "A", outDir, windowSize=500, dotsize=1
)

# This will return an aggregate profile of fraction of bases
# modified centered at features of interest
