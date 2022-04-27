"""
Enrichment Plot Comparison Across BEDs
======================================

Plot overall fraction of methylated bases within multiple sets of regions of interest specified by bed files for a single sample.

"""
# %%
# Create barplot comparing methylation levels in single sample across multiple regions of interest defined in bed files.

# %%
# 1. Python option
# ----------------
import dimelo as dm

bam = "deep_ctcf_mod_mappings_merge.sorted.bam"
beds = ["q10.150.slop.bed", "q10nopeak.bed"]
sampleNames = ["chip_peak", "not_chip_peak"]
outDir = "./out"
dm.plot_enrichment(bam, sampleNames, beds, "A", outDir, threshA=190)

# %%
# 2. Command line option
# ----------------------
# ``dimelo-plot-enrichment -f deep_ctcf_mod_mappings_merge.sorted.bam -s chip_peak not_chip_peak -b q10.150.slop.bed q10nopeak.bed -m A -o ./out -A 190``

# %%
# Output
# ----------------------
# .. figure:: ../auto_examples/images/sample_deep_ctcf_mod_mappings_merge.sorted_A_enrichment_barplot.png
