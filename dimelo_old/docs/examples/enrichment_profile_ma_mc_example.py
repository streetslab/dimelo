"""
Enrichment Profile mA & mCG
=================================

Aggregate and single molecule plots colored by modification and centered at regions of interest defined in bed file.

"""
# %%
# Create (1) aggregate profile plots for mA/A and mCG/CG, (2) single-molecule plots for mA + mCG, and (3) base abundance plots for A and CG.

# %%
# 1. Python option
# ----------------
import dimelo as dm

bam = "deep_ctcf_mod_mappings_merge.sorted.bam"
sampleName = "quartile4"
bed = "quart4.bed"
outDir = "./out"
dm.plot_enrichment_profile(
    bam,
    sampleName,
    bed,
    "A+CG",
    outDir,
    threshA=190,
    threshC=190,
    dotsize=0.05,
)


# %%
# 2. Command line option
# ----------------------
# ``dimelo-plot-enrichment-profile -f deep_ctcf_mod_mappings_merge.sorted.bam -s quartile4 -b quart4.bed -m A+CG -o ./out -A 190 -C 190 -d 0.05``

# %%
# Output
# ----------------------
# .. figure:: ../auto_examples/images/quartile4_A+CG_sm_rolling_avg.png
#     :align: center
# .. figure:: ../auto_examples/images/quartile4_A+CG_sm_scatter.png
# .. figure:: ../auto_examples/images/quartile4_A_base_count.png
# .. figure:: ../auto_examples/images/quartile4_CG_base_count.png
