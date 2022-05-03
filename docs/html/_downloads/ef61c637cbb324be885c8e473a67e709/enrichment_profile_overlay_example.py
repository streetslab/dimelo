"""
Enrichment Profile Overlay
=================================

Aggregate fraction of methylated bases centered at regions of interest defined in bed files.

"""
# %%
# Create (1) aggregate profile plots for mA/A and mCG/CG, (2) single-molecule plots for mA + mCG, and (3) base abundance plots for A and CG.

# %%
# 1. Python option
# ----------------
import dimelo as dm

bam = "deep_ctcf_mod_mappings_merge.sorted.bam"
sampleNames = ["q4", "q3", "q2", "q1"]
beds = ["quart4.bed", "quart3.bed", "quart2.bed", "quart1.bed"]
outDir = "./out"
dm.plot_enrichment_profile(
    bam, sampleNames, beds, "A", outDir, threshA=190, dotsize=0.05
)

# %%
# 2. Command line option
# ----------------------
# ``dimelo-plot-enrichment-profile -f deep_ctcf_mod_mappings_merge.sorted.bam -s q4 q3 q2 q1 -b quart4.bed quart3.bed quart2.bed quart1.bed -m A -o ./out -A 190 -d 0.05``

# %%
# Output
# ----------------------
# .. figure:: ../auto_examples/images/sample_deep_ctcf_mod_mappings_merge.sorted_A_sm_rolling_avg_overlay.png
#     :align: center
