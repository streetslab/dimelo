"""
Browser Plot Example
===========================

This plots single molecules with colored base modifications in region of interest
"""
# %%
# Plotting Browser
# ------------------------
# add description here

import dimelo as dm

input_bam = "../../dimelo/test/data/mod_mappings_subset.bam"
sample_name = "test"
region = "chr1:2907273-2909473"
outDir = "../../dimelo/dimelo_test"
dm.plot_browser(input_bam, sample_name, region, "A+CG", outDir, static=True)

# This will return an HTML file with single molecules displayed o
# ver region of interest. Modified bases are colored according to colorA and colorC.
# To return a PDF file, set static = True