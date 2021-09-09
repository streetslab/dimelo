# helper functions for dimelo_test.py

from dimelo.parse_bam import parse_bam


def create_methylation_objects():

    fileNames = [
        "dimelo/test/data/mod_mappings_subset.bam",
        "dimelo/test/data/winnowmap_guppy_merge_subset.bam",
    ]
    sampleName = "test"
    bedFile = "dimelo/test/data/test.bed"
    basemods = ["A", "CG", "A+CG"]
    centers = [True, False]
    windowSize = 1000

    all_data_combined = []
    for fileName in fileNames:
        for basemod in basemods:
            for center in centers:
                all_data = parse_bam(
                    fileName,
                    sampleName,
                    bedFile,
                    basemod,
                    center,
                    windowSize,
                )
                all_data_combined.append(all_data)

    # returned object is a list of 12 all_data objects
    # 0. mod_mappings       A       center=True
    # 1. mod_mappings       A       center=False
    # 2. mod_mappings       CG      center=True
    # 3. mod_mappings       CG      center=False
    # 4. mod_mappings       A+CG    center=True
    # 5. mod_mappings       A+CG    center=False
    # 6. winnow_guppy       A       center=True
    # 7. winnow_guppy       A       center=False
    # 8. winnow_guppy       CG      center=True
    # 9. winnow_guppy       CG      center=False
    # 10. winnow_guppy      A+CG    center=True
    # 11. winnow_guppy      A+CG    center=False
    return all_data_combined


def extract_methylation_data_subset(data, index, read_name):
    # returns a dataframe subset for given methylation data and read name
    return data[index][data[index]["read_name"] == read_name]
