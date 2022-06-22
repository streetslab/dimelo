import subprocess
import unittest
from pathlib import Path

import dimelo as dm
from dimelo.test import DiMeLoTestCase

"""
Inputs

NOTE: Changing any of these paths or database contents may require recomputing the hard-coded database
hashes in the tests below.
"""
# TODO: Is this a reasonable way to specify input files? Where is this intended to be run from?
input_bams = [
    Path("dimelo/test/data/mod_mappings_subset.bam"),
    Path("dimelo/test/data/winnowmap_guppy_merge_subset.bam"),
]
input_sample_names = ["test1", "test2"]
input_bed = Path("dimelo/test/data/test.bed")
input_region = "chr1:2907273-2909473"

"""
Outputs
"""
# TODO: When implemented elsewhere, replace these explicit database paths with modular calls
output_dbs = [f.with_suffix(".db").name for f in input_bams]


def db_hash(db_path: Path) -> bytes:
    """Computes a hash of the sqlite db at the specified path."""
    # Use the sqlite3 .sha3sum command to compute hash
    sqlite_output = subprocess.run(
        ["sqlite3", db_path, ".sha3sum"], capture_output=True
    )
    return sqlite_output.stdout


class TestParseBam(DiMeLoTestCase):
    def test_parse_bam_bedFile(self):
        """Tests parsing a bam file into a database, specifying windows using a bed file.

        Notes:
            - cores set to 1 to ensure that hashes come out the same each time
            - thresholds set super low to ensure that meaningful rows are inserted for all modifications
        """
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]
        sample_name = input_sample_names[bam_idx]
        bed_file = input_bed

        # Outputs
        db_file = output_dbs[bam_idx]

        dm.parse_bam(
            fileName=str(bam_file),
            sampleName=sample_name,
            outDir=str(self.outDir),
            bedFile=str(bed_file),
            basemod="A+CG",
            center=True,
            windowSize=500,
            threshA=1,
            threshC=1,
            extractAllBases=False,
            cores=1,
        )
        database_path = self.outDir / db_file
        # Check whether database contents are the same as expected
        self.assertEqual(
            db_hash(database_path),
            b"8a8fe3984448ee6f215d1d71c01a4d8edde7691bcdb0bf28f2c01cd2\n",
        )

    def test_parse_bam_region(self):
        """Tests parsing a bam file into a database, specifying windows using a region string.

        Notes:
            - cores set to 1 to ensure that hashes come out the same each time
            - thresholds set super low to ensure that meaningful rows are inserted for all modifications
        """
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]
        sample_name = input_sample_names[bam_idx]

        # Outputs
        db_file = output_dbs[bam_idx]

        dm.parse_bam(
            fileName=str(bam_file),
            sampleName=sample_name,
            outDir=str(self.outDir),
            region=input_region,
            basemod="A+CG",
            threshA=1,
            threshC=1,
            cores=1,
        )
        database_path = self.outDir / db_file
        # Check whether database contents are the same as expected
        self.assertEqual(
            db_hash(database_path),
            b"baf72c009547e8a2e87402db04695d698e648c3856bcc9b8c0c6cf8a\n",
        )

    def test_parse_bam_bedFile_region_mutual_exclusion(self):
        """Verifies that bedFile and region arguments remain mutually exclusive."""
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]
        sample_name = input_sample_names[bam_idx]
        bed_file = input_bed

        with self.assertRaises(RuntimeError):
            dm.parse_bam(
                fileName=str(bam_file),
                sampleName=sample_name,
                outDir=str(self.outDir),
                bedFile=str(bed_file),
                region=input_region,
            )

    def test_parse_bam_region_center_incompatible(self):
        """Verifies that region and center arguments are incompatible."""
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]
        sample_name = input_sample_names[bam_idx]

        with self.assertRaises(RuntimeError):
            dm.parse_bam(
                fileName=str(bam_file),
                sampleName=sample_name,
                outDir=str(self.outDir),
                region=input_region,
                center=True,
            )


# TODO: More robust qc_report tests
class TestQCReport(DiMeLoTestCase):
    def test_qc_report_one_sample(self):
        """Tests generating a single qc report.

        Notes:
            - cores set to 1 to ensure that hashes come out the same each time
        """
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]
        sample_name = input_sample_names[bam_idx]

        # Outputs
        db_file = output_dbs[bam_idx]
        qc_report_file = f"{sample_name}_qc_report.pdf"

        dm.qc_report(
            fileNames=str(bam_file),
            sampleNames=sample_name,
            outDir=str(self.outDir),
            cores=1,
        )

        database_path = self.outDir / db_file
        # Check whether database contents are the same as expected
        self.assertEqual(
            db_hash(database_path),
            b"58ed1ba2ce0c2e0f257ead1d9f2b9239f2777b8aefbf38ea7a40e464\n",
        )

        self.assertOutputFileExists(qc_report_file)

    def test_qc_report_multi_sample(self):
        """Tests generating multiple qc reports at once."""
        dm.qc_report(
            fileNames=[str(f) for f in input_bams],
            sampleNames=input_sample_names,
            outDir=str(self.outDir),
        )

        for db_file in output_dbs:
            self.assertOutputFileExists(db_file)
        for sample_name in input_sample_names:
            self.assertOutputFileExists(f"{sample_name}_qc_report.pdf")


class TestPlotBrowser(DiMeLoTestCase):
    def test_plot_browser_html(self):
        """Tests generation of an html browser."""
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]
        sample_name = input_sample_names[bam_idx]

        # Outputs
        db_file = output_dbs[bam_idx]

        dm.plot_browser(
            fileNames=str(bam_file),
            sampleNames=sample_name,
            region=input_region,
            basemod="A+CG",
            outDir=str(self.outDir),
            static=False,
        )

        self.assertOutputFileExists(db_file)

        # TODO: It's difficult to get the output file name. Find a better way to do this check.
        n_html_files = len(list(self.outDir.glob("*.html")))
        self.assertEqual(n_html_files, 1)

        for basemod in ["A", "C"]:
            for plot_type in ["fraction", "total"]:
                rolling_avg_file = (
                    f"{sample_name}_{basemod}_sm_rolling_avg_{plot_type}.pdf"
                )
                self.assertOutputFileExists(rolling_avg_file)

        # TODO: It's difficult to get the output file name. Find a better way to do this check.
        n_bed_files = len(list(self.outDir.glob("*.bed")))
        self.assertEqual(n_bed_files, 2)

    def test_plot_browser_pdf(self):
        """Tests generation of a pdf browser."""
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]
        sample_name = input_sample_names[bam_idx]

        # Outputs
        db_file = output_dbs[bam_idx]

        dm.plot_browser(
            fileNames=str(bam_file),
            sampleNames=sample_name,
            region=input_region,
            basemod="A+CG",
            outDir=str(self.outDir),
            static=True,
        )

        self.assertOutputFileExists(db_file)

        # TODO: It's difficult to get the output file name. Find a better way to do this check.
        n_pdf_files = len(list(self.outDir.glob("*.pdf")))
        self.assertEqual(n_pdf_files, 5)

        for basemod in ["A", "C"]:
            for plot_type in ["fraction", "total"]:
                rolling_avg_file = (
                    f"{sample_name}_{basemod}_sm_rolling_avg_{plot_type}.pdf"
                )
                self.assertOutputFileExists(rolling_avg_file)

        # TODO: It's difficult to get the output file name. Find a better way to do this check.
        n_bed_files = len(list(self.outDir.glob("*.bed")))
        self.assertEqual(n_bed_files, 2)


class TestPlotEnrichment(DiMeLoTestCase):
    def test_plot_enrichment_2_bams(self):
        """Tests enrichment comparison for the same region over two different bam files."""
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]

        # Outputs
        db_file = output_dbs[bam_idx]
        pdf_file = f"region_{input_bed.stem}_CG_enrichment_barplot.pdf"

        dm.plot_enrichment(
            fileNames=[str(bam_file), str(bam_file)],
            sampleNames=input_sample_names,
            bedFiles=str(input_bed),
            basemod="CG",
            outDir=str(self.outDir),
            threshC=129,
        )

        self.assertOutputFileExists(db_file)
        self.assertOutputFileExists(pdf_file)

    def test_plot_enrichment_2_beds(self):
        """Tests enrichment comparison for two different regions over the same bam file"""
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]

        # Outputs
        db_file = output_dbs[bam_idx]
        pdf_file = f"sample_{bam_file.stem}_CG_enrichment_barplot.pdf"

        dm.plot_enrichment(
            fileNames=str(bam_file),
            sampleNames=input_sample_names,
            bedFiles=[str(input_bed), str(input_bed)],
            basemod="CG",
            outDir=str(self.outDir),
            threshC=129,
        )

        self.assertOutputFileExists(db_file)
        self.assertOutputFileExists(pdf_file)

    def test_plot_enrichment_incompatible_bed_bam(self):
        """Verifies that passing equal numbers of beds and bams remains an error."""
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]

        with self.assertRaises(RuntimeError):
            dm.plot_enrichment(
                fileNames=[str(bam_file), str(bam_file)],
                sampleNames=input_sample_names,
                bedFiles=[str(input_bed), str(input_bed)],
                basemod="CG",
                outDir=str(self.outDir),
            )

    def test_plot_enrichment_incompatible_basemod(self):
        """Verifies that plot_enrichment remains incompatible with multi-basemod options."""
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]

        with self.assertRaises(RuntimeError):
            dm.plot_enrichment(
                fileNames=[str(bam_file), str(bam_file)],
                sampleNames=input_sample_names,
                bedFiles=str(input_bed),
                basemod="A+CG",
                outDir=str(self.outDir),
            )


class TestPlotEnrichmentProfile(DiMeLoTestCase):
    def test_plot_enrichment_profile_single(self):
        """Tests profile plotting for a single sample and region."""
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]
        sample_name = input_sample_names[bam_idx]

        # Outputs
        db_file = output_dbs[bam_idx]
        enrichment_plot = f"{sample_name}_A+CG_sm_rolling_avg.pdf"
        single_molecule_plot = f"{sample_name}_A+CG_sm_scatter.png"
        base_count_plots = [
            f"{sample_name}_{basemod}_base_count.png"
            for basemod in ["A", "CG"]
        ]

        dm.plot_enrichment_profile(
            fileNames=str(bam_file),
            sampleNames=sample_name,
            bedFiles=str(input_bed),
            basemod="A+CG",
            outDir=str(self.outDir),
            windowSize=500,
            dotsize=1,
        )

        self.assertOutputFileExists(db_file)
        self.assertOutputFileExists(enrichment_plot)
        self.assertOutputFileExists(single_molecule_plot)
        for f in base_count_plots:
            self.assertOutputFileExists(f)

    def test_plot_enrichment_profile_sample_overlay(self):
        """Tests profile plotting for multiple regions over a single sample."""
        bam_idx = 0

        # Inputs
        bam_file = input_bams[bam_idx]

        # Outputs
        db_file = output_dbs[bam_idx]
        overlay_plot = f"sample_{bam_file.stem}_A_sm_rolling_avg_overlay.pdf"

        dm.plot_enrichment_profile(
            fileNames=str(bam_file),
            sampleNames=input_sample_names,
            bedFiles=[str(input_bed), str(input_bed)],
            basemod="A",
            outDir=str(self.outDir),
            windowSize=500,
            dotsize=1,
        )

        self.assertOutputFileExists(db_file)
        self.assertOutputFileExists(overlay_plot)

    def test_plot_enrichment_profile_region_overlay(self):
        """Tests profile plotting for multiple samples over a single region."""
        # Outputs
        overlay_plot = f"region_{input_bed.stem}_A_sm_rolling_avg_overlay.pdf"

        dm.plot_enrichment_profile(
            fileNames=[str(f) for f in input_bams],
            sampleNames=input_sample_names,
            bedFiles=str(input_bed),
            basemod="A",
            outDir=str(self.outDir),
            windowSize=500,
            dotsize=1,
        )

        for db_file in output_dbs:
            self.assertOutputFileExists(db_file)
        self.assertOutputFileExists(overlay_plot)

    def test_plot_enrichment_profile_overlay_incompatible_basemod(self):
        """Verifies that plot_enrichment_profile remains incompatible with multi-basemod options."""
        with self.assertRaises(RuntimeError):
            dm.plot_enrichment_profile(
                fileNames=[str(f) for f in input_bams],
                sampleNames=input_sample_names,
                bedFiles=str(input_bed),
                basemod="A+CG",
                outDir=str(self.outDir),
            )

    def test_plot_enrichment_profile_overlay_incompatible_bed_bam(self):
        """Verifies that passing equal numbers of beds and bams remains an error."""
        with self.assertRaises(RuntimeError):
            dm.plot_enrichment_profile(
                fileNames=[str(f) for f in input_bams],
                sampleNames=input_sample_names,
                bedFiles=[str(input_bed), str(input_bed)],
                basemod="A",
                outDir=str(self.outDir),
            )


if __name__ == "__main__":
    unittest.main()
