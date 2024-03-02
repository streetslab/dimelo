import subprocess
import unittest
from pathlib import Path
import urllib.request
import gzip
import filecmp

import numpy as np
from matplotlib.axes import Axes
import pytest

import dimelo as dm
from test import DiMeLoParsingTestCase

script_location = Path(__file__).resolve().parent

data_dir = script_location / 'data'
peaks_bed = data_dir / 'ctcf_demo_peak.bed'
not_peaks_bed = data_dir / 'ctcf_demo_not_peak.bed'
region = 'chr14:17376348-17378348'

targets_dir = data_dir / 'test_targets'
pileup_merged_bedgz_target = targets_dir / 'megalodon_merged_regions' / 'pileup.sorted.bed.gz'
pileup_merged_bedregions_target = targets_dir / 'megalodon_merged_regions' / 'regions.processed.bed'
extract_merged_h5_target = None
pileup_one_bedgz_target = targets_dir / 'megalodon_one_region' / 'pileup.sorted.bed.gz'
pileup_one_bedregions_target = targets_dir / 'megalodon_one_region' / 'regions.processed.bed'
extract_one_h5_target = targets_dir / 'megalodon_one_region' / 'reads.combined_basemods.h5'

@pytest.mark.parametrize("output_name,regions,motifs,thresh,window_size,pileup_target,regions_target,extract_target", [
    ('megalodon_merged_regions',
     [peaks_bed,not_peaks_bed],
     ['A,0','CG,0'],
     190,
     1000,
     pileup_merged_bedgz_target,
     pileup_merged_bedregions_target,
     extract_merged_h5_target), 
    ('megalodon_one_region',
     region,
     ['A,0','CG,0'],
     190,
     None,
     pileup_one_bedgz_target,
     pileup_one_bedregions_target,
     extract_one_h5_target),
    ])
class TestParseBam(DiMeLoParsingTestCase):
    """
    Tests parsing a bam file into a bed.gz pileup and an hdf5 single read file.
    
    This test class requires the output files be bitwise identical, compared to pre-defined reference files.
    This means that interface changes require replacing these files.
    """
    
    def test_pileup(
        cls,
        output_name,
        regions,
        motifs,
        thresh,
        window_size,
        pileup_target,
        extract_target,
        regions_target,
    ):
        pileup_bed,regions_processed = dm.parse_bam.pileup(
            input_file = cls.megalodon_input,
            output_name = output_name,
            ref_genome = cls.reference_genome,
            output_directory = cls.outDir,
            regions = regions,
            motifs = motifs,
            thresh = thresh,
            window_size = window_size,
        )

        if pileup_target is not None:
            assert filecmp.cmp(pileup_bed,pileup_target,shallow=False)
        if regions_target is not None:
            assert filecmp.cmp(regions_processed,regions_target,shallow=False)

    def test_extract(
        cls,
        output_name,
        regions,
        motifs,
        thresh,
        window_size,
        pileup_target,
        extract_target,
        regions_target,           
    ):
        extract_h5,regions_processed = dm.parse_bam.extract(
            input_file = cls.megalodon_input,
            output_name = output_name,
            ref_genome = cls.reference_genome,
            output_directory = cls.outDir,
            regions = regions,
            motifs = motifs,
            thresh = thresh,
            window_size = window_size,
        )

        if extract_target is not None:
            assert filecmp.cmp(extract_h5,extract_target,shallow=False)
        if regions_target is not None:
            assert filecmp.cmp(regions_processed,regions_target,shallow=False)

    


class TestLoadProcessed:
    """
    Tests loading values from bed.gz pileups and hdf5 single read files.
    
    This test class requires that values are identical. It loads from pre-defined reference files.
    This means that interface changes require replacing these files.
    """
    def test_pileup_counts_from_bedmethyl_A0_peaks(self):
        modified_base_count,valid_base_count = dm.load_processed.pileup_counts_from_bedmethyl(
            bedmethyl_file = pileup_merged_bedgz_target,
            motif = 'A,0',
            regions = peaks_bed)
        
        assert modified_base_count == 14242 and valid_base_count == 174078

    def test_pileup_counts_from_bedmethyl_A0_site(self):
        modified_base_count,valid_base_count = dm.load_processed.pileup_counts_from_bedmethyl(
            bedmethyl_file = pileup_merged_bedgz_target,
            motif = 'A,0',
            regions = 'chr1:114358436-114358754')
        
        assert modified_base_count == 213 and valid_base_count == 1619

    def test_pileup_counts_from_bedmethyl_CG0_peaks(self):
        modified_base_count,valid_base_count = dm.load_processed.pileup_counts_from_bedmethyl(
            bedmethyl_file = pileup_merged_bedgz_target,
            motif = 'CG,0',
            regions = peaks_bed)
        
        assert modified_base_count == 66 and valid_base_count == 34790

    def test_pileup_vectors_from_bedmethyl_A0_peaks(self):
        modified_base_counts,valid_base_counts = dm.load_processed.pileup_vectors_from_bedmethyl(
            bedmethyl_file = pileup_merged_bedgz_target,
            motif = 'A,0',
            regions = peaks_bed, # Regions from bed file
            window_size = 10, # Trim/extend regions to same size
        )       
        mod_correct_str = '38 53 32  6  4 54  2 83 90  0 85 39 15 45 11 18 11 16 38 20'
        val_correct_str = '588 679 297 210  70 534  58 445 839   0 891 431  80 486  51 303 286 514 491 589'
        mod_array = np.array([int(n) for n in mod_correct_str.split()])
        val_array = np.array([int(n) for n in val_correct_str.split()])
        assert np.array_equal(modified_base_counts,mod_array) and np.array_equal(valid_base_counts,val_array)


class TestPlotEnrichment:
    """
    Tests plotting functionality in plot_enrichment.

    This test simply checks that we can make plots from real or synthetic data without raising errors. 
    Appearance of plots is not verified.
    """
    def test_plot_enrichment_plot_enrichment_synthetic(self):
        ax = dm.plot_enrichment.plot_enrichment(
            mod_file_names=['test.fake', 'test.fake'],
            regions_list=['test.bed', 'test.bed'],
            motifs=['A', 'A'],
            sample_names=['a', 'b']
        )
        assert isinstance(ax,Axes)

    def test_plot_enrichment_by_modification_synthetic(self):
        ax = dm.plot_enrichment.by_modification(
            mod_file_name='test.fake',
            regions='test.bed',
            motifs=['A', 'C']
        )
        assert isinstance(ax,Axes)

    def test_plot_enrichment_by_regions_synthetic(self):
        ax = dm.plot_enrichment.by_regions(
            mod_file_name='test.fake',
            regions_list=['test1.bed', 'test2.bed'],
            motif='A'
        )
        assert isinstance(ax,Axes)

    def test_plot_enrichment_by_dataset_synthetic(self):
        ax = dm.plot_enrichment.by_dataset(
            mod_file_names=['test1.fake', 'test2.fake'],
            regions='test.bed',
            motif='A'
        )
        assert isinstance(ax,Axes)

    def test_plot_enrichment_by_regions_megalodon(self):
        ax = dm.plot_enrichment.by_regions(
            mod_file_name=pileup_merged_bedgz_target,
            regions_list=[peaks_bed, not_peaks_bed],
            motif='A,0',
            sample_names=['on target', 'off target']
        )
        assert isinstance(ax,Axes)

    def test_plot_enrichment_by_modification_megalodon(self):
        ax = dm.plot_enrichment.by_modification(
            mod_file_name=pileup_merged_bedgz_target,
            regions=peaks_bed,
            motifs=['A,0','CG,0'],)
        assert isinstance(ax,Axes)

class TestPlotEnrichmentProfile:
    """
    Tests plotting functionality in plot_enrichment_profile.

    This test simply checks that we can make plots from real or synthetic data without raising errors. 
    Appearance of plots is not verified.
    """
    def test_plot_enrichment_profile_plot_enrichment_profile_synthetic(self):
        ax = dm.plot_enrichment_profile.plot_enrichment_profile(
            mod_file_names=['test1.fake', 'test2.fake'],
            regions_list=['test1.bed', 'test2.bed'],
            motifs=['A', 'C'],
            window_size=500,
            sample_names=['sample1', 'sample2'],
            smooth_window=50
        )
        assert isinstance(ax,Axes)

    def test_plot_enrichment_profile_by_modification_synthetic(self):
        ax = dm.plot_enrichment_profile.by_modification(
            mod_file_name='test.fake',
            regions='test.bed',
            window_size=500,
            motifs=['A', 'C'],
            smooth_window=50
        )
        assert isinstance(ax,Axes)

    def test_plot_enrichment_profile_by_region_synthetic(self):
        ax = dm.plot_enrichment_profile.by_regions(
            mod_file_name='test.fake',
            regions_list=['test1.bed', 'test2.bed'],
            motif='A',
            window_size=500,
            sample_names=['on target', 'off target'],
            smooth_window=50
        )
        assert isinstance(ax,Axes)

    def test_plot_enrichment_profile_by_dataset_synthetic(self):
        ax = dm.plot_enrichment_profile.by_dataset(
            mod_file_names=['test1.fake', 'test2.fake'],
            regions='test.bed',
            motif='A',
            window_size=500,
            sample_names=['experiment 1', 'experiment 2'],
            smooth_window=50
        )
        assert isinstance(ax,Axes)

    def test_plot_enrichment_profile_by_modification_megalodon(self):
        ax = dm.plot_enrichment_profile.by_modification(
            mod_file_name=pileup_merged_bedgz_target,
            regions=peaks_bed,
            window_size=1000,
            motifs=['A,0','CG,0'],
            smooth_window=50
        )
        assert isinstance(ax,Axes)

    def test_plot_enrichment_profile_by_regions_megalodon(self):
        ax = dm.plot_enrichment_profile.by_regions(
            mod_file_name=pileup_merged_bedgz_target,
            regions_list=[peaks_bed,not_peaks_bed],
            window_size=1000,
            motif='A,0',
            smooth_window=50           
        )

class TestPlotReads:
    """
    Tests plotting functionality in plot_reads.

    This test simply checks that we can make plots from real or synthetic data without raising errors. 
    Appearance of plots is not verified.
    """
    def test_plot_reads_plot_reads_synthetic(self):
        ax = dm.plot_reads.plot_reads(
            mod_file_name = 'test.fake', 
            regions = 'test.bed', 
            motifs = ['A,0', 'CG,0']
        )
        assert isinstance(ax,Axes)

    def test_plot_reads_plot_reads_megalodon(self):
        ax = dm.plot_reads.plot_reads(
            mod_file_name = extract_one_h5_target,
            regions = pileup_one_bedregions_target,
            motifs = ['A,0', 'CG,0'],
        )
        assert isinstance(ax,Axes)

class TestParseToPlot(DiMeLoParsingTestCase):
    """
    Tests that each stage of parse_bam -> load_processed -> plotting works correctly, including values
    where applicable.

    This test does not look at the saved pileup/reads files themselves, simply the values that get pulled out.
    Thus it tests interfaces end-to-end. Test coverage will for the time being be less comprehensive than the 
    unit tests.
    """
    def test_pileup_load_plot(cls):
        pileup_bed,regions_processed = dm.parse_bam.pileup(
            input_file = cls.megalodon_input,
            output_name = 'megalodon_merged_regions',
            ref_genome = cls.reference_genome,
            output_directory = cls.outDir,
            regions = [peaks_bed,not_peaks_bed],
            motifs = ['A,0', 'CG,0'],
            thresh = 190,
            window_size = 1000,
        )
        modified_base_count,valid_base_count = dm.load_processed.pileup_counts_from_bedmethyl(
            bedmethyl_file = pileup_bed,
            motif = 'A,0',
            regions = peaks_bed)
        
        assert modified_base_count == 14242 and valid_base_count == 174078

        ax = dm.plot_enrichment.by_regions(
            mod_file_name=pileup_bed,
            regions_list=[peaks_bed, not_peaks_bed],
            motif='A,0',
            sample_names=['on target', 'off target']
        )
        assert isinstance(ax,Axes)


    def test_extract(cls):
        extract_h5,regions_processed = dm.parse_bam.extract(
            input_file = cls.megalodon_input,
            output_name = 'megalodon_one_region',
            ref_genome = cls.reference_genome,
            output_directory = cls.outDir,
            regions = region,
            motifs = ['A,0', 'CG,0'],
            thresh = 190,
            window_size = None,
        )

        ax = dm.plot_reads.plot_reads(
            mod_file_name=extract_h5,
            regions=region,
            motifs=['A,0', 'CG,0'],   
        )
        assert isinstance(ax,Axes)


    