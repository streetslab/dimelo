import filecmp
import gzip
import pickle
from pathlib import Path

import h5py
import numpy as np
import plotly
import pytest
from matplotlib.axes import Axes

import dimelo as dm
from dimelo.test import DiMeLoParsingTestCase, filter_kwargs_for_func

script_location = Path(__file__).parent

with open(
    script_location / "data" / "test_targets" / "test_matrix.pickle", "rb"
) as file:
    test_matrix = pickle.load(file)


@pytest.mark.parametrize(
    "test_case,kwargs,results",
    [(case, inputs, outputs) for case, (inputs, outputs) in test_matrix.items()],
)
class TestParseToPlot(DiMeLoParsingTestCase):
    """
    Tests parsing a bam file into a bed.gz pileup and an hdf5 single read file, then tests that each stage
    of parse_bam -> load_processed -> plotting works correctly, including comparing where applicable.

    This test class requires the output files be bitwise identical, compared to pre-defined reference files.
    This means that interface changes require replacing these files.

    For integration tests we test interfaces end-to-end.
    """

    def test_unit__pileup(
        cls,
        test_case,
        kwargs,
        results,
    ):
        kwargs_pileup = filter_kwargs_for_func(dm.parse_bam.pileup, kwargs)
        kwargs_pileup["output_directory"] = cls.outDir
        pileup_bed, regions_processed = dm.parse_bam.pileup(
            **kwargs_pileup,
            ref_genome=cls.reference_genome,
        )

        pileup_target, regions_target = results["pileup"]

        if pileup_target is not None and regions_target is not None:
            # This is necessary because the gzipped files are different on mac vs linux, but the contents should be identical (and are, so far)
            # Not sure why the compression ratio is better on Linux when both are using pysam.tabix_compress with pysam 0.22.0 and zlib 1.2.13 but whatcha gonna do
            with (
                gzip.open(pileup_bed, "rt") as f1,
                gzip.open(pileup_target, "rt") as f2,
            ):
                # Read and compare file contents
                file1_contents = f1.read()
                file2_contents = f2.read()
                assert (
                    file1_contents == file2_contents
                ), f"{test_case}: {pileup_bed} does not match {pileup_target}."
            assert filecmp.cmp(
                regions_processed, regions_target, shallow=False
            ), f"{test_case}: {regions_processed} does not match {regions_target}."
        else:
            print(f"{test_case} skipped for pileup.")

    def test_unit__extract(
        cls,
        test_case,
        kwargs,
        results,
    ):
        kwargs_extract = filter_kwargs_for_func(dm.parse_bam.extract, kwargs)
        kwargs_extract["output_directory"] = cls.outDir
        # extract can behave non-deterministically in terms of read output order
        # if run on multiple cores. This is not an issue for sorted read loads,
        # but means the file itself changes even though the meaningful contents don't
        # given that parallelization here is purely within modkit anyway, the design
        # choice was made to always single-core for extract testing
        if "cores" in kwargs_extract:
            del kwargs_extract["cores"]
        extract_h5, regions_processed = dm.parse_bam.extract(
            **kwargs_extract,
            ref_genome=cls.reference_genome,
            cores=1,
        )

        extract_target, regions_target = results["extract"]

        if extract_target is not None and regions_target is not None:
            # The hdf5 files will have a few bits different due to gzip compression timestamps, but comparing the exact size should pass because
            # the timestamps are not themselves compressed inside the vector gzip objects
            h5_test = h5py.File(extract_h5)
            h5_target = h5py.File(extract_target)
            datasets = [
                name for name, obj in h5_target.items() if isinstance(obj, h5py.Dataset)
            ]
            for dataset in datasets:
                if dataset in ["threshold"]:
                    assert h5_test[dataset] == h5_target[dataset] or np.isnan(
                        h5_test[dataset]
                    ) == np.isnan(h5_target[dataset])
                else:
                    test_dataset = list(h5_test[dataset][:])
                    target_dataset = list(h5_target[dataset][:])
                    if dataset in ["mod_vector", "val_vector"]:
                        assert [
                            gzip.decompress(test_item.tobytes())
                            for test_item in test_dataset
                        ] == [
                            gzip.decompress(target_item.tobytes())
                            for target_item in target_dataset
                        ], f"{test_case}: {dataset} does not match."
                    else:
                        assert (
                            test_dataset == target_dataset
                        ), f"{test_case}: {dataset} does not match."
            # assert os.path.getsize(extract_h5) == os.path.getsize(extract_target), f"{test_case}: {extract_h5} does not match {extract_target}."
            assert filecmp.cmp(
                regions_processed, regions_target, shallow=False
            ), f"{test_case}: {regions_processed} does not match {regions_target}."
        else:
            print(f"{test_case} skipped for extract.")

    def test_integration__pileup_load_plot(
        cls,
        test_case,
        kwargs,
        results,
    ):
        # This stuff is commented out because if we run this integration test in the same class as the unit test
        # for parsing, we can cut down total end-to-end testing overhead by about 2x by just using that output
        #
        # kwargs_pileup = filter_kwargs_for_func(dm.parse_bam.pileup,kwargs)
        # kwargs_pileup['output_directory'] = cls.outDir
        # pileup_bed,_ = dm.parse_bam.pileup(
        #     **kwargs_pileup,
        #     ref_genome = cls.reference_genome,
        # )
        # We just grab the output from TestParseBam::test_pileup, wasteful to re-run an identical modkit command
        pileup_bed = cls.outDir / kwargs["output_name"] / "pileup.sorted.bed.gz"

        # If we have results for this pileup, check that the load_processed values are ok out of the output file
        if results["pileup"][0] is not None:
            kwargs_counts_from_bedmethyl = filter_kwargs_for_func(
                dm.load_processed.pileup_counts_from_bedmethyl, kwargs
            )
            for motif in kwargs["motifs"]:
                expected = results["pileup_counts_from_bedmethyl"][motif]
                actual = dm.load_processed.pileup_counts_from_bedmethyl(
                    bedmethyl_file=pileup_bed,
                    motif=motif,
                    **kwargs_counts_from_bedmethyl,
                )
                assert (
                    actual == expected
                ), f"{test_case}: Counts for motif {motif} are not equal"

            kwargs_vectors_from_bedmethyl = filter_kwargs_for_func(
                dm.load_processed.pileup_vectors_from_bedmethyl, kwargs
            )
            for motif in kwargs["motifs"]:
                expected_tuple = results["pileup_vectors_from_bedmethyl"][motif]
                actual_tuple = dm.load_processed.pileup_vectors_from_bedmethyl(
                    bedmethyl_file=results["pileup"][0],
                    motif=motif,
                    **kwargs_vectors_from_bedmethyl,
                )
                assert len(expected_tuple) == len(
                    actual_tuple
                ), f"{test_case}: Unexpected number of arrays returned for {motif}"

                for expected, actual in zip(expected_tuple, actual_tuple):
                    # TODO: The following was the original assertion error message, but it was not written in a functional way. Find a way to make it work as intended.
                    # assert np.array_equal(expected, actual), f"{test_case}: Arrays for motif {motif} are not equal: expected {value} but got {actual[key]}"
                    assert np.array_equal(
                        expected, actual
                    ), f"{test_case}: Arrays for motif {motif} are not equal."
        else:
            print(
                f"{test_case} loading skipped for pileup_load_plot, continuing to plotting."
            )

        kwargs_plot_enrichment_plot_enrichment = filter_kwargs_for_func(
            dm.plot_enrichment.plot_enrichment, kwargs
        )
        for motif in kwargs["motifs"]:
            regions_list = (
                kwargs["regions"]
                if isinstance(kwargs["regions"], list)
                else [kwargs["regions"]]
            )
            kwargs_plot_enrichment_plot_enrichment["motifs"] = [
                motif for _ in regions_list
            ]
            ax = dm.plot_enrichment.plot_enrichment(
                mod_file_names=[pileup_bed for _ in regions_list],
                regions_list=regions_list,
                sample_names=["label" for _ in regions_list],
                **kwargs_plot_enrichment_plot_enrichment,
            )
            assert isinstance(ax, Axes), f"{test_case}: plotting failed for {motif}."
        kwargs_plot_enrichment_profile_plot_enrichment_profile = filter_kwargs_for_func(
            dm.plot_enrichment_profile.plot_enrichment_profile,
            kwargs,
            extra_args=["window_size", "smooth_window"],
        )
        for motif in kwargs["motifs"]:
            regions_list = (
                kwargs["regions"]
                if isinstance(kwargs["regions"], list)
                else [kwargs["regions"]]
            )
            kwargs_plot_enrichment_profile_plot_enrichment_profile["motifs"] = [
                motif for _ in regions_list
            ]
            ax = dm.plot_enrichment_profile.plot_enrichment_profile(
                mod_file_names=[pileup_bed for _ in regions_list],
                regions_list=regions_list,
                sample_names=["label" for _ in regions_list],
                **kwargs_plot_enrichment_profile_plot_enrichment_profile,
            )
            assert isinstance(ax, Axes), f"{test_case}: plotting failed for {motif}."

    def test_integration__extract_load_plot(
        cls,
        test_case,
        kwargs,
        results,
    ):
        # This stuff is commented out because if we run this integration test in the same class as the unit test
        # for parsing, we can cut down total end-to-end testing overhead by about 2x by just using that output
        #
        # if results['extract'][0] is None:
        #     return

        # kwargs_extract = filter_kwargs_for_func(dm.parse_bam.extract,kwargs)
        # kwargs_extract['output_directory'] = cls.outDir
        # extract_h5,_ = dm.parse_bam.extract(
        #     **kwargs_extract,
        #     ref_genome = cls.reference_genome,
        # )
        # We just grab the output from TestParseBam::test_extract, wasteful to re-run an identical modkit command
        extract_h5 = cls.outDir / kwargs["output_name"] / "reads.combined_basemods.h5"

        # If we have results for this extraction, check that the load_processed values are ok out of the output file
        if results["extract"][0] is not None:
            kwargs_read_vectors_from_hdf5 = filter_kwargs_for_func(
                dm.load_processed.read_vectors_from_hdf5, kwargs
            )
            read_data_list, datasets, _ = dm.load_processed.read_vectors_from_hdf5(
                file=extract_h5,
                **kwargs_read_vectors_from_hdf5,
            )
            read_data_dict = {}
            # Pull out the data from the first read
            for idx, dataset in enumerate(datasets):
                for read_data in read_data_list:
                    read_data_dict[dataset] = read_data[idx]
                    break
            expected = results["read_vectors_from_hdf5"]
            actual = read_data_dict
            for key, value in expected.items():
                if isinstance(value, np.ndarray):
                    assert np.allclose(
                        actual[key], expected[key], atol=1e-5
                    ), f"""{test_case}: Arrays for {key} are not equal
mismatch at {np.where(value != actual[key])}
mismatch values expected {value[np.where(value != actual[key])]} vs actual {actual[key][np.where(value != actual[key])]}
{value[np.where(value != actual[key])[0]]} vs {actual[key][np.where(value != actual[key])[0]]}.
                    """
                elif isinstance(value, (str, int, bool)):
                    assert (
                        actual[key] == expected[key]
                    ), f"{test_case}: Values for {key} are not equal: expected {value} but got {actual[key]}."
                else:
                    assert np.isclose(
                        actual[key], value, atol=1e-4
                    ), f"{test_case}: Values for {key} are not equal: expected {value} but got {actual[key]}."
        else:
            print("{test_case} skipped for read_vectors_from_hdf5.")
        kwargs_plot_reads_plot_reads = filter_kwargs_for_func(
            dm.plot_reads.plot_reads, kwargs
        )
        if kwargs["thresh"] is not None:
            ax = dm.plot_reads.plot_reads(
                mod_file_name=extract_h5,
                **kwargs_plot_reads_plot_reads,
            )
            assert isinstance(ax, Axes), f"{test_case}: plotting failed."
        else:  # if the extract parameters did not have a threshold, plot_reads.plot_reads should raise an error
            with pytest.raises(ValueError) as excinfo:
                ax = dm.plot_reads.plot_reads(
                    mod_file_name=extract_h5,
                    **kwargs_plot_reads_plot_reads,
                )
            assert "No threshold has been applied" in str(
                excinfo.value
            ), f"{test_case}: unexpected exception {excinfo.value}"
            # providing a threshold should be enough to run plot_reads.plot_reads without an error
            kwargs_plot_reads_plot_reads["thresh"] = 0.75
            ax = dm.plot_reads.plot_reads(
                mod_file_name=extract_h5,
                **kwargs_plot_reads_plot_reads,
            )
            assert isinstance(ax, Axes), f"{test_case}: plotting failed."


@pytest.mark.parametrize(
    "test_case,kwargs,results",
    [(case, inputs, outputs) for case, (inputs, outputs) in test_matrix.items()],
)
class TestLoadProcessed:
    """
    Tests loading values from bed.gz pileups and hdf5 single read files.

    This test class requires that values are identical. It loads from pre-defined reference files.
    This means that interface changes require replacing these files by re-running the Generate parse_bam
    outputs section of dimelo/test/generate_test_targets.ipynb.
    """

    def test_unit__regions_to_list(
        self,
        test_case,
        kwargs,
        results,
    ):
        """
        This test currently only tests that regions_to_list can run all relevant loaders, and assumes their
        values are correct based on the subsequent tests that verify values.
        """
        if results["pileup"][0] is not None:
            # test pileup loading
            kwargs_counts_from_bedmethyl = filter_kwargs_for_func(
                dm.load_processed.pileup_counts_from_bedmethyl, kwargs
            )
            kwargs_vectors_from_bedmethyl = filter_kwargs_for_func(
                dm.load_processed.pileup_vectors_from_bedmethyl, kwargs
            )
            for motif in kwargs["motifs"]:
                dm.load_processed.regions_to_list(
                    function_handle=dm.load_processed.pileup_counts_from_bedmethyl,
                    bedmethyl_file=results["pileup"][0],
                    motif=motif,
                    **kwargs_counts_from_bedmethyl,
                )
                dm.load_processed.regions_to_list(
                    function_handle=dm.load_processed.pileup_vectors_from_bedmethyl,
                    bedmethyl_file=results["pileup"][0],
                    motif=motif,
                    **kwargs_vectors_from_bedmethyl,
                )
        if results["extract"][0] is not None:
            kwargs_read_vectors_from_hdf5 = filter_kwargs_for_func(
                dm.load_processed.read_vectors_from_hdf5, kwargs
            )
            dm.load_processed.regions_to_list(
                function_handle=dm.load_processed.read_vectors_from_hdf5,
                file=results["extract"][0],
                **kwargs_read_vectors_from_hdf5,
            )

    def test_unit__pileup_counts_from_bedmethyl(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["pileup"][0] is not None:
            kwargs_counts_from_bedmethyl = filter_kwargs_for_func(
                dm.load_processed.pileup_counts_from_bedmethyl, kwargs
            )
            for motif in kwargs["motifs"]:
                expected = results["pileup_counts_from_bedmethyl"][motif]
                actual = dm.load_processed.pileup_counts_from_bedmethyl(
                    bedmethyl_file=results["pileup"][0],
                    motif=motif,
                    **kwargs_counts_from_bedmethyl,
                )
                assert (
                    actual == expected
                ), f"{test_case}: Counts for motif {motif} are not equal"
        else:
            print(f"{test_case} skipped for pileup_counts_from_bedmethyl.")

    def test_unit__pileup_vectors_from_bedmethyl(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["pileup"][0] is not None:
            kwargs_vectors_from_bedmethyl = filter_kwargs_for_func(
                dm.load_processed.pileup_vectors_from_bedmethyl, kwargs
            )
            for motif in kwargs["motifs"]:
                expected_tuple = results["pileup_vectors_from_bedmethyl"][motif]
                actual_tuple = dm.load_processed.pileup_vectors_from_bedmethyl(
                    bedmethyl_file=results["pileup"][0],
                    motif=motif,
                    **kwargs_vectors_from_bedmethyl,
                )
                assert len(expected_tuple) == len(
                    actual_tuple
                ), f"{test_case}: Unexpected number of arrays returned for {motif}"

                for expected, actual in zip(expected_tuple, actual_tuple):
                    assert np.array_equal(
                        expected, actual
                    ), f"{test_case}: Arrays for motif {motif} are not equal"
        else:
            print(f"{test_case} skipped for pileup_vectors_from_bedmethyl.")

    def test_unit__read_vectors_from_hdf5(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["extract"][0] is not None:
            kwargs_read_vectors_from_hdf5 = filter_kwargs_for_func(
                dm.load_processed.read_vectors_from_hdf5, kwargs
            )
            read_data_list, datasets, _ = dm.load_processed.read_vectors_from_hdf5(
                file=results["extract"][0],
                **kwargs_read_vectors_from_hdf5,
            )
            read_data_dict = {}
            # Pull out the data from the first read
            for idx, dataset in enumerate(datasets):
                for read_data in read_data_list:
                    read_data_dict[dataset] = read_data[idx]
                    break
            expected = results["read_vectors_from_hdf5"]
            actual = read_data_dict
            for key, value in expected.items():
                if isinstance(value, np.ndarray):
                    assert np.allclose(
                        actual[key], expected[key], atol=1e-5
                    ), f"""{test_case}: Arrays for {key} are not equal
mismatch at {np.where(value != actual[key])}
mismatch values expected {value[np.where(value != actual[key])]} vs actual {actual[key][np.where(value != actual[key])]}
{value[np.where(value != actual[key])[0]]} vs {actual[key][np.where(value != actual[key])[0]]}.
                    """
                elif isinstance(value, (str, int, bool)):
                    assert (
                        actual[key] == expected[key]
                    ), f"{test_case}: Values for {key} are not equal: expected {value} but got {actual[key]}."
                else:
                    assert np.isclose(
                        actual[key], value, atol=1e-4
                    ), f"{test_case}: Values for {key} are not equal: expected {value} but got {actual[key]}."
        else:
            print("{test_case} skipped for read_vectors_from_hdf5.")


@pytest.mark.parametrize(
    "test_case,kwargs,results",
    [(case, inputs, outputs) for case, (inputs, outputs) in test_matrix.items()],
)
class TestExport(DiMeLoParsingTestCase):
    """
    Tests file export functionality in export module.

    This test currently simply checks that we can make the appropriate output files without raising errors.
    The values stored in the files are not verified. Future work should add test coverage for values, but at
    the moment there is no loading infrastructure in place for bigwig files making such implementation high-overhead.
    """

    def test_unit__pileup_to_bigwig(
        cls,
        test_case,
        kwargs,
        results,
    ):
        kwargs_bigwig = filter_kwargs_for_func(dm.export.pileup_to_bigwig, kwargs)
        kwargs_bigwig["bedmethyl_file"] = results["pileup"][0]
        kwargs_bigwig["bigwig_file"] = (
            cls.outDir / kwargs["output_name"] / "pileup.fractions.bigwig"
        )
        for motif in kwargs["motifs"]:
            dm.export.pileup_to_bigwig(
                **kwargs_bigwig,
                motif=motif,
            )


class TestPlotEnrichmentSynthetic:
    """
    Tests plotting functionality in plot_enrichment.

    This test simply checks that we can make plots from synthetic data without raising errors.
    Appearance of plots is not verified.
    """

    def test_unit__plot_enrichment_plot_enrichment_synthetic(self):
        ax = dm.plot_enrichment.plot_enrichment(
            mod_file_names=["test.fake", "test.fake"],
            regions_list=["test.bed", "test.bed"],
            motifs=["A", "A"],
            sample_names=["a", "b"],
        )
        assert isinstance(ax, Axes)

    def test_unit__plot_enrichment_by_modification_synthetic(self):
        ax = dm.plot_enrichment.by_modification(
            mod_file_name="test.fake", regions="test.bed", motifs=["A", "C"]
        )
        assert isinstance(ax, Axes)

    def test_unit__plot_enrichment_by_regions_synthetic(self):
        ax = dm.plot_enrichment.by_regions(
            mod_file_name="test.fake",
            regions_list=["test1.bed", "test2.bed"],
            motif="A",
        )
        assert isinstance(ax, Axes)

    def test_unit__plot_enrichment_by_dataset_synthetic(self):
        ax = dm.plot_enrichment.by_dataset(
            mod_file_names=["test1.fake", "test2.fake"], regions="test.bed", motif="A"
        )
        assert isinstance(ax, Axes)


@pytest.mark.parametrize(
    "test_case,kwargs,results",
    [(case, inputs, outputs) for case, (inputs, outputs) in test_matrix.items()],
)
class TestPlotEnrichment:
    def test_unit__plot_enrichment_plot_enrichment(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["pileup"][0] is not None:
            kwargs_plot_enrichment_plot_enrichment = filter_kwargs_for_func(
                dm.plot_enrichment.plot_enrichment, kwargs
            )
            for motif in kwargs["motifs"]:
                regions_list = (
                    kwargs["regions"]
                    if isinstance(kwargs["regions"], list)
                    else [kwargs["regions"]]
                )
                kwargs_plot_enrichment_plot_enrichment["motifs"] = [
                    motif for _ in regions_list
                ]
                ax = dm.plot_enrichment.plot_enrichment(
                    mod_file_names=[results["pileup"][0] for _ in regions_list],
                    regions_list=regions_list,
                    sample_names=["label" for _ in regions_list],
                    **kwargs_plot_enrichment_plot_enrichment,
                )
                assert isinstance(
                    ax, Axes
                ), f"{test_case}: plotting failed for {motif}."
        else:
            print(f"{test_case} skipped for plot_enrichment.plot_enrichment.")

    def test_unit__plot_enrichment_by_regions(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["pileup"][0] is not None:
            kwargs_plot_enrichment_by_regions = filter_kwargs_for_func(
                dm.plot_enrichment.by_regions, kwargs
            )
            for motif in kwargs["motifs"]:
                regions_list = (
                    kwargs["regions"]
                    if isinstance(kwargs["regions"], list)
                    else [kwargs["regions"]]
                )
                ax = dm.plot_enrichment.by_regions(
                    mod_file_name=results["pileup"][0],
                    regions_list=regions_list,
                    motif=motif,
                    sample_names=["label" for _ in regions_list],
                    **kwargs_plot_enrichment_by_regions,
                )
                assert isinstance(
                    ax, Axes
                ), f"{test_case}: plotting failed for {motif}."
        else:
            print(f"{test_case} skipped for plot_enrichment.by_regions.")

    def test_unit__plot_enrichment_by_modification(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["pileup"][0] is not None:
            kwargs_plot_enrichment_by_modification = filter_kwargs_for_func(
                dm.plot_enrichment.by_modification, kwargs
            )
            ax = dm.plot_enrichment.by_modification(
                mod_file_name=results["pileup"][0],
                **kwargs_plot_enrichment_by_modification,
            )
            assert isinstance(ax, Axes), f"{test_case}: plotting failed."
        else:
            print(f"{test_case} skipped for plot_enrichment.by_modification.")

    def test_unit__plot_enrichment_by_dataset(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["pileup"][0] is not None:
            kwargs_plot_enrichment_by_dataset = filter_kwargs_for_func(
                dm.plot_enrichment.by_dataset, kwargs
            )
            for motif in kwargs["motifs"]:
                ax = dm.plot_enrichment.by_dataset(
                    mod_file_names=[results["pileup"][0]],
                    motif=motif,
                    **kwargs_plot_enrichment_by_dataset,
                )
                assert isinstance(ax, Axes), f"{test_case}: plotting failed."
        else:
            print(f"{test_case} skipped for plot_enrichment.by_dataset.")


class TestPlotEnrichmentProfileSynthetic:
    """
    Tests plotting functionality in plot_enrichment_profile.

    This test simply checks that we can make plots from synthetic data without raising errors.
    Appearance of plots is not verified.
    """

    def test_unit__plot_enrichment_profile_plot_enrichment_profile_synthetic(self):
        ax = dm.plot_enrichment_profile.plot_enrichment_profile(
            mod_file_names=["test1.fake", "test2.fake"],
            regions_list=["test1.bed", "test2.bed"],
            motifs=["A", "C"],
            window_size=500,
            sample_names=["sample1", "sample2"],
            smooth_window=50,
        )
        assert isinstance(ax, Axes)

    def test_unit__plot_enrichment_profile_by_modification_synthetic(self):
        ax = dm.plot_enrichment_profile.by_modification(
            mod_file_name="test.fake",
            regions="test.bed",
            window_size=500,
            motifs=["A", "C"],
            smooth_window=50,
        )
        assert isinstance(ax, Axes)

    def test_unit__plot_enrichment_profile_by_region_synthetic(self):
        ax = dm.plot_enrichment_profile.by_regions(
            mod_file_name="test.fake",
            regions_list=["test1.bed", "test2.bed"],
            motif="A",
            window_size=500,
            sample_names=["on target", "off target"],
            smooth_window=50,
        )
        assert isinstance(ax, Axes)

    def test_unit__plot_enrichment_profile_by_dataset_synthetic(self):
        ax = dm.plot_enrichment_profile.by_dataset(
            mod_file_names=["test1.fake", "test2.fake"],
            regions="test.bed",
            motif="A",
            window_size=500,
            sample_names=["experiment 1", "experiment 2"],
            smooth_window=50,
        )
        assert isinstance(ax, Axes)


@pytest.mark.parametrize(
    "test_case,kwargs,results",
    [(case, inputs, outputs) for case, (inputs, outputs) in test_matrix.items()],
)
class TestPlotEnrichmentProfile:
    def test_unit__plot_enrichment_profile_plot_enrichment_profile(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["pileup"][0] is not None:
            kwargs_plot_enrichment_profile_plot_enrichment_profile = (
                filter_kwargs_for_func(
                    dm.plot_enrichment_profile.plot_enrichment_profile,
                    kwargs,
                    extra_args=["window_size", "smooth_window"],
                )
            )
            for motif in kwargs["motifs"]:
                regions_list = (
                    kwargs["regions"]
                    if isinstance(kwargs["regions"], list)
                    else [kwargs["regions"]]
                )
                kwargs_plot_enrichment_profile_plot_enrichment_profile["motifs"] = [
                    motif for _ in regions_list
                ]
                ax = dm.plot_enrichment_profile.plot_enrichment_profile(
                    mod_file_names=[results["pileup"][0] for _ in regions_list],
                    regions_list=regions_list,
                    sample_names=["label" for _ in regions_list],
                    **kwargs_plot_enrichment_profile_plot_enrichment_profile,
                )
                assert isinstance(
                    ax, Axes
                ), f"{test_case}: plotting failed for {motif}."
        else:
            print(
                f"{test_case} skipped for plot_enrichment_profile.plot_enrichment_profile."
            )

    def test_unit__plot_enrichment_profile_by_regions(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["pileup"][0] is not None:
            kwargs_plot_enrichment_profile_by_regions = filter_kwargs_for_func(
                dm.plot_enrichment_profile.by_regions,
                kwargs,
                extra_args=["window_size", "smooth_window"],
            )
            for motif in kwargs["motifs"]:
                regions_list = (
                    kwargs["regions"]
                    if isinstance(kwargs["regions"], list)
                    else [kwargs["regions"]]
                )
                ax = dm.plot_enrichment_profile.by_regions(
                    mod_file_name=results["pileup"][0],
                    regions_list=regions_list,
                    motif=motif,
                    sample_names=["label" for _ in regions_list],
                    **kwargs_plot_enrichment_profile_by_regions,
                )
                assert isinstance(
                    ax, Axes
                ), f"{test_case}: plotting failed for {motif}."
        else:
            print(f"{test_case} skipped for plot_enrichment_profile.by_regions.")

    def test_unit__plot_enrichment_profile_by_modification(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["pileup"][0] is not None:
            kwargs_plot_enrichment_profile_by_modification = filter_kwargs_for_func(
                dm.plot_enrichment_profile.by_modification,
                kwargs,
                extra_args=["window_size", "smooth_window"],
            )
            ax = dm.plot_enrichment_profile.by_modification(
                mod_file_name=results["pileup"][0],
                **kwargs_plot_enrichment_profile_by_modification,
            )
            assert isinstance(ax, Axes), f"{test_case}: plotting failed."
        else:
            print(f"{test_case} skipped for plot_enrichment_profile.by_modification.")

    def test_unit__plot_enrichment_by_dataset(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["pileup"][0] is not None:
            kwargs_plot_enrichment_profile_by_dataset = filter_kwargs_for_func(
                dm.plot_enrichment_profile.by_dataset,
                kwargs,
                extra_args=["window_size", "smooth_window"],
            )
            for motif in kwargs["motifs"]:
                ax = dm.plot_enrichment_profile.by_dataset(
                    mod_file_names=[results["pileup"][0]],
                    motif=motif,
                    **kwargs_plot_enrichment_profile_by_dataset,
                )
                assert isinstance(ax, Axes), f"{test_case}: plotting failed."
        else:
            print(f"{test_case} skipped for plot_enrichment_profile.by_dataset.")


class TestPlotDepthProfileSynthetic:
    """
    Tests plotting functionality in plot_depth_profile.

    This test simply checks that we can make plots from synthetic data without raising errors.
    Appearance of plots is not verified.
    """

    def test_unit__plot_depth_profile_plot_depth_profile_synthetic(self):
        ax = dm.plot_depth_profile.plot_depth_profile(
            mod_file_names=["test1.fake", "test2.fake"],
            regions_list=["test1.bed", "test2.bed"],
            motifs=["A", "C"],
            window_size=500,
            sample_names=["sample1", "sample2"],
            smooth_window=50,
        )
        assert isinstance(ax, Axes)

    def test_unit__plot_depth_profile_by_modification_synthetic(self):
        ax = dm.plot_depth_profile.by_modification(
            mod_file_name="test.fake",
            regions="test.bed",
            window_size=500,
            motifs=["A", "C"],
            smooth_window=50,
        )
        assert isinstance(ax, Axes)

    def test_unit__plot_depth_profile_by_region_synthetic(self):
        ax = dm.plot_depth_profile.by_regions(
            mod_file_name="test.fake",
            regions_list=["test1.bed", "test2.bed"],
            motif="A",
            window_size=500,
            sample_names=["on target", "off target"],
            smooth_window=50,
        )
        assert isinstance(ax, Axes)

    def test_unit__plot_depth_profile_by_dataset_synthetic(self):
        ax = dm.plot_depth_profile.by_dataset(
            mod_file_names=["test1.fake", "test2.fake"],
            regions="test.bed",
            motif="A",
            window_size=500,
            sample_names=["experiment 1", "experiment 2"],
            smooth_window=50,
        )
        assert isinstance(ax, Axes)


@pytest.mark.parametrize(
    "test_case,kwargs,results",
    [(case, inputs, outputs) for case, (inputs, outputs) in test_matrix.items()],
)
class TestPlotDepthProfile:
    def test_unit__plot_depth_profile_plot_depth_profile(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["pileup"][0] is not None:
            kwargs_plot_depth_profile_plot_depth_profile = filter_kwargs_for_func(
                dm.plot_depth_profile.plot_depth_profile,
                kwargs,
                extra_args=["window_size", "smooth_window"],
            )
            for motif in kwargs["motifs"]:
                regions_list = (
                    kwargs["regions"]
                    if isinstance(kwargs["regions"], list)
                    else [kwargs["regions"]]
                )
                kwargs_plot_depth_profile_plot_depth_profile["motifs"] = [
                    motif for _ in regions_list
                ]
                ax = dm.plot_depth_profile.plot_depth_profile(
                    mod_file_names=[results["pileup"][0] for _ in regions_list],
                    regions_list=regions_list,
                    sample_names=["label" for _ in regions_list],
                    **kwargs_plot_depth_profile_plot_depth_profile,
                )
                assert isinstance(
                    ax, Axes
                ), f"{test_case}: plotting failed for {motif}."
        else:
            print(f"{test_case} skipped for plot_depth_profile.plot_depth_profile.")

    def test_unit__plot_depth_profile_by_regions(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["pileup"][0] is not None:
            kwargs_plot_depth_profile_by_regions = filter_kwargs_for_func(
                dm.plot_depth_profile.by_regions,
                kwargs,
                extra_args=["window_size", "smooth_window"],
            )
            for motif in kwargs["motifs"]:
                regions_list = (
                    kwargs["regions"]
                    if isinstance(kwargs["regions"], list)
                    else [kwargs["regions"]]
                )
                ax = dm.plot_depth_profile.by_regions(
                    mod_file_name=results["pileup"][0],
                    regions_list=regions_list,
                    motif=motif,
                    sample_names=["label" for _ in regions_list],
                    **kwargs_plot_depth_profile_by_regions,
                )
                assert isinstance(
                    ax, Axes
                ), f"{test_case}: plotting failed for {motif}."
        else:
            print(f"{test_case} skipped for plot_depth_profile.by_regions.")

    def test_unit__plot_depth_profile_by_modification(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["pileup"][0] is not None:
            kwargs_plot_depth_profile_by_modification = filter_kwargs_for_func(
                dm.plot_depth_profile.by_modification,
                kwargs,
                extra_args=["window_size", "smooth_window"],
            )
            ax = dm.plot_depth_profile.by_modification(
                mod_file_name=results["pileup"][0],
                **kwargs_plot_depth_profile_by_modification,
            )
            assert isinstance(ax, Axes), f"{test_case}: plotting failed."
        else:
            print(f"{test_case} skipped for plot_depth_profile.by_modification.")

    def test_unit__plot_depth_by_dataset(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["pileup"][0] is not None:
            kwargs_plot_depth_profile_by_dataset = filter_kwargs_for_func(
                dm.plot_depth_profile.by_dataset,
                kwargs,
                extra_args=["window_size", "smooth_window"],
            )
            for motif in kwargs["motifs"]:
                ax = dm.plot_depth_profile.by_dataset(
                    mod_file_names=[results["pileup"][0]],
                    motif=motif,
                    **kwargs_plot_depth_profile_by_dataset,
                )
                assert isinstance(ax, Axes), f"{test_case}: plotting failed."
        else:
            print(f"{test_case} skipped for plot_depth_profile.by_dataset.")


class TestPlotDepthHistogramSynthetic:
    """
    Tests plotting functionality in plot_depth_histogram.

    This test simply checks that we can make plots from synthetic data without raising errors.
    Appearance of plots is not verified.
    """

    def test_unit__plot_depth_histogram_plot_depth_histogram_synthetic(self):
        ax = dm.plot_depth_histogram.plot_depth_histogram(
            mod_file_names=["test1.fake", "test2.fake"],
            regions_list=["test1.bed", "test2.bed"],
            motifs=["A", "C"],
            window_size=500,
            sample_names=["sample1", "sample2"],
        )
        assert isinstance(ax, Axes)

    def test_unit__plot_depth_histogram_by_modification_synthetic(self):
        ax = dm.plot_depth_histogram.by_modification(
            mod_file_name="test.fake",
            regions="test.bed",
            window_size=500,
            motifs=["A", "C"],
        )
        assert isinstance(ax, Axes)

    def test_unit__plot_depth_histogram_by_region_synthetic(self):
        ax = dm.plot_depth_histogram.by_regions(
            mod_file_name="test.fake",
            regions_list=["test1.bed", "test2.bed"],
            motif="A",
            window_size=500,
            sample_names=["on target", "off target"],
        )
        assert isinstance(ax, Axes)

    def test_unit__plot_depth_histogram_by_dataset_synthetic(self):
        ax = dm.plot_depth_histogram.by_dataset(
            mod_file_names=["test1.fake", "test2.fake"],
            regions="test.bed",
            motif="A",
            window_size=500,
            sample_names=["experiment 1", "experiment 2"],
        )
        assert isinstance(ax, Axes)


@pytest.mark.parametrize(
    "test_case,kwargs,results",
    [(case, inputs, outputs) for case, (inputs, outputs) in test_matrix.items()],
)
class TestPlotDepthHistogram:
    def test_unit__plot_depth_histogram_plot_depth_histogram(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["pileup"][0] is not None:
            kwargs_plot_depth_histogram_plot_depth_histogram = filter_kwargs_for_func(
                dm.plot_depth_histogram.plot_depth_histogram,
                kwargs,
                extra_args=["window_size"],
            )
            for motif in kwargs["motifs"]:
                regions_list = (
                    kwargs["regions"]
                    if isinstance(kwargs["regions"], list)
                    else [kwargs["regions"]]
                )
                kwargs_plot_depth_histogram_plot_depth_histogram["motifs"] = [
                    motif for _ in regions_list
                ]
                ax = dm.plot_depth_histogram.plot_depth_histogram(
                    mod_file_names=[results["pileup"][0] for _ in regions_list],
                    regions_list=regions_list,
                    sample_names=["label" for _ in regions_list],
                    **kwargs_plot_depth_histogram_plot_depth_histogram,
                )
                assert isinstance(
                    ax, Axes
                ), f"{test_case}: plotting failed for {motif}."
        else:
            print(f"{test_case} skipped for plot_depth_histogram.plot_depth_histogram.")

    def test_unit__plot_depth_histogram_by_regions(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["pileup"][0] is not None:
            kwargs_plot_depth_histogram_by_regions = filter_kwargs_for_func(
                dm.plot_depth_histogram.by_regions,
                kwargs,
                extra_args=["window_size"],
            )
            for motif in kwargs["motifs"]:
                regions_list = (
                    kwargs["regions"]
                    if isinstance(kwargs["regions"], list)
                    else [kwargs["regions"]]
                )
                ax = dm.plot_depth_histogram.by_regions(
                    mod_file_name=results["pileup"][0],
                    regions_list=regions_list,
                    motif=motif,
                    sample_names=["label" for _ in regions_list],
                    **kwargs_plot_depth_histogram_by_regions,
                )
                assert isinstance(
                    ax, Axes
                ), f"{test_case}: plotting failed for {motif}."
        else:
            print(f"{test_case} skipped for plot_depth_histogram.by_regions.")

    def test_unit__plot_depth_histogram_by_modification(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["pileup"][0] is not None:
            kwargs_plot_depth_histogram_by_modification = filter_kwargs_for_func(
                dm.plot_depth_histogram.by_modification,
                kwargs,
                extra_args=["window_size"],
            )
            ax = dm.plot_depth_histogram.by_modification(
                mod_file_name=results["pileup"][0],
                **kwargs_plot_depth_histogram_by_modification,
            )
            assert isinstance(ax, Axes), f"{test_case}: plotting failed."
        else:
            print(f"{test_case} skipped for plot_depth_histogram.by_modification.")

    def test_unit__plot_depth_by_dataset(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["pileup"][0] is not None:
            kwargs_plot_depth_histogram_by_dataset = filter_kwargs_for_func(
                dm.plot_depth_histogram.by_dataset,
                kwargs,
                extra_args=["window_size"],
            )
            for motif in kwargs["motifs"]:
                ax = dm.plot_depth_histogram.by_dataset(
                    mod_file_names=[results["pileup"][0]],
                    motif=motif,
                    **kwargs_plot_depth_histogram_by_dataset,
                )
                assert isinstance(ax, Axes), f"{test_case}: plotting failed."
        else:
            print(f"{test_case} skipped for plot_depth_histogram.by_dataset.")


class TestPlotReadsSynthetic:
    """
    Tests plotting functionality in plot_reads.

    This test simply checks that we can make plots from synthetic data without raising errors.
    Appearance of plots is not verified.
    """

    def test_unit__plot_reads_plot_reads_synthetic(self):
        ax = dm.plot_reads.plot_reads(
            mod_file_name="test.fake", regions="test.bed", motifs=["A,0", "CG,0"]
        )
        assert isinstance(ax, Axes)


@pytest.mark.parametrize(
    "test_case,kwargs,results",
    [(case, inputs, outputs) for case, (inputs, outputs) in test_matrix.items()],
)
class TestPlotReads:
    def test_unit__plot_reads_plot_reads(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["extract"][0] is not None:
            kwargs_plot_reads_plot_reads = filter_kwargs_for_func(
                dm.plot_reads.plot_reads, kwargs
            )
            if kwargs["thresh"] is not None:
                ax = dm.plot_reads.plot_reads(
                    mod_file_name=results["extract"][0],
                    **kwargs_plot_reads_plot_reads,
                )
                assert isinstance(ax, Axes), f"{test_case}: plotting failed."
            else:  # if the extract parameters did not have a threshold, plot_reads.plot_reads should raise an error
                with pytest.raises(ValueError) as excinfo:
                    ax = dm.plot_reads.plot_reads(
                        mod_file_name=results["extract"][0],
                        **kwargs_plot_reads_plot_reads,
                    )
                assert "No threshold has been applied" in str(
                    excinfo.value
                ), f"{test_case}: unexpected exception {excinfo.value}"
                # providing a threshold should be enough to run plot_reads.plot_reads without an error
                kwargs_plot_reads_plot_reads["thresh"] = 0.75
                ax = dm.plot_reads.plot_reads(
                    mod_file_name=results["extract"][0],
                    **kwargs_plot_reads_plot_reads,
                )
                assert isinstance(ax, Axes), f"{test_case}: plotting failed."
        else:
            print(f"{test_case} skipped for test_unit__plot_reads_plot_reads.")


@pytest.mark.parametrize(
    "test_case,kwargs,results",
    [(case, inputs, outputs) for case, (inputs, outputs) in test_matrix.items()],
)
class TestPlotReadBrowser:
    def test_unit__plot_read_browser(
        self,
        test_case,
        kwargs,
        results,
    ):
        if results["extract"][0] is not None:
            kwargs_plot_read_browser = filter_kwargs_for_func(
                dm.plot_read_browser.plot_read_browser, kwargs
            )
            if (
                kwargs["thresh"] is None
                and kwargs["thresh"] is None
                and not isinstance(kwargs["regions"], list)
                and Path(kwargs["regions"]).suffix != ".bed"
            ):
                fig = dm.plot_read_browser.plot_read_browser(
                    mod_file_name=results["extract"][0],
                    region=kwargs["regions"],
                    **kwargs_plot_read_browser,
                )
                assert isinstance(
                    fig, plotly.graph_objs.Figure
                ), f"{test_case}: plotting failed."
            else:
                with pytest.raises(ValueError) as excinfo:
                    fig = dm.plot_read_browser.plot_read_browser(
                        mod_file_name=results["extract"][0],
                        region=kwargs["regions"],
                        **kwargs_plot_read_browser,
                    )
                if (
                    isinstance(kwargs["regions"], list)
                    or Path(kwargs["regions"]).suffix == ".bed"
                ) and kwargs["thresh"] is None:
                    assert (
                        "Invalid region" in str(excinfo.value)
                    ), f"{test_case}: unexpected exception for no-threshold bad-region case {excinfo.value}"
                elif (
                    kwargs["thresh"] is not None
                    and not isinstance(kwargs["regions"], list)
                    and Path(kwargs["regions"]).suffix != ".bed"
                ):
                    assert (
                        "A threshold has been applied" in str(excinfo.value)
                    ), f"{test_case}: unexpected exception thresholded valid-region case {excinfo.value}"
                else:
                    assert (
                        "A threshold has been applied" in str(excinfo.value)
                        or "Invalid region" in str(excinfo.value)
                    ), f"{test_case}: unexpected exception thresholded bad-region case {excinfo.value}"

        else:
            print(f"{test_case} skipped for test_unit__plot_read_browser")
