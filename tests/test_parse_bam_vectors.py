import gzip
from pathlib import Path

import h5py
import numpy as np
import pytest

from dimelo import load_processed, parse_bam


def _make_extract_line(
    read_name: str,
    pos_in_read: int,
    pos_in_genome: int,
    chromosome: str,
    strand: str,
    read_len: int,
    prob: float,
    mod_code: str,
    canonical_base: str,
) -> str:
    fields = [
        read_name,
        str(pos_in_read),
        str(pos_in_genome),
        chromosome,
        ".",
        strand,
        ".",
        ".",
        ".",
        str(read_len),
        str(prob),
        mod_code,
        ".",
        ".",
        ".",
        canonical_base,
    ]
    return "\t".join(fields)


def _write_extract_file(path: Path) -> None:
    header = "\t".join(f"col{i}" for i in range(16))
    lines = [
        header,
        _make_extract_line(
            "read_sparse",
            pos_in_read=0,
            pos_in_genome=100,
            chromosome="chr1",
            strand="+",
            read_len=10,
            prob=0.10,
            mod_code="a",
            canonical_base="A",
        ),
        _make_extract_line(
            "read_sparse",
            pos_in_read=8,
            pos_in_genome=108,
            chromosome="chr1",
            strand="+",
            read_len=10,
            prob=0.90,
            mod_code="a",
            canonical_base="A",
        ),
        _make_extract_line(
            "read_no_hits",
            pos_in_read=0,
            pos_in_genome=200,
            chromosome="chr1",
            strand="+",
            read_len=10,
            prob=0.20,
            mod_code="m",
            canonical_base="C",
        ),
        _make_extract_line(
            "read_no_hits",
            pos_in_read=9,
            pos_in_genome=210,
            chromosome="chr1",
            strand="+",
            read_len=10,
            prob=0.30,
            mod_code="m",
            canonical_base="C",
        ),
    ]
    path.write_text("\n".join(lines) + "\n")


def _write_coordinate_stress_extract_file(path: Path) -> None:
    header = "\t".join(f"col{i}" for i in range(16))
    lines = [
        header,
        # plus strand read with inferred start drift across rows (tests coordinate robustness)
        _make_extract_line(
            "read_plus_drift",
            pos_in_read=1,
            pos_in_genome=102,
            chromosome="chr2",
            strand="+",
            read_len=10,
            prob=0.90,
            mod_code="a",
            canonical_base="A",
        ),
        _make_extract_line(
            "read_plus_drift",
            pos_in_read=0,
            pos_in_genome=100,
            chromosome="chr2",
            strand="+",
            read_len=10,
            prob=0.80,
            mod_code="a",
            canonical_base="A",
        ),
        # minus strand read where reference-oriented positions should still map left-to-right
        _make_extract_line(
            "read_minus",
            pos_in_read=1,
            pos_in_genome=208,
            chromosome="chr2",
            strand="-",
            read_len=10,
            prob=0.70,
            mod_code="a",
            canonical_base="A",
        ),
        _make_extract_line(
            "read_minus",
            pos_in_read=9,
            pos_in_genome=200,
            chromosome="chr2",
            strand="-",
            read_len=10,
            prob=0.95,
            mod_code="a",
            canonical_base="A",
        ),
    ]
    path.write_text("\n".join(lines) + "\n")


def _decode_vector(dataset, index: int) -> np.ndarray:
    raw_bytes = np.asarray(dataset[index], dtype=np.uint8).tobytes()
    return np.frombuffer(gzip.decompress(raw_bytes), dtype=np.uint8)


def _assert_read_tuples_equal(left: list[tuple], right: list[tuple]) -> None:
    assert len(left) == len(right)
    for left_tuple, right_tuple in zip(left, right, strict=True):
        assert len(left_tuple) == len(right_tuple)
        for left_value, right_value in zip(left_tuple, right_tuple, strict=True):
            if isinstance(left_value, np.ndarray):
                np.testing.assert_array_equal(left_value, right_value)
            else:
                assert left_value == right_value


def test_read_by_base_txt_to_hdf5_uses_dense_span_vectors(tmp_path):
    input_txt = tmp_path / "extract.txt"
    output_h5 = tmp_path / "reads.h5"
    _write_extract_file(input_txt)

    parse_bam.read_by_base_txt_to_hdf5(
        input_txt=input_txt,
        output_h5=output_h5,
        motif="A,0",
        thresh=None,
        quiet=True,
    )

    with h5py.File(output_h5, "r") as h5:
        sparse_index = 0
        sparse_span = int(h5["read_end"][sparse_index] - h5["read_start"][sparse_index])
        sparse_mod_vector = _decode_vector(h5["mod_vector"], sparse_index)
        sparse_val_vector = _decode_vector(h5["val_vector"], sparse_index)

        no_hits_index = 1
        no_hits_span = int(h5["read_end"][no_hits_index] - h5["read_start"][no_hits_index])
        no_hits_mod_vector = _decode_vector(h5["mod_vector"], no_hits_index)
        no_hits_val_vector = _decode_vector(h5["val_vector"], no_hits_index)

    assert sparse_span == 10
    assert len(sparse_mod_vector) == sparse_span
    assert len(sparse_val_vector) == sparse_span

    assert no_hits_span == 11
    assert len(no_hits_mod_vector) == no_hits_span
    assert len(no_hits_val_vector) == no_hits_span
    assert np.all(no_hits_mod_vector == 0)
    assert np.all(no_hits_val_vector == 0)


def test_read_by_base_txt_to_hdf5_thresholds_at_write_time(tmp_path):
    input_txt = tmp_path / "extract.txt"
    output_h5 = tmp_path / "reads_thresholded.h5"
    _write_extract_file(input_txt)

    parse_bam.read_by_base_txt_to_hdf5(
        input_txt=input_txt,
        output_h5=output_h5,
        motif="A,0",
        thresh=0.5,
        quiet=True,
    )

    with h5py.File(output_h5, "r") as h5:
        assert np.isclose(h5["threshold"][()], 0.5)

        sparse_index = 0
        sparse_mod_vector = _decode_vector(h5["mod_vector"], sparse_index)
        sparse_val_vector = _decode_vector(h5["val_vector"], sparse_index)

    assert set(np.unique(sparse_mod_vector)).issubset({0, 1})
    assert set(np.unique(sparse_val_vector)).issubset({0, 1})
    np.testing.assert_array_equal(
        sparse_mod_vector,
        np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 0], dtype=np.uint8),
    )
    np.testing.assert_array_equal(
        sparse_val_vector,
        np.array([1, 0, 0, 0, 0, 0, 0, 0, 1, 0], dtype=np.uint8),
    )


def test_read_vectors_from_hdf5_loads_thresholded_vectors_as_binary(tmp_path):
    input_txt = tmp_path / "extract.txt"
    output_h5 = tmp_path / "reads_thresholded.h5"
    _write_extract_file(input_txt)

    parse_bam.read_by_base_txt_to_hdf5(
        input_txt=input_txt,
        output_h5=output_h5,
        motif="A,0",
        thresh=0.5,
        quiet=True,
    )

    read_tuples, datasets, _ = load_processed.read_vectors_from_hdf5(
        file=output_h5,
        motifs=["A,0"],
        calculate_mod_fractions=False,
        quiet=True,
    )

    mod_vector_index = datasets.index("mod_vector")
    val_vector_index = datasets.index("val_vector")
    read_name_index = datasets.index("read_name")

    assert len(read_tuples) == 2
    assert all(read_tuple[mod_vector_index].dtype == np.bool_ for read_tuple in read_tuples)
    assert all(read_tuple[val_vector_index].dtype == np.bool_ for read_tuple in read_tuples)

    reads_by_name = {
        read_tuple[read_name_index]: read_tuple for read_tuple in read_tuples
    }

    sparse_read = reads_by_name["read_sparse"]
    no_hits_read = reads_by_name["read_no_hits"]

    np.testing.assert_array_equal(
        sparse_read[mod_vector_index],
        np.array([False, False, False, False, False, False, False, False, True, False]),
    )
    np.testing.assert_array_equal(
        sparse_read[val_vector_index],
        np.array([True, False, False, False, False, False, False, False, True, False]),
    )
    assert not np.any(no_hits_read[mod_vector_index])
    assert not np.any(no_hits_read[val_vector_index])

    mod_positions, _, _, _ = load_processed.readwise_binary_modification_arrays(
        file=output_h5,
        motifs=["A,0"],
        regions="chr1:95-220",
        thresh=None,
        quiet=True,
    )
    assert len(mod_positions) == 1


def test_read_by_base_txt_to_hdf5_reconstructs_coordinates_robustly(tmp_path):
    input_txt = tmp_path / "extract.coordinate_stress.txt"
    output_h5 = tmp_path / "reads.coordinate_stress.h5"
    _write_coordinate_stress_extract_file(input_txt)

    parse_bam.read_by_base_txt_to_hdf5(
        input_txt=input_txt,
        output_h5=output_h5,
        motif="A,0",
        thresh=0.5,
        quiet=True,
    )

    with h5py.File(output_h5, "r") as h5:
        read_names = [
            name.decode() if isinstance(name, bytes) else str(name)
            for name in h5["read_name"][:]
        ]
        read_name_to_index = {name: idx for idx, name in enumerate(read_names)}

        plus_idx = read_name_to_index["read_plus_drift"]
        minus_idx = read_name_to_index["read_minus"]

        assert int(h5["read_start"][plus_idx]) == 100
        assert int(h5["read_end"][plus_idx]) == 111
        plus_mod = _decode_vector(h5["mod_vector"], plus_idx)
        plus_val = _decode_vector(h5["val_vector"], plus_idx)
        assert len(plus_mod) == 11
        np.testing.assert_array_equal(
            plus_mod, np.array([1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], dtype=np.uint8)
        )
        np.testing.assert_array_equal(
            plus_val, np.array([1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], dtype=np.uint8)
        )

        assert int(h5["read_start"][minus_idx]) == 200
        assert int(h5["read_end"][minus_idx]) == 210
        minus_mod = _decode_vector(h5["mod_vector"], minus_idx)
        minus_val = _decode_vector(h5["val_vector"], minus_idx)
        assert len(minus_mod) == 10
        np.testing.assert_array_equal(
            minus_mod, np.array([1, 0, 0, 0, 0, 0, 0, 0, 1, 0], dtype=np.uint8)
        )
        np.testing.assert_array_equal(
            minus_val, np.array([1, 0, 0, 0, 0, 0, 0, 0, 1, 0], dtype=np.uint8)
        )


def test_read_by_base_txt_to_hdf5_rejects_out_of_bounds_pos_in_read(tmp_path):
    input_txt = tmp_path / "extract.bad_position.txt"
    output_h5 = tmp_path / "reads.bad_position.h5"
    header = "\t".join(f"col{i}" for i in range(16))
    lines = [
        header,
        _make_extract_line(
            "read_bad",
            pos_in_read=10,
            pos_in_genome=110,
            chromosome="chr2",
            strand="+",
            read_len=10,
            prob=0.80,
            mod_code="a",
            canonical_base="A",
        ),
    ]
    input_txt.write_text("\n".join(lines) + "\n")

    with pytest.raises(ValueError, match="out of bounds for read length"):
        parse_bam.read_by_base_txt_to_hdf5(
            input_txt=input_txt,
            output_h5=output_h5,
            motif="A,0",
            thresh=0.5,
            quiet=True,
        )


def test_read_by_base_txt_to_hdf5_rejects_invalid_strand(tmp_path):
    input_txt = tmp_path / "extract.bad_strand.txt"
    output_h5 = tmp_path / "reads.bad_strand.h5"
    header = "\t".join(f"col{i}" for i in range(16))
    lines = [
        header,
        _make_extract_line(
            "read_bad_strand",
            pos_in_read=2,
            pos_in_genome=102,
            chromosome="chr2",
            strand="?",
            read_len=10,
            prob=0.80,
            mod_code="a",
            canonical_base="A",
        ),
    ]
    input_txt.write_text("\n".join(lines) + "\n")

    with pytest.raises(ValueError, match="Unexpected strand"):
        parse_bam.read_by_base_txt_to_hdf5(
            input_txt=input_txt,
            output_h5=output_h5,
            motif="A,0",
            thresh=0.5,
            quiet=True,
        )


def test_read_by_base_txt_to_hdf5_rejects_non_positive_read_length(tmp_path):
    input_txt = tmp_path / "extract.bad_read_length.txt"
    output_h5 = tmp_path / "reads.bad_read_length.h5"
    header = "\t".join(f"col{i}" for i in range(16))
    lines = [
        header,
        _make_extract_line(
            "read_bad_length",
            pos_in_read=0,
            pos_in_genome=100,
            chromosome="chr2",
            strand="+",
            read_len=0,
            prob=0.80,
            mod_code="a",
            canonical_base="A",
        ),
    ]
    input_txt.write_text("\n".join(lines) + "\n")

    with pytest.raises(ValueError, match="read length must be positive"):
        parse_bam.read_by_base_txt_to_hdf5(
            input_txt=input_txt,
            output_h5=output_h5,
            motif="A,0",
            thresh=0.5,
            quiet=True,
        )


def test_readwise_binary_modification_arrays_returns_empty_for_empty_region(tmp_path):
    input_txt = tmp_path / "extract.txt"
    output_h5 = tmp_path / "reads_thresholded.h5"
    _write_extract_file(input_txt)

    parse_bam.read_by_base_txt_to_hdf5(
        input_txt=input_txt,
        output_h5=output_h5,
        motif="A,0",
        thresh=0.5,
        quiet=True,
    )

    mod_positions, read_ids, motifs, regions_dict = (
        load_processed.readwise_binary_modification_arrays(
            file=output_h5,
            motifs=["A,0"],
            regions="chr1:1000-1010",
            thresh=None,
            quiet=True,
        )
    )

    assert mod_positions.size == 0
    assert read_ids.size == 0
    assert motifs.size == 0
    assert regions_dict == {"chr1": [(1000, 1010, ".")]}


def test_read_vectors_from_hdf5_subset_empty_region_returns_empty(tmp_path):
    input_txt = tmp_path / "extract.txt"
    output_h5 = tmp_path / "reads_thresholded.h5"
    _write_extract_file(input_txt)

    parse_bam.read_by_base_txt_to_hdf5(
        input_txt=input_txt,
        output_h5=output_h5,
        motif="A,0",
        thresh=0.5,
        quiet=True,
    )

    read_tuples, datasets, regions_dict = load_processed.read_vectors_from_hdf5(
        file=output_h5,
        motifs=["A,0"],
        regions="chr1:1000-1010",
        subset_parameters={"n": 1},
        calculate_mod_fractions=False,
        quiet=True,
    )

    assert read_tuples == []
    assert "mod_vector" in datasets
    assert regions_dict == {"chr1": [(1000, 1010, ".")]}


def test_readwise_binary_modification_arrays_splits_duplicate_read_names_by_region(
    tmp_path,
):
    input_txt = tmp_path / "extract.txt"
    output_h5 = tmp_path / "reads_thresholded.h5"
    _write_extract_file(input_txt)

    parse_bam.read_by_base_txt_to_hdf5(
        input_txt=input_txt,
        output_h5=output_h5,
        motif="A,0",
        thresh=0.5,
        quiet=True,
    )

    _, read_ids, motifs, regions_dict = load_processed.readwise_binary_modification_arrays(
        file=output_h5,
        motifs=["A,0"],
        regions=["chr1:95-105", "chr1:100-110"],
        thresh=None,
        quiet=True,
    )

    assert set(read_ids.tolist()) == {0, 1}
    assert motifs.size == read_ids.size == 2
    assert regions_dict == {"chr1": [(95, 105, "."), (100, 110, ".")]}


def test_read_vectors_from_hdf5_rejects_subset_without_n_or_frac(tmp_path):
    input_txt = tmp_path / "extract.txt"
    output_h5 = tmp_path / "reads_thresholded.h5"
    _write_extract_file(input_txt)

    parse_bam.read_by_base_txt_to_hdf5(
        input_txt=input_txt,
        output_h5=output_h5,
        motif="A,0",
        thresh=0.5,
        quiet=True,
    )

    with pytest.raises(ValueError, match="subset_parameters must include"):
        load_processed.read_vectors_from_hdf5(
            file=output_h5,
            motifs=["A,0"],
            regions="chr1:95-110",
            subset_parameters={},
            calculate_mod_fractions=False,
            quiet=True,
        )


def test_read_vectors_from_hdf5_rejects_non_dict_subset_parameters(tmp_path):
    input_txt = tmp_path / "extract.txt"
    output_h5 = tmp_path / "reads_thresholded.h5"
    _write_extract_file(input_txt)

    parse_bam.read_by_base_txt_to_hdf5(
        input_txt=input_txt,
        output_h5=output_h5,
        motif="A,0",
        thresh=0.5,
        quiet=True,
    )

    with pytest.raises(ValueError, match="provided as a dictionary"):
        load_processed.read_vectors_from_hdf5(
            file=output_h5,
            motifs=["A,0"],
            regions="chr1:95-110",
            subset_parameters=1,
            calculate_mod_fractions=False,
            quiet=True,
        )


def test_read_vectors_from_hdf5_rejects_subset_with_array_key(tmp_path):
    input_txt = tmp_path / "extract.txt"
    output_h5 = tmp_path / "reads_thresholded.h5"
    _write_extract_file(input_txt)

    parse_bam.read_by_base_txt_to_hdf5(
        input_txt=input_txt,
        output_h5=output_h5,
        motif="A,0",
        thresh=0.5,
        quiet=True,
    )

    with pytest.raises(ValueError, match="cannot include 'array'"):
        load_processed.read_vectors_from_hdf5(
            file=output_h5,
            motifs=["A,0"],
            regions="chr1:95-110",
            subset_parameters={"array": np.array([0]), "n": 1},
            calculate_mod_fractions=False,
            quiet=True,
        )


def test_read_vectors_from_hdf5_accepts_tuple_sort_by(tmp_path):
    input_txt = tmp_path / "extract.txt"
    output_h5 = tmp_path / "reads_thresholded.h5"
    _write_extract_file(input_txt)

    parse_bam.read_by_base_txt_to_hdf5(
        input_txt=input_txt,
        output_h5=output_h5,
        motif="A,0",
        thresh=0.5,
        quiet=True,
    )

    read_tuples, datasets, _ = load_processed.read_vectors_from_hdf5(
        file=output_h5,
        motifs=["A,0"],
        sort_by=("read_start", "desc"),
        calculate_mod_fractions=False,
        quiet=True,
    )

    read_name_index = datasets.index("read_name")
    assert [read_tuple[read_name_index] for read_tuple in read_tuples] == [
        "read_no_hits",
        "read_sparse",
    ]


def test_readwise_binary_modification_arrays_passes_subset_and_quiet(monkeypatch):
    captured_kwargs = {}

    def _fake_read_vectors_from_hdf5(**kwargs):
        captured_kwargs.update(kwargs)
        return (
            [
                (
                    "read1",
                    np.array([True, False], dtype=np.bool_),
                    "A,0",
                    100,
                    110,
                    100,
                    ".",
                )
            ],
            [
                "read_name",
                "mod_vector",
                "motif",
                "region_start",
                "region_end",
                "read_start",
                "region_strand",
            ],
            {"chr1": [(100, 110, ".")]},
        )

    monkeypatch.setattr(
        load_processed, "read_vectors_from_hdf5", _fake_read_vectors_from_hdf5
    )

    mod_positions, read_ids, motifs, regions_dict = (
        load_processed.readwise_binary_modification_arrays(
            file=Path("dummy.h5"),
            motifs=["A,0"],
            regions="chr1:100-110",
            subset_parameters={"n": 1},
            quiet=False,
        )
    )

    assert captured_kwargs["quiet"] is False
    assert captured_kwargs["subset_parameters"] == {"n": 1}
    assert mod_positions.tolist() == [-5]
    assert read_ids.tolist() == [0]
    assert motifs.tolist() == ["A,0"]
    assert regions_dict == {"chr1": [(100, 110, ".")]}


def test_readwise_binary_modification_arrays_parallel_matches_serial(tmp_path):
    input_txt = tmp_path / "extract.txt"
    output_h5 = tmp_path / "reads_thresholded.h5"
    _write_extract_file(input_txt)

    parse_bam.read_by_base_txt_to_hdf5(
        input_txt=input_txt,
        output_h5=output_h5,
        motif="A,0",
        thresh=0.5,
        quiet=True,
    )

    serial = load_processed.readwise_binary_modification_arrays(
        file=output_h5,
        motifs=["A,0"],
        regions="chr1:95-220",
        thresh=None,
        quiet=True,
        cores=1,
    )
    parallel = load_processed.readwise_binary_modification_arrays(
        file=output_h5,
        motifs=["A,0"],
        regions="chr1:95-220",
        thresh=None,
        quiet=True,
        cores=2,
    )

    np.testing.assert_array_equal(serial[0], parallel[0])
    np.testing.assert_array_equal(serial[1], parallel[1])
    np.testing.assert_array_equal(serial[2], parallel[2])
    assert serial[3] == parallel[3]


def test_read_vectors_from_hdf5_parallel_region_selection_matches_serial(tmp_path):
    input_txt = tmp_path / "extract.txt"
    output_h5 = tmp_path / "reads_thresholded.h5"
    _write_extract_file(input_txt)

    parse_bam.read_by_base_txt_to_hdf5(
        input_txt=input_txt,
        output_h5=output_h5,
        motif="A,0",
        thresh=0.5,
        quiet=True,
    )

    serial = load_processed.read_vectors_from_hdf5(
        file=output_h5,
        motifs=["A,0"],
        regions=["chr1:95-105", "chr1:100-110"],
        quiet=True,
        calculate_mod_fractions=False,
        cores=1,
    )
    parallel = load_processed.read_vectors_from_hdf5(
        file=output_h5,
        motifs=["A,0"],
        regions=["chr1:95-105", "chr1:100-110"],
        quiet=True,
        calculate_mod_fractions=False,
        cores=2,
    )

    _assert_read_tuples_equal(serial[0], parallel[0])
    assert serial[1] == parallel[1]
    assert serial[2] == parallel[2]
