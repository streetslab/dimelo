import gzip
from pathlib import Path

import h5py
import numpy as np

from dimelo import parse_bam


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


def _decode_vector(dataset, index: int) -> np.ndarray:
    raw_bytes = np.asarray(dataset[index], dtype=np.uint8).tobytes()
    return np.frombuffer(gzip.decompress(raw_bytes), dtype=np.uint8)


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
