from __future__ import annotations

import gzip

import pandas as pd
import pytest

from dimelo import encode


def test_parse_encode_search():
    payload = {
        "@graph": [
            {
                "accession": "ENCSR000AAA",
                "target": {"label": "CTCF"},
                "assay_title": "TF ChIP-seq",
                "biosample_ontology": {"term_name": "GM12878"},
                "assembly": ["GRCh38"],
                "description": "CTCF ChIP",
            },
            {
                "accession": "ENCSR000BBB",
                "target": {"label": "H3K27ac"},
                "assay_title": "Histone ChIP-seq",
                "biosample_ontology": {"term_name": "K562"},
                "assembly": "GRCh38",
            },
        ]
    }
    table = encode.parse_encode_search(payload).set_index("accession")
    assert table.loc["ENCSR000AAA", "target"] == "CTCF"
    assert table.loc["ENCSR000AAA", "assembly"] == "GRCh38"  # list unwrapped
    assert table.loc["ENCSR000BBB", "biosample"] == "K562"
    assert pd.isna(table.loc["ENCSR000BBB", "description"])


def test_parse_encode_files_filters_and_builds_url():
    payload = {
        "files": [
            {
                "accession": "ENCFF001",
                "file_format": "bigWig",
                "output_type": "fold change over control",
                "assembly": "GRCh38",
                "file_size": 100,
                "href": "/files/ENCFF001/@@download/ENCFF001.bigWig",
            },
            {
                "accession": "ENCFF002",
                "file_format": "bed",
                "output_type": "IDR thresholded peaks",
                "assembly": "GRCh38",
                "href": "/files/ENCFF002/@@download/ENCFF002.bed.gz",
            },
            {
                "accession": "ENCFF003",
                "file_format": "bigWig",
                "output_type": "signal",
                "assembly": "hg19",
                "href": "/files/ENCFF003/@@download/ENCFF003.bigWig",
            },
        ]
    }
    # unfiltered: absolute urls built from href
    all_files = encode.parse_encode_files(payload)
    assert (
        all_files.set_index("accession").loc["ENCFF001", "url"]
        == "https://www.encodeproject.org/files/ENCFF001/@@download/ENCFF001.bigWig"
    )
    # filter by format + assembly
    bigwig_grch38 = encode.parse_encode_files(
        payload, file_format="bigWig", assembly="GRCh38"
    )
    assert list(bigwig_grch38["accession"]) == ["ENCFF001"]
    # filter by output_type
    peaks = encode.parse_encode_files(payload, output_type="IDR thresholded peaks")
    assert list(peaks["accession"]) == ["ENCFF002"]


def test_load_encode_peaks_plain_and_gzip(tmp_path):
    narrowpeak = (
        "chr1\t100\t200\tpeak1\t0\t.\t5.0\t3.0\t2.0\t50\n"
        "chr2\t300\t400\tpeak2\t0\t.\t4.0\t2.5\t1.5\t60\n"
    )
    plain = tmp_path / "peaks.bed"
    plain.write_text(narrowpeak)
    table = encode.load_encode_peaks(plain)
    assert list(table.columns) == ["chromosome", "start", "end"]
    assert list(table["chromosome"]) == ["chr1", "chr2"]
    assert list(table["start"]) == [100, 300]

    gz = tmp_path / "peaks.bed.gz"
    with gzip.open(gz, "wt") as handle:
        handle.write(narrowpeak)
    gz_table = encode.load_encode_peaks(gz, max_peaks=1)
    assert len(gz_table) == 1
    assert gz_table.loc[0, "end"] == 200


def test_download_encode_file_cache_hit_no_network(tmp_path):
    # pre-populate the cache so no network request is attempted
    cached = tmp_path / "ENCFF001.bigWig"
    cached.write_bytes(b"stub")
    url = "https://www.encodeproject.org/files/ENCFF001/@@download/ENCFF001.bigWig"
    result = encode.download_encode_file(url, cache_dir=tmp_path)
    assert result == cached
    assert result.read_bytes() == b"stub"


def test_encode_peaks_feed_cluster_reference(tmp_path):
    # end-to-end wiring: ENCODE peaks -> cluster_reference association
    from dimelo import cluster_reference

    (tmp_path / "peaks.bed").write_text("chr1\t0\t100\nchr1\t300\t400\n")
    peaks = encode.load_encode_peaks(tmp_path / "peaks.bed")
    sites = pd.DataFrame(
        {
            "cluster": ["C0", "C0", "C1", "C1"],
            "chromosome": ["chr1", "chr1", "chr1", "chr1"],
            "start": [0, 300, 5000, 6000],
            "end": [100, 400, 5100, 6100],
        }
    )
    enrichment = cluster_reference.annotate_clusters_by_site_groups(
        sites, {"encode_peaks": peaks}
    ).set_index("cluster")
    assert enrichment.loc["C0", "feature_fraction"] == pytest.approx(1.0)
    assert enrichment.loc["C1", "feature_fraction"] == pytest.approx(0.0)
