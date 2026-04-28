import pandas as pd
import pytest

from dimelo import chip_atlas


def test_region_ids_to_bed_dataframe_parses_region_ids():
    frame = chip_atlas.region_ids_to_bed_dataframe(["chr1:100-200,+", "chr2:500-700"])
    assert frame.shape[0] == 2
    assert frame.loc[0, "chrom"] == "chr1"
    assert int(frame.loc[0, "start"]) == 100
    assert int(frame.loc[0, "end"]) == 200
    assert frame.loc[0, "strand"] == "+"
    assert frame.loc[1, "strand"] == "."


def test_submit_enrichment_extracts_request_id_from_json(monkeypatch):
    captured = {}

    def fake_request(**kwargs):
        captured.update(kwargs)
        return 200, "application/json", '{"id":"REQ123"}'

    monkeypatch.setattr(chip_atlas, "_http_request", fake_request)
    result = chip_atlas.submit_enrichment(
        regions=pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [100]}),
        genome="hg38",
    )
    assert result["request_id"] == "REQ123"
    assert captured["url"] == chip_atlas.DEFAULT_SUBMIT_URL
    payload = captured["data"].decode("utf-8")
    assert "format=text" in payload
    assert "result=www" in payload
    assert "typeA=bed" in payload
    assert "bedAFile=" in payload
    assert "typeB=rnd" in payload
    assert "bedBFile=empty" in payload
    assert "antigenClass=TFs+and+others" in payload
    assert "cellClass=No+description" in payload
    assert "threshold=100" in payload


def test_get_status_uses_request_path_endpoint(monkeypatch):
    captured = {}

    def fake_request(**kwargs):
        captured.update(kwargs)
        return 200, "text/plain", "current-state: running"

    monkeypatch.setattr(chip_atlas, "_http_request", fake_request)
    status = chip_atlas.get_status("REQ-PATH")
    assert status["status"] == "running"
    assert captured["url"] == f"{chip_atlas.DEFAULT_STATUS_URL}REQ-PATH"


def test_poll_request_collects_history_until_finished(monkeypatch):
    states = iter(
        [
            {"status": "running", "raw_status": "running", "status_code": 200},
            {"status": "finished", "raw_status": "finished", "status_code": 200},
        ]
    )

    def fake_get_status(*args, **kwargs):
        payload = next(states)
        return {"request_id": "REQ1", **payload}

    monkeypatch.setattr(chip_atlas, "get_status", fake_get_status)
    poll = chip_atlas.poll_request(
        "REQ1", poll_interval_seconds=0.0, timeout_seconds=5.0
    )
    assert poll["status"] == "finished"
    assert [step["status"] for step in poll["history"]] == ["running", "finished"]


def test_fetch_result_parses_tsv(monkeypatch):
    captured = {}

    def fake_request(**kwargs):
        captured.update(kwargs)
        return 200, "text/tab-separated-values", "target\tcount\nCTCF\t42\n"

    monkeypatch.setattr(chip_atlas, "_http_request", fake_request)
    result = chip_atlas.fetch_result("REQ2")
    assert isinstance(result, pd.DataFrame)
    assert result.loc[0, "target"] == "CTCF"
    assert int(result.loc[0, "count"]) == 42
    assert (
        captured["url"] == f"{chip_atlas.DEFAULT_RESULT_URL}REQ2?info=result&format=tsv"
    )


def test_run_enrichment_raises_on_terminal_failure(monkeypatch):
    monkeypatch.setattr(
        chip_atlas,
        "submit_enrichment",
        lambda **kwargs: {"request_id": "REQ-ERR", "query": {}, "submit_url": "u"},
    )
    monkeypatch.setattr(
        chip_atlas,
        "poll_request",
        lambda *args, **kwargs: {"status": "error", "history": []},
    )
    with pytest.raises(RuntimeError, match="status 'error'"):
        chip_atlas.run_enrichment(
            regions=["chr1:0-100"],
            wait=True,
            fetch_results=False,
            raise_on_failure=True,
        )


def test_submit_enrichment_runs_crossmap_when_regions_genome_differs(monkeypatch):
    captured = {}

    def fake_convert(**kwargs):
        captured.update(kwargs)
        return pd.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [10],
                "end": [60],
                "name": ["region_0"],
                "score": [0],
                "strand": ["+"],
            }
        )

    def fake_resolve(**kwargs):
        return "/tmp/hs1_to_hg38.over.chain.gz"

    def fake_request(**kwargs):
        return 200, "application/json", '{"id":"REQX"}'

    monkeypatch.setattr(chip_atlas, "_resolve_chain_file", fake_resolve)
    monkeypatch.setattr(chip_atlas, "convert_regions_with_crossmap", fake_convert)
    monkeypatch.setattr(chip_atlas, "_http_request", fake_request)

    result = chip_atlas.submit_enrichment(
        regions=["chr1:100-200,+"],
        genome="hg38",
        regions_genome="chm13",
    )

    assert result["request_id"] == "REQX"
    assert captured["source_genome"] == "hs1"
    assert captured["target_genome"] == "hg38"
    assert captured["chain_file"] == "/tmp/hs1_to_hg38.over.chain.gz"
    assert result["query"]["crossmap"]["applied"] is True


def test_cache_chain_files_uses_resolver(monkeypatch):
    calls = []

    def fake_resolve(**kwargs):
        calls.append((kwargs["source_genome"], kwargs["target_genome"]))
        return (
            f"/tmp/{kwargs['source_genome']}_to_{kwargs['target_genome']}.over.chain.gz"
        )

    monkeypatch.setattr(chip_atlas, "_resolve_chain_file", fake_resolve)
    cached = chip_atlas.cache_chain_files(
        source_genome="chm13", target_genomes=("hg38", "hg19")
    )

    assert calls == [("chm13", "hg38"), ("chm13", "hg19")]
    assert str(cached["hg38"]).endswith("chm13_to_hg38.over.chain.gz")
    assert str(cached["hg19"]).endswith("chm13_to_hg19.over.chain.gz")


def test_resolve_crossmap_executable_falls_back_to_crossmap(monkeypatch):
    def fake_which(cmd):
        return "/usr/local/bin/CrossMap" if cmd == "CrossMap" else None

    monkeypatch.setattr(chip_atlas.shutil, "which", fake_which)
    resolved = chip_atlas._resolve_crossmap_executable("CrossMap.py")
    assert resolved == "CrossMap"


def test_search_peak_datasets_filters_metadata_rows(monkeypatch, tmp_path):
    rows = [
        {
            "Experimental ID": "SRX1",
            "Genome assembly": "hg38",
            "Antigen class": "TFs and others",
            "Antigen": "CTCF",
            "Cell type class": "Blood",
            "Cell type": "K562",
            "Cell type description": "",
            "Processing logs": "",
            "Title": "CTCF K562",
            "Meta data": "",
            "BigWig": "bw1",
            "Peak-call (BED) (q < 1E-05)": "http://example.org/srx1.bed",
            "Peak-call (BigBed) (q < 1E-05)": "http://example.org/srx1.bb",
        },
        {
            "Experimental ID": "SRX2",
            "Genome assembly": "mm10",
            "Antigen class": "TFs and others",
            "Antigen": "CTCF",
            "Cell type class": "Blood",
            "Cell type": "MEL",
            "Cell type description": "",
            "Processing logs": "",
            "Title": "CTCF mouse",
            "Meta data": "",
            "BigWig": "bw2",
            "Peak-call (BED) (q < 1E-05)": "http://example.org/srx2.bed",
            "Peak-call (BigBed) (q < 1E-05)": "http://example.org/srx2.bb",
        },
    ]

    monkeypatch.setattr(
        chip_atlas,
        "_ensure_experiment_list_zip",
        lambda **kwargs: tmp_path / "fake.zip",
    )
    monkeypatch.setattr(chip_atlas, "_iter_experiment_rows", lambda path: iter(rows))

    found = chip_atlas.search_peak_datasets(
        antigen="CTCF",
        genome="hg38",
        cell_type="k562",
        threshold="05",
    )
    assert found.shape[0] == 1
    assert found.loc[0, "dataset_id"] == "SRX1"
    assert found.loc[0, "bed_url"] == "http://example.org/srx1.bed"


def test_download_peak_datasets_writes_variants_and_crossmapped(monkeypatch, tmp_path):
    datasets = pd.DataFrame(
        [
            {
                "dataset_id": "SRX3",
                "bed_url": "http://example.org/srx3.bed",
                "genome_assembly": "hg38",
            }
        ]
    )
    bed_payload = (
        b"chr1\t10\t20\tp1\t100\t.\n"
        b"chr1\t30\t40\tp2\t10\t.\n"
        b"chr1\t50\t60\tp3\t50\t.\n"
        b"chr1\t70\t80\tp4\t5\t.\n"
    )
    monkeypatch.setattr(chip_atlas, "_download_bytes", lambda **kwargs: bed_payload)
    monkeypatch.setattr(
        chip_atlas,
        "convert_regions_with_crossmap",
        lambda **kwargs: pd.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [100],
                "end": [200],
                "name": ["mapped"],
                "score": [1],
                "strand": ["+"],
            }
        ),
    )

    manifest = chip_atlas.download_peak_datasets(
        datasets=datasets,
        output_dir=tmp_path,
        include_complete_sorted=True,
        include_top_n=2,
        include_bottom_n=2,
        stratify="quartiles",
        crossmap_target_genome="hg19",
    )
    variants = set(manifest["variant"].tolist())
    assert "full_sorted" in variants
    assert "top_2" in variants
    assert "bottom_2" in variants
    assert "quantile_1_of_4" in variants
    assert "full_sorted_crossmapped" in variants
    assert manifest.loc[manifest["crossmapped"] == True].shape[0] > 0  # noqa: E712
