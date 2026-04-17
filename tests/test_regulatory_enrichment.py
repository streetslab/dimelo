from pathlib import Path

import pytest

from dimelo import regulatory_enrichment, workflows
from dimelo.models import UniBindJobResult


def test_regulatory_spec_defaults_to_human_and_enables_screen_and_unibind():
    spec = regulatory_enrichment.RegulatoryEnrichmentSpec()

    assert spec.species == "homo_sapiens"
    assert spec.target_genome == "hg38"
    assert spec.enabled_providers == ("screen", "unibind")
    assert spec.provider_notes == {}


def test_regulatory_spec_infers_mouse_from_target_genome():
    spec = regulatory_enrichment.RegulatoryEnrichmentSpec(
        species=None,
        target_genome="mm10",
    )

    assert spec.species == "mus_musculus"
    assert spec.enabled_providers == ("screen", "unibind")


def test_regulatory_spec_filters_screen_for_non_screen_species():
    spec = regulatory_enrichment.RegulatoryEnrichmentSpec(
        species="rat",
        providers=("screen", "unibind"),
    )

    assert spec.species == "rattus_norvegicus"
    assert spec.enabled_providers == ("unibind",)
    assert "screen" in spec.provider_notes
    assert "supports only homo_sapiens and mus_musculus" in spec.provider_notes["screen"]


def test_regulatory_spec_can_raise_on_filtered_screen_when_strict():
    with pytest.raises(regulatory_enrichment.RegulatoryEnrichmentSpecError, match="SCREEN disabled"):
        regulatory_enrichment.RegulatoryEnrichmentSpec(
            species="rattus norvegicus",
            providers=("screen",),
            strict_provider_support=True,
        )


def test_regulatory_spec_rejects_species_not_supported_by_unibind():
    with pytest.raises(regulatory_enrichment.RegulatoryEnrichmentSpecError, match="Unsupported species"):
        regulatory_enrichment.RegulatoryEnrichmentSpec(
            species="bos_taurus",
            providers=("unibind",),
        )


def test_regulatory_spec_preserves_arbitrary_genome_and_chain_fields():
    spec = regulatory_enrichment.RegulatoryEnrichmentSpec(
        species="human",
        providers=("unibind",),
        reference_genome="custom_build_v2",
        target_genome="custom_build_v3",
        crossmap_chain_file=Path("cache/chains/custom.over.chain.gz"),
        crossmap_chain_url="https://example.org/custom.chain.gz",
        crossmap_chain_cache_dir=Path("cache/chains/custom"),
        crossmap_executable="CrossMap",
    )

    serialized = spec.as_dict()
    assert serialized["reference_genome"] == "custom_build_v2"
    assert serialized["target_genome"] == "custom_build_v3"
    assert serialized["crossmap_chain_file"] == "cache/chains/custom.over.chain.gz"
    assert serialized["crossmap_chain_url"] == "https://example.org/custom.chain.gz"
    assert serialized["crossmap_chain_cache_dir"] == "cache/chains/custom"
    assert serialized["crossmap_executable"] == "CrossMap"


def test_workflow_wrapper_returns_regulatory_spec():
    spec = workflows.resolve_regulatory_enrichment_spec(
        species="mouse",
        providers=("screen", "unibind"),
        reference_genome="mm39",
    )

    assert isinstance(spec, regulatory_enrichment.RegulatoryEnrichmentSpec)
    assert spec.species == "mus_musculus"
    assert spec.enabled_providers == ("screen", "unibind")


def test_provider_alias_normalization_and_deduplication():
    spec = regulatory_enrichment.RegulatoryEnrichmentSpec(
        species="human",
        providers=("SCREEN2.0", "screen", "UniBind2021", "unibind"),
    )

    assert spec.providers == ("screen", "unibind")
    assert spec.enabled_providers == ("screen", "unibind")


def test_unibind_trackhub_url_supports_named_collections():
    robust = regulatory_enrichment.unibind_trackhub_url("robust")
    permissive = regulatory_enrichment.unibind_trackhub_url("permissive")
    assert robust.endswith("/UniBind_hubs_Robust/UCSC/hub.txt")
    assert permissive.endswith("/UniBind_hubs_Permissive/UCSC/hub.txt")


def test_search_unibind_trackhub_tracks_and_resolve_to_cache(monkeypatch, tmp_path):
    hub_url = "https://example.org/hub.txt"
    genomes_url = "https://example.org/genomes.txt"
    trackdb_url = "https://example.org/hg38/trackDb.txt"
    bigbed_url = "https://example.org/hg38/CTCF_track.bb"
    other_url = "https://example.org/hg38/ATAC_track.bb"

    payloads = {
        hub_url: (
            "hub test\n"
            "shortLabel Test\n"
            "longLabel Test\n"
            "genomesFile genomes.txt\n"
            "email test@example.org\n"
        ).encode(),
        genomes_url: (
            "genome hg38\n"
            "trackDb hg38/trackDb.txt\n\n"
            "genome mm10\n"
            "trackDb mm10/trackDb.txt\n"
        ).encode(),
        trackdb_url: (
            "track CTCF_track\n"
            "shortLabel CTCF\n"
            "longLabel CTCF example\n"
            "type bigBed 9 .\n"
            "bigDataUrl CTCF_track.bb\n\n"
            "track ATAC_track\n"
            "shortLabel ATAC\n"
            "longLabel ATAC example\n"
            "type bigBed 9 .\n"
            "bigDataUrl ATAC_track.bb\n"
        ).encode(),
        bigbed_url: b"fake-ctcf-track",
        other_url: b"fake-atac-track",
    }

    class FakeResponse:
        def __init__(self, payload: bytes) -> None:
            self._payload = payload

        def read(self) -> bytes:
            return self._payload

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    def fake_urlopen(req, timeout=60.0):  # noqa: ARG001
        url = req.full_url if hasattr(req, "full_url") else str(req)
        if url not in payloads:
            raise RuntimeError(f"unexpected url: {url}")
        return FakeResponse(payloads[url])

    monkeypatch.setattr(regulatory_enrichment.request, "urlopen", fake_urlopen)

    rows = regulatory_enrichment.search_unibind_trackhub_tracks(
        trackhub_url=hub_url,
        assembly="hg38",
        search_terms=["ctcf"],
    )
    assert len(rows) == 1
    assert rows[0]["track"] == "CTCF_track"
    assert rows[0]["url"] == bigbed_url

    cached = regulatory_enrichment.resolve_unibind_track_paths(
        trackhub_url=hub_url,
        assembly="hg38",
        search_terms=["ctcf"],
        cache_dir=tmp_path,
    )
    assert len(cached) == 1
    assert cached[0].exists()
    assert cached[0].read_bytes() == b"fake-ctcf-track"


def test_resolve_unibind_track_paths_prefers_explicit_local_paths(tmp_path):
    p = tmp_path / "local_track.bed"
    p.write_text("chr1\\t0\\t10\\n", encoding="utf-8")

    resolved = regulatory_enrichment.resolve_unibind_track_paths(
        track_paths=[p],
    )
    assert len(resolved) == 1
    assert resolved[0] == p.resolve()


def test_workflow_wrapper_resolve_unibind_track_paths_uses_regulatory_spec(monkeypatch):
    captured = {}

    def fake_resolve(**kwargs):
        captured.update(kwargs)
        return [Path("/tmp/fake.bb")]

    monkeypatch.setattr(regulatory_enrichment, "resolve_unibind_track_paths", fake_resolve)
    spec = workflows.resolve_regulatory_enrichment_spec(
        species="human",
        target_genome="hg38",
        providers=("unibind",),
    )

    paths = workflows.resolve_unibind_track_paths(
        search_terms=["ctcf"],
        regulatory_spec=spec,
    )
    assert len(paths) == 1
    assert captured["assembly"] == "hg38"


def test_submit_unibind_tfbs_extraction_parses_job_page(monkeypatch):
    endpoint = "https://example.org/TFBS_extraction/"
    landing_html = (
        '<form><input type="hidden" name="csrfmiddlewaretoken" value="token123"></form>'
    )
    submitted_html = (
        "<table>"
        "<tr><th>Job status:</th><td>Queued</td></tr>"
        "<tr><th>Results URL:</th><td><a href=\"https://example.org/TFBS_extraction/jobA/\">"
        "https://example.org/TFBS_extraction/jobA/"
        "</a></td></tr>"
        "</table>"
    )

    class FakeResponse:
        def __init__(self, payload: str, url: str) -> None:
            self._payload = payload.encode("utf-8")
            self._url = url

        def read(self) -> bytes:
            return self._payload

        def geturl(self) -> str:
            return self._url

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    class FakeOpener:
        def open(self, req, timeout=60.0):  # noqa: ARG002
            url = req.full_url if hasattr(req, "full_url") else str(req)
            data = getattr(req, "data", None)
            if url == endpoint and data is None:
                return FakeResponse(landing_html, endpoint)
            if url == endpoint and data is not None:
                assert b"performTFBSextraction" in data
                assert b"collection" in data
                return FakeResponse(submitted_html, "https://example.org/TFBS_extraction/jobA/")
            raise RuntimeError(f"Unexpected request url={url!r}")

    monkeypatch.setattr(regulatory_enrichment.request, "build_opener", lambda *args: FakeOpener())

    result = regulatory_enrichment.submit_unibind_tfbs_extraction(
        regions=["chr1:10-20"],
        species="human",
        collection="robust",
        endpoint_url=endpoint,
    )
    assert isinstance(result, UniBindJobResult)
    assert result.job_id == "jobA"
    assert result.status == "queued"
    assert result.job_url == "https://example.org/TFBS_extraction/jobA/"
    assert result.query["species"] == "homo_sapiens"


def test_submit_unibind_enrichment_requires_required_inputs():
    with pytest.raises(ValueError, match="background_regions"):
        regulatory_enrichment.submit_unibind_enrichment(
            regions=["chr1:10-20"],
            analysis_type="oneSetBg",
            endpoint_url="https://example.org/enrichment/",
        )

    with pytest.raises(ValueError, match="comparison_regions"):
        regulatory_enrichment.submit_unibind_enrichment(
            regions=["chr1:10-20"],
            analysis_type="twoSets",
            endpoint_url="https://example.org/enrichment/",
        )


def test_poll_unibind_job_waits_until_terminal(monkeypatch):
    pages = [
        "<tr><th>Job status:</th><td>Running</td></tr>",
        (
            "<tr><th>Job status:</th><td>Completed</td></tr>"
            '<a href="/temp/jobB/extraction_results.bed">Download</a>'
        ),
    ]
    state = {"i": 0}

    def fake_fetch_text(url: str, timeout_seconds: float = 60.0):  # noqa: ARG001
        index = state["i"]
        if index < len(pages) - 1:
            state["i"] += 1
        return pages[index]

    monkeypatch.setattr(regulatory_enrichment, "_fetch_text", fake_fetch_text)

    start = UniBindJobResult(
        job_id="jobB",
        status="queued",
        job_url="https://example.org/TFBS_extraction/jobB/",
        endpoint_url="https://example.org/TFBS_extraction/",
    )
    result = regulatory_enrichment.poll_unibind_job(
        start,
        poll_interval_seconds=0.01,
        timeout_seconds=1.0,
    )
    assert result.status == "completed"
    assert "https://example.org/temp/jobB/extraction_results.bed" in result.download_urls
    assert len(result.status_history) >= 2


def test_download_unibind_job_outputs_writes_files(monkeypatch, tmp_path):
    payloads = {
        "https://example.org/temp/jobC/a.tsv": b"a",
        "https://example.org/temp/jobC/b.bed": b"b",
    }

    def fake_fetch_bytes(url: str, timeout_seconds: float = 60.0):  # noqa: ARG001
        return payloads[url]

    monkeypatch.setattr(regulatory_enrichment, "_fetch_bytes", fake_fetch_bytes)
    job = UniBindJobResult(
        job_id="jobC",
        status="completed",
        job_url="https://example.org/enrichment/jobC/",
        endpoint_url="https://example.org/enrichment/",
        download_urls=list(payloads.keys()),
    )
    outputs = regulatory_enrichment.download_unibind_job_outputs(
        job=job,
        output_dir=tmp_path,
    )
    assert len(outputs) == 2
    assert outputs[0].exists()
    assert outputs[1].exists()


def test_workflow_wrapper_unibind_tfbs_uses_regulatory_spec(monkeypatch):
    captured = {}

    def fake_run(**kwargs):
        captured.update(kwargs)
        return UniBindJobResult(
            job_id="jobD",
            status="completed",
            job_url="https://example.org/TFBS_extraction/jobD/",
            endpoint_url="https://example.org/TFBS_extraction/",
        )

    monkeypatch.setattr(regulatory_enrichment, "run_unibind_tfbs_extraction", fake_run)
    spec = workflows.resolve_regulatory_enrichment_spec(
        species="mouse",
        providers=("unibind",),
    )
    _ = workflows.unibind_tfbs_extraction_workflow(
        regions=["chr1:1-2"],
        regulatory_spec=spec,
        wait=False,
    )
    assert captured["species"] == "mus_musculus"


def test_workflow_wrapper_unibind_enrichment_uses_regulatory_spec(monkeypatch):
    captured = {}

    def fake_run(**kwargs):
        captured.update(kwargs)
        return UniBindJobResult(
            job_id="jobE",
            status="completed",
            job_url="https://example.org/enrichment/jobE/",
            endpoint_url="https://example.org/enrichment/",
        )

    monkeypatch.setattr(regulatory_enrichment, "run_unibind_enrichment", fake_run)
    spec = workflows.resolve_regulatory_enrichment_spec(
        species="human",
        providers=("unibind",),
    )
    _ = workflows.unibind_enrichment_workflow(
        regions=["chr1:1-2"],
        background_regions=["chr1:1-3"],
        regulatory_spec=spec,
        wait=False,
    )
    assert captured["species"] == "homo_sapiens"
