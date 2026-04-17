from pathlib import Path

import pytest

from dimelo import run_modkit


def _fake_capabilities() -> run_modkit.ModkitCapabilities:
    return run_modkit.ModkitCapabilities(
        executable="modkit",
        version_raw="modkit 0.6.1",
        version="0.6.1",
        version_tuple=(0, 6, 1),
        supports_mod_threshold=True,
        supports_mod_thresholds=False,
        supports_modified_bases=True,
        supports_force_allow_implicit=False,
        supports_extract_subcommands=True,
        extract_supports_reference_long=True,
        extract_supports_reference_short=True,
    )


def test_ensure_modkit_available_forwards_executable_override(monkeypatch):
    captured: dict[str, object] = {}
    capabilities = _fake_capabilities()

    def fake_get_modkit_capabilities(quiet=False, executable=None):
        captured["quiet"] = quiet
        captured["executable"] = executable
        return capabilities

    monkeypatch.setattr(run_modkit, "get_modkit_capabilities", fake_get_modkit_capabilities)

    result = run_modkit._ensure_modkit_available(
        quiet=True,
        executable=Path("/opt/modkit-0.6.1/bin/modkit"),
    )

    assert result == capabilities
    assert captured == {
        "quiet": True,
        "executable": "/opt/modkit-0.6.1/bin/modkit",
    }


def test_resolve_modkit_executable_prefers_env_override(monkeypatch):
    monkeypatch.setenv(run_modkit.MODKIT_EXECUTABLE_ENV, "modkit-0.6")
    monkeypatch.setattr(run_modkit.shutil, "which", lambda name: f"/usr/local/bin/{name}")

    resolved = run_modkit._resolve_modkit_executable(None)

    assert resolved == "/usr/local/bin/modkit-0.6"


def test_resolve_modkit_executable_raises_for_missing_path():
    with pytest.raises(FileNotFoundError, match="does not exist"):
        run_modkit._resolve_modkit_executable("/tmp/does-not-exist/modkit")


def test_configure_modkit_executable_sets_and_clears_env(monkeypatch):
    monkeypatch.delenv(run_modkit.MODKIT_EXECUTABLE_ENV, raising=False)

    run_modkit.configure_modkit_executable("/opt/modkit-custom/bin/modkit")
    assert (
        run_modkit.os.environ.get(run_modkit.MODKIT_EXECUTABLE_ENV)
        == "/opt/modkit-custom/bin/modkit"
    )

    run_modkit.configure_modkit_executable(None)
    assert run_modkit.MODKIT_EXECUTABLE_ENV not in run_modkit.os.environ


def test_modkit_cache_fingerprint_changes_when_binary_changes(tmp_path):
    fake_modkit = tmp_path / "modkit"
    fake_modkit.write_text("v1")
    fp1 = run_modkit._modkit_cache_fingerprint(str(fake_modkit))

    # Ensure metadata changes (mtime and size) to trigger new fingerprint.
    fake_modkit.write_text("v2-updated")
    fp2 = run_modkit._modkit_cache_fingerprint(str(fake_modkit))

    assert fp1 != fp2


def test_get_modkit_capabilities_uses_resolved_path_and_fingerprint(monkeypatch):
    captured: dict[str, object] = {}
    capabilities = _fake_capabilities()

    def fake_prepare_modkit_path(quiet=False):
        captured["prepared_quiet"] = quiet

    def fake_resolve_modkit_executable(executable):
        captured["resolved_from"] = executable
        return "/resolved/modkit"

    def fake_modkit_cache_fingerprint(executable_path):
        captured["fingerprint_path"] = executable_path
        return "fingerprint-123"

    def fake_cached(*, executable_path, executable_fingerprint, quiet=False):
        captured["cached_executable_path"] = executable_path
        captured["cached_executable_fingerprint"] = executable_fingerprint
        captured["cached_quiet"] = quiet
        return capabilities

    monkeypatch.setattr(run_modkit, "_prepare_modkit_path", fake_prepare_modkit_path)
    monkeypatch.setattr(run_modkit, "_resolve_modkit_executable", fake_resolve_modkit_executable)
    monkeypatch.setattr(run_modkit, "_modkit_cache_fingerprint", fake_modkit_cache_fingerprint)
    monkeypatch.setattr(run_modkit, "_get_modkit_capabilities_cached", fake_cached)

    result = run_modkit.get_modkit_capabilities(
        quiet=True,
        executable="modkit-0.6.1",
    )

    assert result == capabilities
    assert captured == {
        "prepared_quiet": True,
        "resolved_from": "modkit-0.6.1",
        "fingerprint_path": "/resolved/modkit",
        "cached_executable_path": "/resolved/modkit",
        "cached_executable_fingerprint": "fingerprint-123",
        "cached_quiet": True,
    }
