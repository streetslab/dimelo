#!/usr/bin/env python3
"""
Verify that the active Python/kernel environment can run the DiMeLo notebooks.

Usage:
  python scripts/ensure_dimelo_kernel.py
  python scripts/ensure_dimelo_kernel.py --fix
  python scripts/ensure_dimelo_kernel.py --fix --skip-modkit
  python scripts/ensure_dimelo_kernel.py --expected-env dimelo-toolkit
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODKIT_REQUIREMENT = "supported"
DEFAULT_MODKIT_INSTALL_VERSION = "0.2.4"
SUPPORTED_MODKIT_PREFIXES = ("0.2.4", "0.6.")

# setup.py install_requires
CORE_MODULES = {
    "numpy": "numpy",
    "seaborn": "seaborn",
    "pysam": "pysam",
    "h5py": "h5py",
    "pyBigWig": "pyBigWig",
    "notebook": "notebook",
    "ipykernel": "ipykernel",
    "ipywidgets": "ipywidgets",
    "tqdm": "tqdm",
    "plotly": "plotly",
    "kaleido": "kaleido",
}

# setup.py extras_require["clustering"]
CLUSTERING_MODULES = {
    "scikit-learn": "sklearn",
    "scipy": "scipy",
    "hdbscan": "hdbscan",
    "umap-learn": "umap",
    "pyranges": "pyranges",
    "xgboost": "xgboost",
}


@dataclass
class CheckReport:
    missing_core: list[str]
    missing_clustering: list[str]
    runtime_failures: dict[str, str]
    missing_modkit: bool
    modkit_version: str | None
    missing_crossmap: bool
    crossmap_executable: str | None
    conda_env: str | None
    expected_env: str | None

    @property
    def ok(self) -> bool:
        env_ok = self.expected_env is None or self.conda_env == self.expected_env
        return (
            len(self.missing_core) == 0
            and len(self.missing_clustering) == 0
            and len(self.runtime_failures) == 0
            and not self.missing_modkit
            and not self.missing_crossmap
            and env_ok
        )


def _module_installed(module_name: str) -> bool:
    return importlib.util.find_spec(module_name) is not None


def _check_python_modules() -> tuple[list[str], list[str]]:
    missing_core = [
        package for package, module in CORE_MODULES.items() if not _module_installed(module)
    ]
    missing_clustering = [
        package
        for package, module in CLUSTERING_MODULES.items()
        if not _module_installed(module)
    ]
    return missing_core, missing_clustering


def _modkit_version_ok(version_text: str, required_version: str) -> bool:
    normalized = required_version.strip().lower()
    if normalized in {"supported", "dual", "auto"}:
        return any(prefix in version_text for prefix in SUPPORTED_MODKIT_PREFIXES)
    return required_version in version_text


def _check_modkit(required_version: str) -> tuple[bool, str | None]:
    modkit_path = shutil.which("modkit")
    if modkit_path is None:
        return True, None

    try:
        result = subprocess.run(
            ["modkit", "--version"], check=True, capture_output=True, text=True
        )
    except Exception:
        return True, None

    version_text = (result.stdout or result.stderr).strip().splitlines()[0]
    version_ok = _modkit_version_ok(version_text, required_version)
    return (not version_ok), version_text


def _check_crossmap() -> tuple[bool, str | None]:
    for candidate in ("CrossMap.py", "CrossMap"):
        path = shutil.which(candidate)
        if path is not None:
            return False, path
    return True, None


def _check_runtime_imports() -> dict[str, str]:
    runtime_failures: dict[str, str] = {}
    runtime_modules = ("numpy", "pysam", "h5py", "pyBigWig")
    for module_name in runtime_modules:
        try:
            __import__(module_name)
        except Exception as exc:  # pragma: no cover - exception type varies by platform
            runtime_failures[module_name] = f"{type(exc).__name__}: {exc}"
    return runtime_failures


def run_checks(
    required_modkit_version: str,
    skip_modkit: bool,
    expected_env: str | None,
) -> CheckReport:
    missing_core, missing_clustering = _check_python_modules()
    runtime_failures = _check_runtime_imports()
    if skip_modkit:
        missing_modkit = False
        modkit_version = None
    else:
        missing_modkit, modkit_version = _check_modkit(required_modkit_version)
    missing_crossmap, crossmap_executable = _check_crossmap()

    return CheckReport(
        missing_core=missing_core,
        missing_clustering=missing_clustering,
        runtime_failures=runtime_failures,
        missing_modkit=missing_modkit,
        modkit_version=modkit_version,
        missing_crossmap=missing_crossmap,
        crossmap_executable=crossmap_executable,
        conda_env=os.environ.get("CONDA_DEFAULT_ENV"),
        expected_env=expected_env,
    )


def _run_command(cmd: list[str], cwd: Path | None = None) -> None:
    print(f"$ {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=cwd)


def fix_python_packages(report: CheckReport) -> None:
    # Install core and clustering extras in one shot; safe if already installed.
    if report.missing_core or report.missing_clustering:
        _run_command(
            [sys.executable, "-m", "pip", "install", "-e", ".[clustering]"],
            cwd=REPO_ROOT,
        )


def fix_modkit(required_modkit_version: str) -> None:
    requested_version = (
        DEFAULT_MODKIT_INSTALL_VERSION
        if required_modkit_version.strip().lower() in {"supported", "dual", "auto"}
        else required_modkit_version
    )
    conda_path = shutil.which("conda")
    if conda_path is None:
        print(
            "modkit is missing/incompatible and `conda` is not on PATH.\n"
            "Install manually with: conda install -y nanoporetech::modkit=="
            f"{requested_version}"
        )
        return

    _run_command(
        ["conda", "install", "-y", f"nanoporetech::modkit=={requested_version}"]
    )


def fix_crossmap() -> None:
    conda_path = shutil.which("conda")
    if conda_path is None:
        print(
            "CrossMap is missing and `conda` is not on PATH.\n"
            "Install manually with: conda install -y bioconda::crossmap=0.7.3"
        )
        return
    _run_command(["conda", "install", "-y", "bioconda::crossmap=0.7.3"])


def print_report(
    report: CheckReport,
    required_modkit_version: str,
    skip_modkit: bool,
) -> None:
    print("=== DiMeLo Kernel Environment Check ===")
    print(f"Python executable: {sys.executable}")
    print(f"Conda env: {report.conda_env or '(not set)'}")
    if report.expected_env is not None:
        print(f"Expected env: {report.expected_env}")
    print(f"Repo root: {REPO_ROOT}")
    print("")

    if report.missing_core:
        print("Missing core Python packages:")
        for pkg in report.missing_core:
            print(f"  - {pkg}")
    else:
        print("Core Python packages: OK")

    if report.missing_clustering:
        print("Missing clustering Python packages:")
        for pkg in report.missing_clustering:
            print(f"  - {pkg}")
    else:
        print("Clustering Python packages: OK")

    if report.runtime_failures:
        print("Runtime import failures:")
        for module_name, message in report.runtime_failures.items():
            print(f"  - {module_name}: {message}")
    else:
        print("Runtime imports: OK")

    if report.expected_env is not None:
        if report.conda_env == report.expected_env:
            print("Conda env match: OK")
        else:
            print(
                "Conda env match: FAIL "
                f"(active='{report.conda_env}', expected='{report.expected_env}')"
            )

    if skip_modkit:
        print("modkit check: skipped")
    else:
        requirement_text = (
            "any supported version (0.2.4 or 0.6.x)"
            if required_modkit_version.strip().lower() in {"supported", "dual", "auto"}
            else f"version containing '{required_modkit_version}'"
        )
        if report.missing_modkit:
            if report.modkit_version is None:
                print(f"modkit: MISSING (required {requirement_text})")
            else:
                print(
                    "modkit: version mismatch "
                    f"(found '{report.modkit_version}', required {requirement_text})"
                )
        else:
            print(f"modkit: OK ({report.modkit_version})")

    if report.missing_crossmap:
        print("CrossMap: MISSING (expected CrossMap.py or CrossMap on PATH)")
    else:
        print(f"CrossMap: OK ({report.crossmap_executable})")

    print("")
    print("Overall status:", "PASS" if report.ok else "FAIL")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fix",
        action="store_true",
        help="Attempt to install/fix missing dependencies automatically.",
    )
    parser.add_argument(
        "--skip-modkit",
        action="store_true",
        help="Skip modkit checks and install attempts.",
    )
    parser.add_argument(
        "--modkit-version",
        default=DEFAULT_MODKIT_REQUIREMENT,
        help=(
            "Required modkit version marker. "
            "Use 'supported' (default) to allow 0.2.4 or 0.6.x, or pass an explicit marker."
        ),
    )
    parser.add_argument(
        "--expected-env",
        default=None,
        help=(
            "Expected active CONDA_DEFAULT_ENV name (e.g. dimelo-toolkit). "
            "When provided, mismatch is treated as failure."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    report = run_checks(
        required_modkit_version=args.modkit_version,
        skip_modkit=args.skip_modkit,
        expected_env=args.expected_env,
    )
    print_report(
        report,
        required_modkit_version=args.modkit_version,
        skip_modkit=args.skip_modkit,
    )

    if report.ok:
        return 0

    if not args.fix:
        print("\nRun with `--fix` to install missing dependencies.")
        if report.expected_env and report.conda_env != report.expected_env:
            print(
                f"Activate expected env first: `conda activate {report.expected_env}`"
            )
        if report.runtime_failures:
            print("For runtime linker/import issues, run scripts/bootstrap_dimelo_env.sh")
        return 1

    try:
        fix_python_packages(report)
        if not args.skip_modkit and report.missing_modkit:
            fix_modkit(args.modkit_version)
        if report.missing_crossmap:
            fix_crossmap()
    except subprocess.CalledProcessError as exc:
        print(f"\nInstall step failed: {exc}")
        return exc.returncode or 1

    print("\nRe-checking environment after attempted fixes...\n")
    post = run_checks(
        required_modkit_version=args.modkit_version,
        skip_modkit=args.skip_modkit,
        expected_env=args.expected_env,
    )
    print_report(
        post,
        required_modkit_version=args.modkit_version,
        skip_modkit=args.skip_modkit,
    )
    return 0 if post.ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
