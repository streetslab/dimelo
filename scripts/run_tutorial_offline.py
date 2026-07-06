#!/usr/bin/env python3
"""
Execute the deterministic offline tutorial notebook and capture logs/artifacts.
"""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_NOTEBOOK = REPO_ROOT / "tutorial_offline.ipynb"
DEFAULT_ARTIFACT_DIR = REPO_ROOT / "artifacts" / "tutorial_offline"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--notebook",
        type=Path,
        default=DEFAULT_NOTEBOOK,
        help=f"Notebook path to execute (default: {DEFAULT_NOTEBOOK}).",
    )
    parser.add_argument(
        "--artifact-dir",
        type=Path,
        default=DEFAULT_ARTIFACT_DIR,
        help=f"Directory for executed notebook and logs (default: {DEFAULT_ARTIFACT_DIR}).",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=0,
        help="Per-cell timeout for nbconvert execution (0 means no timeout).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    notebook = args.notebook.resolve()
    artifact_dir = args.artifact_dir.resolve()
    artifact_dir.mkdir(parents=True, exist_ok=True)

    if not notebook.exists():
        raise FileNotFoundError(f"Notebook not found: {notebook}")

    executed_name = notebook.with_suffix(".executed.ipynb").name
    log_path = artifact_dir / "tutorial_offline_execution.log"

    command = [
        "python3.11",
        "-m",
        "jupyter",
        "nbconvert",
        "--to",
        "notebook",
        "--execute",
        str(notebook),
        "--output",
        executed_name,
        "--output-dir",
        str(artifact_dir),
        "--ExecutePreprocessor.timeout",
        str(args.timeout_seconds),
    ]

    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")

    result = subprocess.run(
        command,
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
    )
    combined_log = "\n".join(
        part.strip() for part in (result.stdout, result.stderr) if part and part.strip()
    )
    log_path.write_text(combined_log + ("\n" if combined_log else ""))

    print(f"Wrote execution log to {log_path}")
    print(f"Executed notebook target: {artifact_dir / executed_name}")

    if result.returncode != 0:
        print(combined_log)
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
