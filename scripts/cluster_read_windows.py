#!/usr/bin/env python3
"""Cluster DiMeLo single-read windows from one or more extract HDF5 files.

Example:
    python scripts/cluster_read_windows.py \
        --sample-sheet samples.tsv \
        --motif A,0 --regions sites.bed --window-size 1000 \
        --clusters 6 --output-dir read_clusters

All samples are pooled before scaling and clustering so cluster labels have the same
meaning across samples. The input HDF5 files must be outputs of parse_bam.extract().
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dimelo import cluster


def parse_sample(value: str) -> tuple[str, Path]:
    try:
        sample_id, path = value.split("=", 1)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("use SAMPLE_ID=PATH") from exc
    if not sample_id or not path:
        raise argparse.ArgumentTypeError("use SAMPLE_ID=PATH")
    return sample_id, Path(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--extract",
        action="append",
        type=parse_sample,
        metavar="SAMPLE_ID=PATH",
        help="Repeat for each sample; clusters are learned jointly.",
    )
    parser.add_argument(
        "--sample-sheet",
        type=Path,
        help="TSV with sample_id and extract_h5 columns; other columns join assignments.",
    )
    parser.add_argument("--motif", action="append", default=None)
    parser.add_argument("--regions", type=Path, required=True)
    parser.add_argument(
        "--window-size",
        type=int,
        default=1000,
        help="Half-window in bp; 1000 produces a 2-kb window (default: 1000).",
    )
    parser.add_argument(
        "--low-memory",
        action="store_true",
        help="Load regions in batches during read-window extraction.",
    )
    parser.add_argument(
        "--region-batch-size",
        type=int,
        default=64,
        help="Number of regions per extraction batch when --low-memory is set.",
    )
    parser.add_argument("--clusters", type=int, default=6)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--output-dir", type=Path, default=Path("read_clusters"))
    return parser.parse_args()


def clean_for_csv(frame: pd.DataFrame) -> pd.DataFrame:
    cleaned = frame.copy()
    for column in cleaned.columns:
        if pd.api.types.is_object_dtype(cleaned[column].dtype) or pd.api.types.is_string_dtype(
            cleaned[column].dtype
        ):
            cleaned[column] = cleaned[column].map(
                lambda value: value.decode() if isinstance(value, bytes) else value
            )
    return cleaned


def load_samples(args: argparse.Namespace) -> pd.DataFrame:
    frames = []
    if args.sample_sheet is not None:
        sheet = pd.read_csv(args.sample_sheet, sep="\t", dtype={"sample_id": str})
        required = {"sample_id", "extract_h5"}
        missing = sorted(required - set(sheet.columns))
        if missing:
            raise ValueError(f"Sample sheet is missing columns: {', '.join(missing)}")
        frames.append(sheet)
    if args.extract:
        frames.append(
            pd.DataFrame(
                [{"sample_id": sample_id, "extract_h5": str(path)} for sample_id, path in args.extract]
            )
        )
    if not frames:
        raise ValueError("Provide --sample-sheet or at least one --extract SAMPLE_ID=PATH.")
    samples = pd.concat(frames, ignore_index=True)
    if samples["sample_id"].duplicated().any():
        duplicates = sorted(samples.loc[samples["sample_id"].duplicated(), "sample_id"].unique())
        raise ValueError(f"Sample IDs must be unique; duplicates: {', '.join(duplicates)}")
    return samples


def main() -> None:
    args = parse_args()
    motifs = args.motif or ["A,0"]
    samples = load_samples(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.low_memory and len(motifs) > 1:
        raise ValueError(
            "--low-memory is currently supported only for single-motif extraction; "
            "multi-motif clustering uses build_multimotif_read_windows()."
        )

    results = []
    sample_ids = []
    sample_metadata = samples.set_index("sample_id").drop(columns="extract_h5").to_dict("index")
    for sample in samples.to_dict("records"):
        sample_id = str(sample["sample_id"])
        h5_path = Path(sample["extract_h5"])
        if not h5_path.is_file():
            raise FileNotFoundError(f"Missing extract HDF5 for {sample_id}: {h5_path}")
        if len(motifs) == 1:
            result = cluster.extract_read_windows(
                hdf5_file=h5_path,
                motifs=motifs,
                regions=args.regions,
                config=cluster.ReadWindowExtractionConfig(
                    window_size=args.window_size,
                    orientation_aware=True,
                    filter_multi_region_reads=True,
                ),
                span_full_window=True,
                low_memory=args.low_memory,
                region_batch_size=args.region_batch_size,
                random_state=args.random_state,
            )
        else:
            result = cluster.build_multimotif_read_windows(
                hdf5_file=h5_path,
                motifs=motifs,
                regions=args.regions,
                window_size=args.window_size,
                orientation_aware=True,
                span_full_window=True,
                require_all_motifs=True,
            )
        results.append(result)
        sample_ids.append(sample_id)

    reads = cluster.merge_read_window_results(results, source_labels=sample_ids)
    features, feature_names, feature_table, scale_table = (
        cluster.read_window_feature_matrix(
            reads,
            scale_features=True,
            scaling_method="robust",
            family_weighting="equal_family",
            return_feature_table=True,
            return_scale_table=True,
        )
    )
    if features.shape[0] < args.clusters:
        raise ValueError(
            f"Only {features.shape[0]} reads passed filtering; cannot fit "
            f"{args.clusters} clusters."
        )

    fit = cluster.cluster_read_windows(
        features,
        method="minibatch_kmeans",
        n_clusters=args.clusters,
        random_state=args.random_state,
        n_init=20,
        batch_size=min(4096, max(256, features.shape[0])),
    )
    labels = fit.labels_size_ordered

    assignments = clean_for_csv(pd.DataFrame(reads.metadata))
    assignments.insert(0, "sample_id", assignments.pop("source_label"))
    for column in samples.columns:
        if column not in {"sample_id", "extract_h5"} and column not in assignments.columns:
            assignments[column] = assignments["sample_id"].map(
                {sample_id: values.get(column) for sample_id, values in sample_metadata.items()}
            )
    assignments["cluster"] = labels
    assignments.to_csv(args.output_dir / "read_cluster_assignments.tsv", sep="\t", index=False)

    sample_group_columns = ["sample_id"]
    if "condition" in assignments.columns:
        sample_group_columns.append("condition")
    counts = (
        assignments.groupby([*sample_group_columns, "cluster"], observed=True)
        .size()
        .rename("read_count")
        .reset_index()
    )
    totals = counts.groupby("sample_id")["read_count"].transform("sum")
    counts["read_fraction"] = counts["read_count"] / totals
    counts.to_csv(args.output_dir / "sample_cluster_fractions.tsv", sep="\t", index=False)

    if "condition" in assignments.columns:
        condition_counts = (
            counts.groupby(["condition", "cluster"], observed=True)
            .agg(read_count=("read_count", "sum"), replicate_n=("sample_id", "nunique"))
            .reset_index()
        )
        condition_totals = condition_counts.groupby("condition")["read_count"].transform("sum")
        condition_counts["read_fraction"] = condition_counts["read_count"] / condition_totals
        condition_counts.to_csv(
            args.output_dir / "condition_cluster_fractions.tsv", sep="\t", index=False
        )

    region_columns = ["chromosome", "region_start", "region_end", "region_strand"]
    if set(region_columns).issubset(assignments.columns):
        region_group_columns = [*region_columns, *sample_group_columns]
        region_counts = (
            assignments.dropna(subset=["chromosome", "region_start", "region_end"])
            .groupby([*region_group_columns, "cluster"], observed=True, dropna=False)
            .size()
            .rename("read_count")
            .reset_index()
        )
        region_totals = region_counts.groupby(region_group_columns, dropna=False)[
            "read_count"
        ].transform("sum")
        region_counts["read_fraction"] = region_counts["read_count"] / region_totals
        region_counts.to_csv(
            args.output_dir / "region_cluster_occupancy.tsv", sep="\t", index=False
        )

    profiles = []
    if len(motifs) == 1:
        motif_width = reads.data_matrix.shape[1]
    elif reads.data_matrix.shape[1] % len(motifs) != 0:
        raise ValueError("Read-window width is not divisible by the number of motifs.")
    else:
        motif_width = reads.data_matrix.shape[1] // len(motifs)
    x = np.arange(motif_width) - motif_width // 2
    for label in np.unique(labels):
        cluster_mean = reads.data_matrix[labels == label].mean(axis=0)
        for motif_index, motif in enumerate(motifs):
            start = motif_index * motif_width
            stop = start + motif_width
            profiles.append(
                pd.DataFrame(
                    {
                        "position_bp": x,
                        "cluster": label,
                        "motif": motif,
                        "mean_signal": cluster_mean[start:stop],
                    }
                )
            )
    profile_table = pd.concat(profiles, ignore_index=True)
    profile_table.to_csv(args.output_dir / "cluster_profiles.tsv", sep="\t", index=False)
    feature_table.to_csv(args.output_dir / "feature_manifest.tsv", sep="\t", index=False)
    scale_table.to_csv(args.output_dir / "feature_scaling.tsv", sep="\t", index=False)

    qc = {
        "samples": samples.to_dict("records"),
        "motifs": motifs,
        "window_size_half_bp": args.window_size,
        "n_reads": int(features.shape[0]),
        "n_features": int(features.shape[1]),
        "n_clusters": args.clusters,
        "random_state": args.random_state,
        "low_memory": bool(args.low_memory),
        "region_batch_size": int(args.region_batch_size),
        "assignment_level": "read_window",
        "multi_motif": len(motifs) > 1,
        "metrics": {key: None if value is None else float(value) for key, value in fit.metrics.items()},
        "feature_names": feature_names,
    }
    (args.output_dir / "clustering_qc.json").write_text(json.dumps(qc, indent=2) + "\n")

    pivot = counts.pivot(index="sample_id", columns="cluster", values="read_fraction").fillna(0)
    ax = pivot.plot.bar(stacked=True, figsize=(max(6, len(pivot) * 1.2), 4))
    ax.set_ylabel("Fraction of reads")
    ax.legend(title="Cluster", bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.figure.tight_layout()
    ax.figure.savefig(args.output_dir / "sample_cluster_fractions.png", dpi=200)
    plt.close(ax.figure)

    fig, ax = plt.subplots(figsize=(8, 4))
    for (label, motif), group in profile_table.groupby(["cluster", "motif"], observed=True):
        ax.plot(
            group["position_bp"],
            group["mean_signal"],
            label=f"cluster {label}, {motif}",
        )
    ax.set(xlabel="Position relative to region center (bp)", ylabel="Mean thresholded signal")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(args.output_dir / "cluster_profiles.png", dpi=200)
    plt.close(fig)

    print(f"Wrote {features.shape[0]:,} read assignments to {args.output_dir}")


if __name__ == "__main__":
    main()
