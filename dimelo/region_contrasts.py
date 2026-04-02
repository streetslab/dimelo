from __future__ import annotations

from functools import partial
import math
from typing import Iterable

import pandas as pd

from . import load_processed, utils
from .models import ContrastSpec, RegionContrastResult


def validate_region_contrast_request(
    *,
    analysis_unit: str,
    representation: str,
    signal_source: str,
    test: str,
) -> None:
    if analysis_unit == "ensemble_region":
        if signal_source != "pileup_counts":
            raise ValueError(
                "V1 region_contrasts inference requires signal_source='pileup_counts'."
            )
        if representation not in {"modified_fraction", "modified_count"}:
            raise ValueError(
                "V1 region_contrasts inference requires representation to be "
                "'modified_fraction' or 'modified_count'."
            )
        if test not in {"effect_size_only", "beta_binomial"}:
            raise ValueError(
                "Current region_contrasts scoring support requires test in "
                "{'effect_size_only', 'beta_binomial'}."
            )
        return

    if analysis_unit == "cluster_occupancy":
        if signal_source != "cluster_occupancy":
            raise ValueError(
                "V1 region_contrasts inference requires signal_source='cluster_occupancy'."
            )
        if representation not in {
            "cluster_fraction",
            "dominant_cluster",
            "cluster_entropy",
        }:
            raise ValueError(
                "V1 region_contrasts inference requires representation to be "
                "'cluster_fraction', 'dominant_cluster', or 'cluster_entropy'."
            )
        if test not in {"effect_size_only", "fraction_test"}:
            raise ValueError(
                "Current cluster_occupancy scoring support requires test in "
                "{'effect_size_only', 'fraction_test'}."
            )
        return

    if analysis_unit != "ensemble_region":
        raise ValueError(
            "V1 region_contrasts inference requires analysis_unit='ensemble_region'."
        )


def build_region_evidence_table(
    *,
    samples,
    regions,
    motifs: Iterable[str],
    window_size: int | None = None,
    quiet: bool = True,
    cores: int | None = None,
) -> pd.DataFrame:
    motifs = list(motifs)
    if len(motifs) != 1:
        raise ValueError("build_region_evidence_table currently supports exactly one motif.")

    motif = motifs[0]
    regions_dict = utils.regions_dict_from_input(regions, window_size)
    ordered_regions = [
        (chromosome, start, end, strand)
        for chromosome, region_list in regions_dict.items()
        for start, end, strand in region_list
    ]

    rows = []
    for sample in samples:
        metadata = sample.metadata or {}
        if "pileup_path" not in metadata:
            raise ValueError(
                f"Sample {sample.sample_id} is missing metadata['pileup_path']."
            )

        pileup_loader = partial(
            load_processed.pileup_counts_from_bedmethyl,
            bedmethyl_file=metadata["pileup_path"],
            motif=motif,
        )
        counts_by_region = load_processed.regions_to_list(
            function_handle=pileup_loader,
            regions=regions,
            window_size=window_size,
            quiet=quiet,
            cores=cores,
        )

        if len(counts_by_region) != len(ordered_regions):
            raise ValueError(
                "Pileup evidence count length did not match the number of parsed regions."
            )

        for (chromosome, start, end, strand), (
            modified_count,
            valid_count,
        ) in zip(ordered_regions, counts_by_region):
            mod_fraction = 0.0 if valid_count == 0 else modified_count / valid_count
            rows.append(
                {
                    "region_id": f"{chromosome}:{start}-{end},{strand}",
                    "chromosome": chromosome,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "sample_id": sample.sample_id,
                    "condition": sample.condition,
                    "replicate": sample.replicate,
                    "modified_count": modified_count,
                    "valid_count": valid_count,
                    "mod_fraction": mod_fraction,
                }
            )

    return pd.DataFrame(rows)


def build_cluster_occupancy_evidence_table(
    *,
    region_summaries: pd.DataFrame,
) -> pd.DataFrame:
    required_columns = {"region_id", "sample_id", "condition", "cluster", "count"}
    missing_columns = required_columns - set(region_summaries.columns)
    if missing_columns:
        missing_display = ", ".join(sorted(missing_columns))
        raise ValueError(
            "build_cluster_occupancy_evidence_table requires columns: "
            f"{missing_display}."
        )

    evidence = region_summaries.copy()
    grouped_keys = ["region_id", "sample_id", "condition"]
    totals = evidence.groupby(grouped_keys, dropna=False)["count"].transform("sum")
    evidence["fraction"] = (
        evidence["count"].div(totals.where(totals != 0), fill_value=0).fillna(0.0)
    )

    summary = (
        evidence.sort_values(
            by=grouped_keys + ["fraction", "cluster"],
            ascending=[True, True, True, False, True],
            kind="mergesort",
        )
        .drop_duplicates(subset=grouped_keys, keep="first")
        .loc[:, grouped_keys + ["cluster"]]
        .rename(columns={"cluster": "dominant_cluster"})
    )

    def _cluster_entropy(values: pd.Series) -> float:
        probabilities = values.astype(float)
        if probabilities.empty:
            return 0.0
        probabilities = probabilities[probabilities > 0]
        if probabilities.empty:
            return 0.0
        return float(-(probabilities * probabilities.map(math.log2)).sum())

    entropy = (
        evidence.groupby(grouped_keys, dropna=False)["fraction"]
        .apply(_cluster_entropy)
        .reset_index(name="cluster_entropy")
    )

    return evidence.merge(summary, on=grouped_keys, how="left").merge(
        entropy,
        on=grouped_keys,
        how="left",
    )


def _zero_safe_divide(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    return numerator.div(denominator.where(denominator != 0), fill_value=0).fillna(0.0)


def _estimate_beta_binomial_prior(
    denominator_modified_count: int,
    denominator_valid_count: int,
) -> tuple[float, float]:
    if denominator_modified_count < 0:
        raise ValueError("beta-binomial denominator modified_count must be >= 0.")
    if denominator_valid_count < 0:
        raise ValueError("beta-binomial denominator valid_count must be >= 0.")
    if denominator_modified_count > denominator_valid_count:
        raise ValueError(
            "beta-binomial denominator modified_count cannot exceed valid_count."
        )

    denominator_unmodified_count = denominator_valid_count - denominator_modified_count
    return float(denominator_modified_count + 1), float(denominator_unmodified_count + 1)


def _log_beta_function(alpha: float, beta: float) -> float:
    return math.lgamma(alpha) + math.lgamma(beta) - math.lgamma(alpha + beta)


def _beta_binomial_logpmf(k: int, n: int, alpha: float, beta: float) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return (
        math.lgamma(n + 1)
        - math.lgamma(k + 1)
        - math.lgamma(n - k + 1)
        + _log_beta_function(k + alpha, n - k + beta)
        - _log_beta_function(alpha, beta)
    )


def _beta_binomial_two_sided_p_value(k: int, n: int, alpha: float, beta: float) -> float:
    if k < 0:
        raise ValueError("beta-binomial modified_count must be >= 0.")
    if n < 0:
        raise ValueError("beta-binomial valid_count must be >= 0.")
    if k > n:
        raise ValueError("beta-binomial modified_count cannot exceed valid_count.")
    if n <= 0:
        return 1.0

    support = list(range(n + 1))
    logpmf = [_beta_binomial_logpmf(x, n, alpha, beta) for x in support]
    observed_logpmf = logpmf[k]
    tail_probability = sum(
        math.exp(value)
        for value in logpmf
        if value <= observed_logpmf + 1e-12
    )
    return min(max(tail_probability, 0.0), 1.0)


def _adjust_p_values_bh(p_values: pd.Series) -> pd.Series:
    if p_values.empty:
        return pd.Series(dtype=float, index=p_values.index)

    ranked = sorted(enumerate(p_values.tolist()), key=lambda item: item[1], reverse=True)
    total = len(ranked)
    adjusted = [1.0] * total
    running_min = 1.0

    for rank_from_end, (original_index, p_value) in enumerate(ranked, start=1):
        rank = total - rank_from_end + 1
        candidate = min(1.0, p_value * total / rank)
        running_min = min(running_min, candidate)
        adjusted[original_index] = running_min

    return pd.Series(adjusted, index=p_values.index, dtype=float)


def _add_beta_binomial_scores(
    regions_table: pd.DataFrame,
    *,
    multiple_testing: str,
) -> pd.DataFrame:
    if multiple_testing != "fdr_bh":
        raise ValueError(
            "Current beta-binomial region contrast scoring requires "
            "multiple_testing='fdr_bh'."
        )

    scored = regions_table.copy()
    scored["p_value"] = [
        _beta_binomial_two_sided_p_value(
            int(modified_count),
            int(valid_count),
            *_estimate_beta_binomial_prior(
                int(reference_modified_count),
                int(reference_valid_count),
            ),
        )
        for modified_count, valid_count, reference_modified_count, reference_valid_count in zip(
            scored["numerator_modified_count"],
            scored["numerator_valid_count"],
            scored["denominator_modified_count"],
            scored["denominator_valid_count"],
        )
    ]
    scored["adjusted_p_value"] = _adjust_p_values_bh(scored["p_value"])
    return scored


def _pool_region_groups(evidence: pd.DataFrame, contrast: ContrastSpec) -> pd.DataFrame:
    if contrast.mode not in {"pairwise", "group_vs_group"}:
        raise NotImplementedError(
            f"Effect-size pooling is not implemented for contrast mode '{contrast.mode}'."
        )

    pooled_frames = []
    side_specs = {
        "numerator": contrast.numerator or [],
        "denominator": contrast.denominator or [],
    }

    for side, conditions in side_specs.items():
        side_evidence = evidence.loc[evidence["condition"].isin(conditions)].copy()
        available_conditions = set(side_evidence["condition"].dropna().unique())
        missing_conditions = sorted(set(conditions) - available_conditions)
        if missing_conditions:
            missing_display = ", ".join(missing_conditions)
            raise ValueError(
                f"Missing {side} evidence for requested condition(s): {missing_display}."
            )
        pooled = (
            side_evidence.groupby(
                ["region_id", "chromosome", "start", "end", "strand"],
                dropna=False,
                sort=False,
                as_index=False,
            )
            .agg(
                modified_count=("modified_count", "sum"),
                valid_count=("valid_count", "sum"),
                replicate_n=("sample_id", "nunique"),
            )
            .assign(contrast_side=side)
        )
        pooled["fraction"] = _zero_safe_divide(
            pooled["modified_count"], pooled["valid_count"]
        )
        pooled_frames.append(pooled)

    return pd.concat(pooled_frames, ignore_index=True)


def score_regions(
    *,
    samples,
    regions,
    motifs,
    contrast: ContrastSpec,
    analysis_unit: str = "ensemble_region",
    representation: str = "modified_fraction",
    signal_source: str = "pileup_counts",
    test: str = "effect_size_only",
    multiple_testing: str = "fdr_bh",
) -> RegionContrastResult:
    validate_region_contrast_request(
        analysis_unit=analysis_unit,
        representation=representation,
        signal_source=signal_source,
        test=test,
    )
    evidence = build_region_evidence_table(
        samples=samples,
        regions=regions,
        motifs=motifs,
    )

    pooled = _pool_region_groups(evidence=evidence, contrast=contrast)
    numerator = (
        pooled.loc[pooled["contrast_side"] == "numerator"]
        .drop(columns=["contrast_side"])
        .rename(
            columns={
                "modified_count": "numerator_modified_count",
                "valid_count": "numerator_valid_count",
                "replicate_n": "numerator_replicate_n",
                "fraction": "fraction",
            }
        )
    )
    denominator = (
        pooled.loc[pooled["contrast_side"] == "denominator"]
        .drop(columns=["contrast_side"])
        .rename(
            columns={
                "modified_count": "denominator_modified_count",
                "valid_count": "denominator_valid_count",
                "replicate_n": "denominator_replicate_n",
                "fraction": "reference_fraction",
            }
        )
    )

    merged = numerator.merge(
        denominator,
        on=["region_id", "chromosome", "start", "end", "strand"],
        how="outer",
        sort=False,
    )

    for column in [
        "numerator_modified_count",
        "numerator_valid_count",
        "numerator_replicate_n",
        "denominator_modified_count",
        "denominator_valid_count",
        "denominator_replicate_n",
    ]:
        merged[column] = merged[column].fillna(0)

    merged["fraction"] = merged["fraction"].fillna(0.0)
    merged["reference_fraction"] = merged["reference_fraction"].fillna(0.0)
    merged["delta_fraction"] = merged["fraction"] - merged["reference_fraction"]
    pseudocount = 1e-6
    merged["log2_fc"] = (
        (merged["fraction"] + pseudocount) / (merged["reference_fraction"] + pseudocount)
    ).map(math.log2)

    integer_columns = [
        "numerator_modified_count",
        "numerator_valid_count",
        "numerator_replicate_n",
        "denominator_modified_count",
        "denominator_valid_count",
        "denominator_replicate_n",
    ]
    merged[integer_columns] = merged[integer_columns].astype(int)

    if representation == "modified_count":
        merged["count"] = merged["numerator_modified_count"]
        merged["reference_count"] = merged["denominator_modified_count"]
        merged["delta_count"] = merged["count"] - merged["reference_count"]
        merged["log2_fc_count"] = (
            (merged["count"] + pseudocount) / (merged["reference_count"] + pseudocount)
        ).map(math.log2)
        merged["effect_size"] = merged["delta_count"].abs()
    else:
        merged["effect_size"] = merged["delta_fraction"].abs()

    if test == "beta_binomial":
        regions_table = _add_beta_binomial_scores(
            merged,
            multiple_testing=multiple_testing,
        )
        regions_table = regions_table.sort_values(
            by=["adjusted_p_value", "p_value"],
            ascending=[True, True],
            kind="mergesort",
        ).reset_index(drop=True)
    else:
        regions_table = merged.sort_values(
            by="effect_size",
            ascending=False,
            kind="mergesort",
        ).reset_index(drop=True)
    regions_table["rank"] = range(1, len(regions_table) + 1)

    summary_columns = [
        "region_id",
        "fraction",
        "reference_fraction",
        "delta_fraction",
        "log2_fc",
        "rank",
        "numerator_modified_count",
        "numerator_valid_count",
        "numerator_replicate_n",
        "denominator_modified_count",
        "denominator_valid_count",
        "denominator_replicate_n",
    ]
    if test == "beta_binomial":
        summary_columns.extend(["p_value", "adjusted_p_value"])
    if representation == "modified_count":
        summary_columns.extend(["count", "reference_count", "delta_count", "log2_fc_count"])
    summary = regions_table.loc[:, summary_columns].copy()

    contrast_metadata = contrast.metadata or {}
    metadata = {
        "contrast_mode": contrast.mode,
        "analysis_unit": analysis_unit,
        "representation": representation,
        "signal_source": signal_source,
        "test": test,
        "multiple_testing": multiple_testing,
        "normalization_mode": contrast_metadata.get("normalization_mode", "none"),
        "biological_interpretation": contrast_metadata.get(
            "biological_interpretation",
            "region-level difference in pooled modified fraction",
        ),
        "renderer": contrast_metadata.get("renderer", "region_effect_sizes"),
    }

    return RegionContrastResult(
        regions=regions_table,
        summary=summary,
        contrast=contrast,
        plot_data={"region_effect_sizes": summary.copy()},
        metadata=metadata,
        figures={},
    )
