from __future__ import annotations

from functools import partial
import math
from typing import Iterable

import pandas as pd
from scipy import stats

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
        if representation == "cluster_fraction" and test not in {
            "effect_size_only",
            "fraction_test",
        }:
            raise ValueError(
                "Current cluster_occupancy cluster_fraction scoring support requires "
                "test in "
                "{'effect_size_only', 'fraction_test'}."
            )
        if representation in {"dominant_cluster", "cluster_entropy"} and test != (
            "effect_size_only"
        ):
            raise ValueError(
                "Current cluster_occupancy descriptive scoring requires "
                "representation='cluster_fraction' for inferential tests."
            )
        return

    if analysis_unit not in {"ensemble_region", "cluster_occupancy"}:
        raise ValueError(
            "V1 region_contrasts inference requires analysis_unit='ensemble_region' "
            "or analysis_unit='cluster_occupancy'."
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


def _require_occupancy_columns(
    occupancy_table: pd.DataFrame,
    *,
    required_columns: set[str],
) -> None:
    missing_columns = required_columns - set(occupancy_table.columns)
    if missing_columns:
        missing_display = ", ".join(sorted(missing_columns))
        raise ValueError(f"occupancy_table requires columns: {missing_display}.")


def _pool_cluster_occupancy_groups(
    evidence: pd.DataFrame,
    contrast: ContrastSpec,
    *,
    value_column: str,
    group_columns: list[str],
) -> pd.DataFrame:
    if contrast.mode not in {"pairwise", "group_vs_group"}:
        raise NotImplementedError(
            f"Cluster occupancy scoring is not implemented for contrast mode '{contrast.mode}'."
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
            side_evidence.groupby(group_columns, dropna=False, sort=False)
            .agg(
                value=(value_column, "mean"),
                replicate_n=("sample_id", "nunique"),
                sample_values=(value_column, lambda values: tuple(float(v) for v in values)),
            )
            .reset_index()
            .assign(contrast_side=side)
        )
        pooled_frames.append(pooled)

    return pd.concat(pooled_frames, ignore_index=True)


def _add_fraction_test_scores(
    regions_table: pd.DataFrame,
    *,
    multiple_testing: str,
) -> pd.DataFrame:
    if multiple_testing != "fdr_bh":
        raise ValueError(
            "Current fraction_test region contrast scoring requires "
            "multiple_testing='fdr_bh'."
        )

    def _welch_p_value(row: pd.Series) -> float:
        numerator_values = list(row["numerator_sample_values"])
        denominator_values = list(row["denominator_sample_values"])
        if len(numerator_values) < 2 or len(denominator_values) < 2:
            return 1.0

        test_result = stats.ttest_ind(
            numerator_values,
            denominator_values,
            equal_var=False,
        )
        p_value = float(test_result.pvalue)
        if math.isnan(p_value):
            return 1.0
        return min(max(p_value, 0.0), 1.0)

    scored = regions_table.copy()
    scored["p_value"] = scored.apply(_welch_p_value, axis=1)
    scored["adjusted_p_value"] = _adjust_p_values_bh(scored["p_value"])
    return scored


def _dominant_label(values: pd.Series) -> str | None:
    non_null = values.dropna()
    if non_null.empty:
        return None
    counts = non_null.value_counts()
    max_count = counts.max()
    return sorted(counts[counts == max_count].index.tolist())[0]


def _score_cluster_occupancy(
    *,
    evidence: pd.DataFrame,
    contrast: ContrastSpec,
    representation: str,
    test: str,
    multiple_testing: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pseudocount = 1e-6

    if representation == "cluster_fraction":
        _require_occupancy_columns(
            evidence,
            required_columns={"region_id", "sample_id", "condition", "cluster", "fraction"},
        )
        pooled = _pool_cluster_occupancy_groups(
            evidence,
            contrast,
            value_column="fraction",
            group_columns=["region_id", "cluster"],
        )
        numerator = (
            pooled.loc[pooled["contrast_side"] == "numerator"]
            .drop(columns=["contrast_side"])
            .rename(
                columns={
                    "value": "fraction",
                    "replicate_n": "numerator_replicate_n",
                    "sample_values": "numerator_sample_values",
                }
            )
        )
        denominator = (
            pooled.loc[pooled["contrast_side"] == "denominator"]
            .drop(columns=["contrast_side"])
            .rename(
                columns={
                    "value": "reference_fraction",
                    "replicate_n": "denominator_replicate_n",
                    "sample_values": "denominator_sample_values",
                }
            )
        )
        regions_table = numerator.merge(
            denominator,
            on=["region_id", "cluster"],
            how="outer",
            sort=False,
        )
        regions_table["fraction"] = regions_table["fraction"].fillna(0.0)
        regions_table["reference_fraction"] = regions_table["reference_fraction"].fillna(0.0)
        regions_table["numerator_replicate_n"] = (
            regions_table["numerator_replicate_n"].fillna(0).astype(int)
        )
        regions_table["denominator_replicate_n"] = (
            regions_table["denominator_replicate_n"].fillna(0).astype(int)
        )
        regions_table["numerator_sample_values"] = regions_table[
            "numerator_sample_values"
        ].apply(lambda values: values if isinstance(values, tuple) else tuple())
        regions_table["denominator_sample_values"] = regions_table[
            "denominator_sample_values"
        ].apply(lambda values: values if isinstance(values, tuple) else tuple())
        regions_table["delta_fraction"] = (
            regions_table["fraction"] - regions_table["reference_fraction"]
        )
        regions_table["log2_fc"] = (
            (regions_table["fraction"] + pseudocount)
            / (regions_table["reference_fraction"] + pseudocount)
        ).map(math.log2)
        regions_table["effect_size"] = regions_table["delta_fraction"].abs()

        if test == "fraction_test":
            regions_table = _add_fraction_test_scores(
                regions_table,
                multiple_testing=multiple_testing,
            )
            regions_table = regions_table.sort_values(
                by=["adjusted_p_value", "p_value", "effect_size", "region_id", "cluster"],
                ascending=[True, True, False, True, True],
                kind="mergesort",
            ).reset_index(drop=True)
        else:
            regions_table = regions_table.sort_values(
                by=["effect_size", "delta_fraction", "region_id", "cluster"],
                ascending=[False, False, True, True],
                kind="mergesort",
            ).reset_index(drop=True)
        regions_table["rank"] = range(1, len(regions_table) + 1)

        summary_columns = [
            "region_id",
            "cluster",
            "fraction",
            "reference_fraction",
            "delta_fraction",
            "log2_fc",
            "rank",
            "numerator_replicate_n",
            "denominator_replicate_n",
        ]
        if test == "fraction_test":
            summary_columns.extend(["p_value", "adjusted_p_value"])
        summary = regions_table.loc[:, summary_columns].copy()
        return regions_table, summary

    if representation == "dominant_cluster":
        _require_occupancy_columns(
            evidence,
            required_columns={"region_id", "sample_id", "condition", "dominant_cluster"},
        )
        pooled = (
            evidence.loc[:, ["region_id", "sample_id", "condition", "dominant_cluster"]]
            .drop_duplicates()
            .groupby(["region_id", "condition"], dropna=False, sort=False)
            .agg(
                dominant_cluster=("dominant_cluster", _dominant_label),
                replicate_n=("sample_id", "nunique"),
            )
            .reset_index()
        )
        numerator = (
            pooled.loc[pooled["condition"].isin(contrast.numerator or [])]
            .groupby("region_id", dropna=False, sort=False)
            .agg(
                dominant_cluster=("dominant_cluster", _dominant_label),
                numerator_replicate_n=("replicate_n", "sum"),
            )
            .reset_index()
        )
        denominator = (
            pooled.loc[pooled["condition"].isin(contrast.denominator or [])]
            .groupby("region_id", dropna=False, sort=False)
            .agg(
                reference_dominant_cluster=("dominant_cluster", _dominant_label),
                denominator_replicate_n=("replicate_n", "sum"),
            )
            .reset_index()
        )
        regions_table = numerator.merge(
            denominator,
            on="region_id",
            how="outer",
            sort=False,
        )
        regions_table["dominant_cluster_changed"] = (
            regions_table["dominant_cluster"] != regions_table["reference_dominant_cluster"]
        )
        regions_table["effect_size"] = regions_table["dominant_cluster_changed"].astype(int)
        regions_table = regions_table.sort_values(
            by=["effect_size", "region_id"],
            ascending=[False, True],
            kind="mergesort",
        ).reset_index(drop=True)
        regions_table["rank"] = range(1, len(regions_table) + 1)
        summary = regions_table.loc[
            :,
            [
                "region_id",
                "dominant_cluster",
                "reference_dominant_cluster",
                "dominant_cluster_changed",
                "rank",
                "numerator_replicate_n",
                "denominator_replicate_n",
            ],
        ].copy()
        return regions_table, summary

    _require_occupancy_columns(
        evidence,
        required_columns={"region_id", "sample_id", "condition", "cluster_entropy"},
    )
    pooled = _pool_cluster_occupancy_groups(
        evidence.loc[:, ["region_id", "sample_id", "condition", "cluster_entropy"]]
        .drop_duplicates(),
        contrast,
        value_column="cluster_entropy",
        group_columns=["region_id"],
    )
    numerator = (
        pooled.loc[pooled["contrast_side"] == "numerator"]
        .drop(columns=["contrast_side", "sample_values"])
        .rename(
            columns={
                "value": "cluster_entropy",
                "replicate_n": "numerator_replicate_n",
            }
        )
    )
    denominator = (
        pooled.loc[pooled["contrast_side"] == "denominator"]
        .drop(columns=["contrast_side", "sample_values"])
        .rename(
            columns={
                "value": "reference_cluster_entropy",
                "replicate_n": "denominator_replicate_n",
            }
        )
    )
    regions_table = numerator.merge(
        denominator,
        on="region_id",
        how="outer",
        sort=False,
    )
    regions_table["cluster_entropy"] = regions_table["cluster_entropy"].fillna(0.0)
    regions_table["reference_cluster_entropy"] = regions_table[
        "reference_cluster_entropy"
    ].fillna(0.0)
    regions_table["numerator_replicate_n"] = (
        regions_table["numerator_replicate_n"].fillna(0).astype(int)
    )
    regions_table["denominator_replicate_n"] = (
        regions_table["denominator_replicate_n"].fillna(0).astype(int)
    )
    regions_table["delta_cluster_entropy"] = (
        regions_table["cluster_entropy"] - regions_table["reference_cluster_entropy"]
    )
    regions_table["effect_size"] = regions_table["delta_cluster_entropy"].abs()
    regions_table = regions_table.sort_values(
        by=["effect_size", "region_id"],
        ascending=[False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    regions_table["rank"] = range(1, len(regions_table) + 1)
    summary = regions_table.loc[
        :,
        [
            "region_id",
            "cluster_entropy",
            "reference_cluster_entropy",
            "delta_cluster_entropy",
            "rank",
            "numerator_replicate_n",
            "denominator_replicate_n",
        ],
    ].copy()
    return regions_table, summary


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
    occupancy_table: pd.DataFrame | None = None,
) -> RegionContrastResult:
    validate_region_contrast_request(
        analysis_unit=analysis_unit,
        representation=representation,
        signal_source=signal_source,
        test=test,
    )
    contrast_metadata = contrast.metadata or {}

    if analysis_unit == "cluster_occupancy":
        if occupancy_table is None:
            raise ValueError(
                "cluster_occupancy scoring currently requires occupancy_table."
            )
        regions_table, summary = _score_cluster_occupancy(
            evidence=occupancy_table.copy(),
            contrast=contrast,
            representation=representation,
            test=test,
            multiple_testing=multiple_testing,
        )
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
                "region-level difference in sample-level cluster occupancy",
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
