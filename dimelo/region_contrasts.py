from __future__ import annotations

import math
from collections.abc import Iterable
from functools import partial

import pandas as pd
from scipy import stats

from . import cluster, load_processed, utils
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

    if analysis_unit == "single_read":
        if representation == "read_mod_fraction":
            if signal_source != "extract_reads":
                raise ValueError(
                    "Current single_read read_mod_fraction support requires "
                    "signal_source='extract_reads'."
                )
            if test not in {"effect_size_only", "sample_distribution_shift"}:
                raise ValueError(
                    "Current single_read read_mod_fraction scoring support requires "
                    "test in {'effect_size_only', 'sample_distribution_shift'}."
                )
            return

        if representation == "read_window_features":
            if signal_source != "extract_features":
                raise ValueError(
                    "Current single_read read_window_features support requires "
                    "signal_source='extract_features'."
                )
            if test not in {"effect_size_only", "feature_summary_shift"}:
                raise ValueError(
                    "Current single_read read_window_features support requires "
                    "test in {'effect_size_only', 'feature_summary_shift'}."
                )
            return

        raise ValueError(
            "Current single_read support requires representation='read_mod_fraction' "
            "or 'read_window_features'."
        )

    if analysis_unit not in {"ensemble_region", "cluster_occupancy", "single_read"}:
        raise ValueError(
            "V1 region_contrasts inference requires analysis_unit='ensemble_region' "
            "or analysis_unit='cluster_occupancy' "
            "or analysis_unit='single_read'."
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
        raise ValueError(
            "build_region_evidence_table currently supports exactly one motif."
        )

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
        ) in zip(ordered_regions, counts_by_region, strict=False):
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


def build_single_read_mod_fraction_evidence_table(
    *,
    extract_table: pd.DataFrame,
    pairing_key: str | None = None,
) -> pd.DataFrame:
    required_columns = {
        "region_id",
        "sample_id",
        "condition",
        "read_id",
        "modified_count",
        "valid_count",
    }
    if pairing_key:
        required_columns.add(pairing_key)
    missing_columns = required_columns - set(extract_table.columns)
    if missing_columns:
        missing_display = ", ".join(sorted(missing_columns))
        raise ValueError(
            "build_single_read_mod_fraction_evidence_table requires columns: "
            f"{missing_display}."
        )

    selected_columns = [
        "region_id",
        "sample_id",
        "condition",
        "read_id",
        "modified_count",
        "valid_count",
    ]
    if pairing_key:
        selected_columns.append(pairing_key)
    evidence = extract_table.loc[
        :,
        selected_columns,
    ].copy()
    evidence["modified_count"] = pd.to_numeric(
        evidence["modified_count"], errors="coerce"
    )
    evidence["valid_count"] = pd.to_numeric(evidence["valid_count"], errors="coerce")
    invalid_counts = (
        evidence["modified_count"].isna()
        | evidence["valid_count"].isna()
        | ~evidence["modified_count"].map(math.isfinite)
        | ~evidence["valid_count"].map(math.isfinite)
        | (evidence["modified_count"] % 1 != 0)
        | (evidence["valid_count"] % 1 != 0)
        | (evidence["modified_count"] < 0)
        | (evidence["valid_count"] < 0)
        | (evidence["modified_count"] > evidence["valid_count"])
    )
    if invalid_counts.any():
        raise ValueError(
            "build_single_read_mod_fraction_evidence_table modified_count and "
            "valid_count must be finite integers, >= 0, and modified_count <= valid_count."
        )
    evidence["modified_count"] = evidence["modified_count"].astype(int)
    evidence["valid_count"] = evidence["valid_count"].astype(int)
    evidence["read_mod_fraction"] = (
        evidence["modified_count"]
        .div(
            evidence["valid_count"].where(evidence["valid_count"] != 0),
            fill_value=0,
        )
        .fillna(0.0)
    )
    return evidence


def build_single_read_feature_evidence_table(
    *,
    feature_table: pd.DataFrame,
    pairing_key: str | None = None,
) -> pd.DataFrame:
    required_columns = {"region_id", "sample_id", "condition", "read_id"}
    if pairing_key:
        required_columns.add(pairing_key)
    missing_columns = required_columns - set(feature_table.columns)
    if missing_columns:
        missing_display = ", ".join(sorted(missing_columns))
        raise ValueError(
            "build_single_read_feature_evidence_table requires columns: "
            f"{missing_display}."
        )

    feature_columns = [
        column for column in feature_table.columns if column not in required_columns
    ]
    if not feature_columns:
        raise ValueError(
            "build_single_read_feature_evidence_table requires at least one feature column."
        )
    selected_columns = ["region_id", "sample_id", "condition", "read_id"]
    if pairing_key:
        selected_columns.append(pairing_key)
    selected_columns.extend(feature_columns)
    evidence = feature_table.loc[:, selected_columns].copy()
    numeric_features = evidence.loc[:, feature_columns].apply(
        pd.to_numeric,
        errors="coerce",
    )
    finite_feature_values = numeric_features.apply(
        lambda column: column.map(math.isfinite)
    )
    invalid_feature_values = numeric_features.isna() | ~finite_feature_values
    if invalid_feature_values.to_numpy().any():
        raise ValueError(
            "build_single_read_feature_evidence_table feature columns must be numeric "
            "and finite."
        )
    evidence.loc[:, feature_columns] = numeric_features
    return evidence


def _require_single_read_pairing_key(frame: pd.DataFrame, pairing_key: str) -> None:
    if pairing_key not in frame.columns:
        raise ValueError(
            f"single_read matched_pairwise scoring requires column {pairing_key!r}."
        )
    if frame[pairing_key].isna().any():
        raise ValueError(
            f"single_read matched_pairwise scoring requires non-null values in column {pairing_key!r}."
        )


def _require_single_read_matched_binary_sides(contrast: ContrastSpec) -> None:
    if len(contrast.numerator or []) != 1 or len(contrast.denominator or []) != 1:
        raise ValueError(
            "single_read matched_pairwise scoring currently requires exactly one "
            "numerator and one denominator condition."
        )


def _validate_single_read_matched_sample_summary(
    sample_summary: pd.DataFrame,
    *,
    pairing_key: str,
) -> None:
    sample_pair_counts = (
        sample_summary.loc[:, ["sample_id", pairing_key]]
        .drop_duplicates()
        .groupby("sample_id", dropna=False)[pairing_key]
        .nunique(dropna=False)
    )
    if (sample_pair_counts > 1).any():
        raise ValueError(
            "single_read matched_pairwise scoring requires each sample_id to map "
            "to exactly one pairing key."
        )

    side_counts = sample_summary.groupby(
        ["region_id", pairing_key, "condition"],
        dropna=False,
        sort=False,
    )["sample_id"].nunique()
    if (side_counts > 1).any():
        raise ValueError(
            "single_read matched_pairwise scoring requires exactly one sample "
            "per region, pair, and condition."
        )


def _load_builtin_single_read_feature_table(*, samples, regions, motifs):
    if not samples:
        raise ValueError(
            "Built-in single_read feature loading requires at least one sample."
        )

    selected_feature_names: list[str] | None = None
    sample_frames: list[pd.DataFrame] = []
    for sample in samples:
        extracted = cluster.extract_read_windows(
            hdf5_file=sample.extract_h5,
            motifs=motifs,
            regions=regions,
        )
        feature_matrix, feature_names = cluster.read_window_feature_matrix(extracted)
        feature_rows = list(feature_matrix)
        feature_names = list(feature_names)
        if selected_feature_names is None:
            selected_feature_names = feature_names
        elif selected_feature_names != feature_names:
            raise ValueError(
                "Built-in single_read feature loading requires matching feature names "
                "across samples."
            )

        if len(extracted.metadata) != len(feature_rows):
            raise ValueError(
                "Built-in single_read feature loading requires one metadata row per "
                "feature row."
            )

        rows = []
        for metadata, feature_values in zip(
            extracted.metadata, feature_rows, strict=False
        ):
            chromosome = metadata.get("chromosome")
            start = metadata.get("region_start")
            end = metadata.get("region_end")
            read_id = metadata.get("read_name")
            if chromosome is None or start is None or end is None:
                raise ValueError(
                    "Built-in single_read feature loading requires metadata with "
                    "chromosome, region_start, and region_end."
                )
            if read_id is None:
                raise ValueError(
                    "Built-in single_read feature loading requires metadata with read_name."
                )
            strand = metadata.get("region_strand") or "+"
            row = {
                "region_id": f"{chromosome}:{int(start)}-{int(end)},{strand}",
                "sample_id": sample.sample_id,
                "condition": sample.condition,
                "read_id": read_id,
            }
            row.update(
                {
                    feature_name: float(feature_value)
                    for feature_name, feature_value in zip(
                        feature_names, feature_values, strict=False
                    )
                }
            )
            rows.append(row)
        sample_frames.append(pd.DataFrame(rows))

    if not sample_frames:
        return pd.DataFrame(
            columns=[
                "region_id",
                "sample_id",
                "condition",
                "read_id",
                *(selected_feature_names or []),
            ]
        )
    return pd.concat(sample_frames, ignore_index=True)


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
    counts = pd.to_numeric(evidence["count"], errors="raise")
    invalid_counts = counts.isna() | (counts < 0) | ~counts.map(math.isfinite)
    if invalid_counts.any():
        raise ValueError(
            "build_cluster_occupancy_evidence_table count values must be finite and >= 0."
        )

    grouped_keys = ["region_id", "sample_id", "condition"]
    totals = evidence.groupby(grouped_keys, dropna=False)["count"].transform("sum")
    evidence["_group_total_count"] = totals
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
        .loc[:, grouped_keys + ["cluster", "_group_total_count"]]
        .rename(columns={"cluster": "dominant_cluster"})
    )
    summary.loc[summary["_group_total_count"] == 0, "dominant_cluster"] = None
    summary = summary.drop(columns="_group_total_count")

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

    return (
        evidence.merge(summary, on=grouped_keys, how="left")
        .merge(
            entropy,
            on=grouped_keys,
            how="left",
        )
        .drop(columns="_group_total_count")
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
    return float(denominator_modified_count + 1), float(
        denominator_unmodified_count + 1
    )


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


def _beta_binomial_two_sided_p_value(
    k: int, n: int, alpha: float, beta: float
) -> float:
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
        math.exp(value) for value in logpmf if value <= observed_logpmf + 1e-12
    )
    return min(max(tail_probability, 0.0), 1.0)


# STATISTICS CHANGE (review fix #2): the previous "differential" test conditioned
# a beta-binomial prior on the denominator's own counts and treated the reference
# rate as exact. That was asymmetric (A-vs-B != B-vs-A) and produced false
# positives when the reference had low coverage (its sampling uncertainty was
# ignored). Replace it with a symmetric two-sample test (Fisher's exact on the
# 2x2 table of modified / unmodified counts for numerator and denominator), which
# accounts for both sides' sampling uncertainty and is order-invariant.
def _two_sample_proportion_p_value(
    numerator_modified_count: int,
    numerator_valid_count: int,
    denominator_modified_count: int,
    denominator_valid_count: int,
) -> float:
    """Symmetric two-sample test on (modified, unmodified) counts.

    Uses Fisher's exact test on the 2x2 contingency table
        [[num_mod, num_unmod], [den_mod, den_unmod]].
    Returns 1.0 when either side has zero valid coverage (structurally
    untestable). This is symmetric in numerator/denominator and does not
    treat either rate as known without error.
    """
    for name, value in (
        ("numerator_modified_count", numerator_modified_count),
        ("numerator_valid_count", numerator_valid_count),
        ("denominator_modified_count", denominator_modified_count),
        ("denominator_valid_count", denominator_valid_count),
    ):
        if value < 0:
            raise ValueError(f"two-sample proportion {name} must be >= 0.")
    if numerator_modified_count > numerator_valid_count:
        raise ValueError(
            "two-sample proportion numerator modified_count cannot exceed valid_count."
        )
    if denominator_modified_count > denominator_valid_count:
        raise ValueError(
            "two-sample proportion denominator modified_count cannot exceed valid_count."
        )

    if numerator_valid_count <= 0 or denominator_valid_count <= 0:
        return 1.0

    numerator_unmodified_count = numerator_valid_count - numerator_modified_count
    denominator_unmodified_count = denominator_valid_count - denominator_modified_count
    _odds_ratio, p_value = stats.fisher_exact(
        [
            [numerator_modified_count, numerator_unmodified_count],
            [denominator_modified_count, denominator_unmodified_count],
        ],
        alternative="two-sided",
    )
    p_value = float(p_value)
    if math.isnan(p_value):
        return 1.0
    return min(max(p_value, 0.0), 1.0)


def _two_sample_beta_binomial_lrt(
    numerator_pairs: list[tuple[int, int]],
    denominator_pairs: list[tuple[int, int]],
) -> float:
    """Symmetric two-sample beta-binomial likelihood-ratio test (maintainer choice).

    Each side is a list of per-replicate ``(modified_count, valid_count)`` pairs. Both
    sides share an overdispersion parameter ``rho`` (intra-class correlation); the null
    is equal methylation mean across sides, the alternative allows different means. The
    statistic ``2*(loglik_full - loglik_null)`` is referred to chi-square with 1 df.

    Unlike a pooled binomial / Fisher's-exact test, this accounts for BOTH sides'
    sampling uncertainty AND replicate-to-replicate overdispersion, and is symmetric in
    numerator/denominator. Returns a two-sided p-value in [0, 1]; returns 1.0 when either
    side has no valid coverage (structurally untestable). With a single replicate per
    side ``rho`` is weakly identified and the test is conservative.
    """
    from scipy.optimize import minimize

    num = [(int(k), int(n)) for k, n in numerator_pairs if int(n) > 0]
    den = [(int(k), int(n)) for k, n in denominator_pairs if int(n) > 0]
    if not num or not den:
        return 1.0

    def _sig(x: float) -> float:  # logistic, guarded off 0/1
        if x >= 0:
            z = math.exp(-x)
            p = 1.0 / (1.0 + z)
        else:
            z = math.exp(x)
            p = z / (1.0 + z)
        return min(max(p, 1e-9), 1.0 - 1e-9)

    def _logit(p: float) -> float:
        p = min(max(p, 1e-6), 1.0 - 1e-6)
        return math.log(p / (1.0 - p))

    def _bb_ll(pairs: list[tuple[int, int]], mu: float, rho: float) -> float:
        # (mu, rho) -> (alpha, beta); rho -> 0 recovers the binomial (no overdispersion)
        alpha = mu * (1.0 - rho) / rho
        beta = (1.0 - mu) * (1.0 - rho) / rho
        return math.fsum(_beta_binomial_logpmf(k, n, alpha, beta) for k, n in pairs)

    def _nll_full(theta: np.ndarray) -> float:
        return -(
            _bb_ll(num, _sig(theta[0]), _sig(theta[2]))
            + _bb_ll(den, _sig(theta[1]), _sig(theta[2]))
        )

    def _nll_null(theta: np.ndarray) -> float:
        mu, rho = _sig(theta[0]), _sig(theta[1])
        return -(_bb_ll(num, mu, rho) + _bb_ll(den, mu, rho))

    num_k = sum(k for k, _ in num); num_n = sum(n for _, n in num)
    den_k = sum(k for k, _ in den); den_n = sum(n for _, n in den)
    fn = num_k / num_n
    fd = den_k / den_n
    fp = (num_k + den_k) / (num_n + den_n)
    rho_init = _logit(0.05)

    full = minimize(_nll_full, [_logit(fn), _logit(fd), rho_init], method="Nelder-Mead")
    null = minimize(_nll_null, [_logit(fp), rho_init], method="Nelder-Mead")
    stat = 2.0 * (float(null.fun) - float(full.fun))
    if not math.isfinite(stat) or stat < 0.0:
        stat = 0.0
    return float(stats.chi2.sf(stat, df=1))


def _replicate_counts_by_region(
    evidence: pd.DataFrame, contrast: ContrastSpec
) -> dict[object, tuple[list[tuple[int, int]], list[tuple[int, int]]]]:
    """Per-region, per-side lists of (modified_count, valid_count) across replicates.

    Feeds the two-sample beta-binomial LRT with replicate-level counts (before pooling),
    so overdispersion can be estimated. Sides are defined by the contrast's numerator /
    denominator condition lists.
    """
    num_conditions = set(contrast.numerator or [])
    den_conditions = set(contrast.denominator or [])
    out: dict[object, tuple[list[tuple[int, int]], list[tuple[int, int]]]] = {}
    for region_id, region_rows in evidence.groupby("region_id", sort=False):
        num_pairs: list[tuple[int, int]] = []
        den_pairs: list[tuple[int, int]] = []
        for condition, modified_count, valid_count in zip(
            region_rows["condition"],
            region_rows["modified_count"],
            region_rows["valid_count"],
            strict=False,
        ):
            pair = (int(modified_count), int(valid_count))
            if condition in num_conditions:
                num_pairs.append(pair)
            if condition in den_conditions:
                den_pairs.append(pair)
        out[region_id] = (num_pairs, den_pairs)
    return out


def _adjust_p_values_bh(
    p_values: pd.Series,
    *,
    testable: pd.Series | None = None,
) -> pd.Series:
    """Benjamini-Hochberg adjustment.

    STATISTICS CHANGE (review fix #1): only genuinely testable rows are counted
    in the multiple-testing total ``m``. Structurally-untestable rows
    (fabricated/degenerate, e.g. p forced to 1.0 because a side has no valid
    coverage or too few replicates) are kept in the output with
    ``adjusted_p_value = NaN`` and excluded from ``m`` so they no longer inflate
    ``m`` and deflate everyone else's adjusted p-values.
    """
    if p_values.empty:
        return pd.Series(dtype=float, index=p_values.index)

    # Work positionally so a non-unique index cannot cause label-based ambiguity.
    p_array = p_values.to_numpy(dtype=float)
    if testable is None:
        testable_array = [True] * len(p_array)
    else:
        aligned = testable.reindex(p_values.index).fillna(False)
        testable_array = aligned.astype(bool).to_numpy().tolist()

    adjusted_array = [float("nan")] * len(p_array)

    testable_positions = [
        position
        for position in range(len(p_array))
        # NaN p-values are never testable.
        if testable_array[position] and not math.isnan(p_array[position])
    ]
    total = len(testable_positions)
    if total == 0:
        return pd.Series(adjusted_array, index=p_values.index, dtype=float)

    ranked = sorted(
        testable_positions, key=lambda position: p_array[position], reverse=True
    )
    running_min = 1.0
    for rank_from_end, position in enumerate(ranked, start=1):
        rank = total - rank_from_end + 1
        candidate = min(1.0, p_array[position] * total / rank)
        running_min = min(running_min, candidate)
        adjusted_array[position] = running_min

    return pd.Series(adjusted_array, index=p_values.index, dtype=float)


def _add_beta_binomial_scores(
    regions_table: pd.DataFrame,
    *,
    multiple_testing: str,
    replicate_counts: dict[object, tuple[list[tuple[int, int]], list[tuple[int, int]]]],
) -> pd.DataFrame:
    if multiple_testing != "fdr_bh":
        raise ValueError(
            "Current beta-binomial region contrast scoring requires "
            "multiple_testing='fdr_bh'."
        )

    scored = regions_table.copy()
    # STATISTICS CHANGE (review fix #2, maintainer choice): symmetric two-sample
    # BETA-BINOMIAL likelihood-ratio test over per-replicate counts, replacing both the
    # old asymmetric conditional beta-binomial and the interim Fisher's-exact stopgap.
    # This models replicate overdispersion and both sides' sampling uncertainty.
    scored["p_value"] = [
        _two_sample_beta_binomial_lrt(*replicate_counts.get(region_id, ([], [])))
        for region_id in scored["region_id"]
    ]
    # A row is only testable when both sides carry valid coverage; fabricated
    # zero-coverage rows (see fix #3) are structurally untestable and excluded
    # from the BH total.
    testable = (scored["numerator_valid_count"] > 0) & (
        scored["denominator_valid_count"] > 0
    )
    scored["adjusted_p_value"] = _adjust_p_values_bh(
        scored["p_value"], testable=testable
    )
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


def _require_supported_cluster_occupancy_mode(contrast: ContrastSpec) -> None:
    if contrast.mode not in {"pairwise", "group_vs_group"}:
        raise NotImplementedError(
            f"Cluster occupancy scoring is not implemented for contrast mode '{contrast.mode}'."
        )


def _require_requested_conditions_present(
    evidence: pd.DataFrame,
    *,
    side_specs: dict[str, list[str]],
) -> None:
    for side, conditions in side_specs.items():
        side_evidence = evidence.loc[evidence["condition"].isin(conditions)].copy()
        available_conditions = set(side_evidence["condition"].dropna().unique())
        missing_conditions = sorted(set(conditions) - available_conditions)
        if missing_conditions:
            missing_display = ", ".join(missing_conditions)
            raise ValueError(
                f"Missing {side} evidence for requested condition(s): {missing_display}."
            )


def _validate_occupancy_table_numeric_columns(occupancy_table: pd.DataFrame) -> None:
    numeric_specs = [
        (
            "fraction",
            "occupancy_table fraction values must be finite and between 0 and 1.",
            lambda values: (
                values.isna() | ~values.map(math.isfinite) | (values < 0) | (values > 1)
            ),
        ),
        (
            "cluster_entropy",
            "occupancy_table cluster_entropy values must be finite and >= 0.",
            lambda values: values.isna() | ~values.map(math.isfinite) | (values < 0),
        ),
        (
            "count",
            "occupancy_table count values must be finite and >= 0.",
            lambda values: values.isna() | ~values.map(math.isfinite) | (values < 0),
        ),
    ]

    for column, message, invalid_mask_builder in numeric_specs:
        if column not in occupancy_table.columns:
            continue
        values = pd.to_numeric(occupancy_table[column], errors="coerce")
        invalid_values = invalid_mask_builder(values)
        if invalid_values.any():
            raise ValueError(message)


def _zero_fill_missing_cluster_fraction_rows(
    side_evidence: pd.DataFrame,
    *,
    cluster_universe: pd.DataFrame,
) -> pd.DataFrame:
    sample_columns = ["region_id", "sample_id", "condition"]
    sample_level = side_evidence.loc[:, sample_columns].drop_duplicates()
    if sample_level.empty:
        return side_evidence

    complete_index = sample_level.merge(cluster_universe, on="region_id", how="inner")
    merged = complete_index.merge(
        side_evidence,
        on=sample_columns + ["cluster"],
        how="left",
        sort=False,
    )
    merged["fraction"] = merged["fraction"].fillna(0.0)
    return merged


def _pool_cluster_occupancy_groups(
    evidence: pd.DataFrame,
    contrast: ContrastSpec,
    *,
    value_column: str,
    group_columns: list[str],
) -> pd.DataFrame:
    _require_supported_cluster_occupancy_mode(contrast)

    pooled_frames = []
    side_specs = {
        "numerator": contrast.numerator or [],
        "denominator": contrast.denominator or [],
    }
    _require_requested_conditions_present(
        evidence,
        side_specs=side_specs,
    )
    cluster_universe = None
    if value_column == "fraction":
        cluster_universe = evidence.loc[:, ["region_id", "cluster"]].drop_duplicates()

    for side, conditions in side_specs.items():
        side_evidence = evidence.loc[evidence["condition"].isin(conditions)].copy()
        if value_column == "fraction":
            side_evidence = _zero_fill_missing_cluster_fraction_rows(
                side_evidence,
                cluster_universe=cluster_universe,
            )

        pooled = (
            side_evidence.groupby(group_columns, dropna=False, sort=False)
            .agg(
                value=(value_column, "mean"),
                replicate_n=("sample_id", "nunique"),
                sample_values=(
                    value_column,
                    lambda values: tuple(float(v) for v in values),
                ),
            )
            .reset_index()
            .assign(contrast_side=side)
        )
        pooled_frames.append(pooled)

    return pd.concat(pooled_frames, ignore_index=True)


def _summarize_single_read_mod_fraction_by_sample(
    evidence: pd.DataFrame,
    *,
    pairing_key: str | None = None,
) -> pd.DataFrame:
    group_columns = ["region_id", "sample_id", "condition"]
    if pairing_key:
        _require_single_read_pairing_key(evidence, pairing_key)
        group_columns.append(pairing_key)

    return evidence.groupby(group_columns, as_index=False, sort=False).agg(
        read_n=("read_id", "nunique"),
        sample_summary_mean=("read_mod_fraction", "mean"),
    )


def _summarize_single_read_mod_fraction_by_pair(
    sample_summary: pd.DataFrame,
    *,
    pairing_key: str,
) -> pd.DataFrame:
    pair_columns = ["region_id", pairing_key, "condition"]
    return sample_summary.groupby(pair_columns, as_index=False, sort=False).agg(
        sample_summary_mean=("sample_summary_mean", "mean"),
        read_n=("read_n", "sum"),
        sample_n=("sample_id", "nunique"),
    )


def _select_complete_matched_pairs(
    pair_summary: pd.DataFrame,
    *,
    pairing_key: str,
    contrast: ContrastSpec,
) -> pd.DataFrame:
    required_conditions = list(
        dict.fromkeys((contrast.numerator or []) + (contrast.denominator or []))
    )
    required_condition_set = set(required_conditions)
    pair_conditions = (
        pair_summary.loc[:, ["region_id", pairing_key, "condition"]]
        .drop_duplicates()
        .groupby(["region_id", pairing_key])["condition"]
        .agg(lambda values: set(values))
    )
    complete_pair_index = [
        pair_index
        for pair_index, conditions in pair_conditions.items()
        if required_condition_set.issubset(conditions)
    ]
    if not complete_pair_index:
        raise ValueError(
            "single_read matched_pairwise scoring found no complete matched pairs."
        )

    complete_pairs = pd.MultiIndex.from_tuples(
        complete_pair_index,
        names=["region_id", pairing_key],
    )
    paired_index = pd.MultiIndex.from_frame(
        pair_summary.loc[:, ["region_id", pairing_key]]
    )
    return pair_summary.loc[paired_index.isin(complete_pairs)].copy()


def _score_single_read_mod_fraction(
    *,
    evidence: pd.DataFrame,
    contrast: ContrastSpec,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if contrast.mode not in {"pairwise", "group_vs_group", "matched_pairwise"}:
        raise NotImplementedError(
            "Single-read mod-fraction scoring is not implemented for contrast mode "
            f"'{contrast.mode}'."
        )

    pairing_key = contrast.pairing_key if contrast.mode == "matched_pairwise" else None
    if contrast.mode == "matched_pairwise":
        _require_single_read_matched_binary_sides(contrast)
    sample_summary = _summarize_single_read_mod_fraction_by_sample(
        evidence,
        pairing_key=pairing_key,
    )
    side_specs = {
        "numerator": contrast.numerator or [],
        "denominator": contrast.denominator or [],
    }
    _require_requested_conditions_present(
        sample_summary,
        side_specs=side_specs,
    )

    if contrast.mode == "matched_pairwise":
        _validate_single_read_matched_sample_summary(
            sample_summary,
            pairing_key=pairing_key or "",
        )
        pair_summary = _summarize_single_read_mod_fraction_by_pair(
            sample_summary,
            pairing_key=pairing_key or "",
        )
        pair_summary = _select_complete_matched_pairs(
            pair_summary,
            pairing_key=pairing_key or "",
            contrast=contrast,
        )

        numerator = (
            pair_summary.loc[pair_summary["condition"].isin(side_specs["numerator"])]
            .drop(columns=["condition"])
            .rename(
                columns={
                    "sample_summary_mean": "sample_summary_numerator_mean",
                    "read_n": "numerator_read_n",
                    "sample_n": "numerator_sample_n",
                }
            )
        )
        denominator = (
            pair_summary.loc[pair_summary["condition"].isin(side_specs["denominator"])]
            .drop(columns=["condition"])
            .rename(
                columns={
                    "sample_summary_mean": "sample_summary_denominator_mean",
                    "read_n": "denominator_read_n",
                    "sample_n": "denominator_sample_n",
                }
            )
        )
        regions_table = numerator.merge(
            denominator,
            on=["region_id", pairing_key],
            how="inner",
            sort=False,
        )
        if regions_table.empty:
            raise ValueError(
                "single_read matched_pairwise scoring found no complete matched pairs."
            )
        regions_table["delta_summary_mean"] = (
            regions_table["sample_summary_numerator_mean"]
            - regions_table["sample_summary_denominator_mean"]
        )
        regions_table["log2_fc"] = (
            (regions_table["sample_summary_numerator_mean"] + 1e-6)
            / (regions_table["sample_summary_denominator_mean"] + 1e-6)
        ).map(math.log2)
        regions_table = regions_table.groupby(
            "region_id", as_index=False, sort=False
        ).agg(
            sample_summary_numerator_mean=("sample_summary_numerator_mean", "mean"),
            sample_summary_denominator_mean=("sample_summary_denominator_mean", "mean"),
            delta_summary_mean=("delta_summary_mean", "mean"),
            numerator_read_n=("numerator_read_n", "sum"),
            denominator_read_n=("denominator_read_n", "sum"),
            numerator_replicate_n=("numerator_sample_n", "size"),
            denominator_replicate_n=("denominator_sample_n", "size"),
        )
        regions_table["log2_fc"] = (
            (regions_table["sample_summary_numerator_mean"] + 1e-6)
            / (regions_table["sample_summary_denominator_mean"] + 1e-6)
        ).map(math.log2)
        regions_table["effect_size"] = regions_table["delta_summary_mean"].abs()
        regions_table = regions_table.sort_values(
            by=["effect_size", "delta_summary_mean", "region_id"],
            ascending=[False, False, True],
            kind="mergesort",
        ).reset_index(drop=True)
        regions_table["rank"] = range(1, len(regions_table) + 1)
        summary = regions_table.loc[
            :,
            [
                "region_id",
                "sample_summary_numerator_mean",
                "sample_summary_denominator_mean",
                "delta_summary_mean",
                "log2_fc",
                "rank",
                "numerator_read_n",
                "denominator_read_n",
                "numerator_replicate_n",
                "denominator_replicate_n",
            ],
        ].copy()
        return regions_table, summary

    pooled_frames = []
    for side, conditions in side_specs.items():
        side_evidence = sample_summary.loc[
            sample_summary["condition"].isin(conditions)
        ].copy()
        pooled = (
            side_evidence.groupby("region_id", dropna=False, sort=False)
            .agg(
                sample_summary_mean=("sample_summary_mean", "mean"),
                read_n=("read_n", "sum"),
                replicate_n=("sample_id", "nunique"),
                sample_values=(
                    "sample_summary_mean",
                    lambda values: tuple(float(v) for v in values),
                ),
            )
            .reset_index()
            .assign(contrast_side=side)
        )
        pooled_frames.append(pooled)

    pooled = pd.concat(pooled_frames, ignore_index=True)
    numerator = (
        pooled.loc[pooled["contrast_side"] == "numerator"]
        .drop(columns=["contrast_side", "sample_values"])
        .rename(
            columns={
                "sample_summary_mean": "sample_summary_numerator_mean",
                "read_n": "numerator_read_n",
                "replicate_n": "numerator_replicate_n",
            }
        )
    )
    denominator = (
        pooled.loc[pooled["contrast_side"] == "denominator"]
        .drop(columns=["contrast_side", "sample_values"])
        .rename(
            columns={
                "sample_summary_mean": "sample_summary_denominator_mean",
                "read_n": "denominator_read_n",
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
    regions_table["sample_summary_numerator_mean"] = regions_table[
        "sample_summary_numerator_mean"
    ].fillna(0.0)
    regions_table["sample_summary_denominator_mean"] = regions_table[
        "sample_summary_denominator_mean"
    ].fillna(0.0)
    regions_table["numerator_read_n"] = (
        regions_table["numerator_read_n"].fillna(0).astype(int)
    )
    regions_table["denominator_read_n"] = (
        regions_table["denominator_read_n"].fillna(0).astype(int)
    )
    regions_table["numerator_replicate_n"] = (
        regions_table["numerator_replicate_n"].fillna(0).astype(int)
    )
    regions_table["denominator_replicate_n"] = (
        regions_table["denominator_replicate_n"].fillna(0).astype(int)
    )
    regions_table["delta_summary_mean"] = (
        regions_table["sample_summary_numerator_mean"]
        - regions_table["sample_summary_denominator_mean"]
    )
    pseudocount = 1e-6
    regions_table["log2_fc"] = (
        (regions_table["sample_summary_numerator_mean"] + pseudocount)
        / (regions_table["sample_summary_denominator_mean"] + pseudocount)
    ).map(math.log2)
    regions_table["effect_size"] = regions_table["delta_summary_mean"].abs()
    regions_table = regions_table.sort_values(
        by=["effect_size", "delta_summary_mean", "region_id"],
        ascending=[False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    regions_table["rank"] = range(1, len(regions_table) + 1)

    summary = regions_table.loc[
        :,
        [
            "region_id",
            "sample_summary_numerator_mean",
            "sample_summary_denominator_mean",
            "delta_summary_mean",
            "log2_fc",
            "rank",
            "numerator_read_n",
            "denominator_read_n",
            "numerator_replicate_n",
            "denominator_replicate_n",
        ],
    ].copy()
    return regions_table, summary


def _summarize_single_read_features_by_sample(
    evidence: pd.DataFrame,
    *,
    pairing_key: str | None = None,
) -> pd.DataFrame:
    group_columns = ["region_id", "sample_id", "condition"]
    excluded_columns = {"region_id", "sample_id", "condition", "read_id"}
    if pairing_key:
        _require_single_read_pairing_key(evidence, pairing_key)
        group_columns.append(pairing_key)
        excluded_columns.add(pairing_key)
    feature_columns = [
        column for column in evidence.columns if column not in excluded_columns
    ]
    return evidence.groupby(group_columns, as_index=False, sort=False)[
        feature_columns
    ].mean()


def _score_single_read_features(
    *,
    evidence: pd.DataFrame,
    contrast: ContrastSpec,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if contrast.mode not in {"pairwise", "group_vs_group", "matched_pairwise"}:
        raise NotImplementedError(
            "Single-read feature scoring is not implemented for contrast mode "
            f"'{contrast.mode}'."
        )

    pairing_key = contrast.pairing_key if contrast.mode == "matched_pairwise" else None
    if contrast.mode == "matched_pairwise":
        _require_single_read_matched_binary_sides(contrast)
    sample_summary = _summarize_single_read_features_by_sample(
        evidence,
        pairing_key=pairing_key,
    )
    side_specs = {
        "numerator": contrast.numerator or [],
        "denominator": contrast.denominator or [],
    }
    _require_requested_conditions_present(
        sample_summary,
        side_specs=side_specs,
    )

    feature_columns = [
        column
        for column in sample_summary.columns
        if column not in {"region_id", "sample_id", "condition", pairing_key}
    ]
    if contrast.mode == "matched_pairwise":
        _validate_single_read_matched_sample_summary(
            sample_summary,
            pairing_key=pairing_key or "",
        )
        pair_columns = ["region_id", pairing_key, "condition"]
        pair_summary = sample_summary.groupby(pair_columns, as_index=False, sort=False)[
            feature_columns
        ].mean()
        pair_summary = _select_complete_matched_pairs(
            pair_summary,
            pairing_key=pairing_key or "",
            contrast=contrast,
        )

        numerator = (
            pair_summary.loc[pair_summary["condition"].isin(side_specs["numerator"])]
            .drop(columns=["condition"])
            .rename(
                columns={
                    feature_name: f"{feature_name}_numerator_mean"
                    for feature_name in feature_columns
                }
            )
        )
        denominator = (
            pair_summary.loc[pair_summary["condition"].isin(side_specs["denominator"])]
            .drop(columns=["condition"])
            .rename(
                columns={
                    feature_name: f"{feature_name}_denominator_mean"
                    for feature_name in feature_columns
                }
            )
        )
        regions_table = numerator.merge(
            denominator,
            on=["region_id", pairing_key],
            how="inner",
            sort=False,
        )
        if regions_table.empty:
            raise ValueError(
                "single_read matched_pairwise scoring found no complete matched pairs."
            )

        summary_rows = []
        for region_id, region_table in regions_table.groupby("region_id", sort=False):
            row = {"region_id": region_id}
            effect_sizes = []
            for feature_name in feature_columns:
                numerator_mean = float(
                    region_table[f"{feature_name}_numerator_mean"].mean()
                )
                denominator_mean = float(
                    region_table[f"{feature_name}_denominator_mean"].mean()
                )
                delta_mean = numerator_mean - denominator_mean
                row[f"{feature_name}_numerator_mean"] = numerator_mean
                row[f"{feature_name}_denominator_mean"] = denominator_mean
                row[f"{feature_name}_delta_mean"] = delta_mean
                effect_sizes.append(abs(delta_mean))
            row["effect_size"] = max(effect_sizes, default=0.0)
            summary_rows.append(row)

        summary = pd.DataFrame(summary_rows)
        if summary.empty:
            summary = pd.DataFrame(columns=["region_id", "effect_size"])
        else:
            summary = summary.sort_values(
                by=["effect_size", "region_id"],
                ascending=[False, True],
                kind="mergesort",
            ).reset_index(drop=True)
        summary["rank"] = range(1, len(summary) + 1)
        return evidence, summary.drop(columns="effect_size", errors="ignore")

    summary_rows = []
    for region_id, region_table in sample_summary.groupby("region_id", sort=False):
        numerator_table = region_table[
            region_table["condition"].isin(contrast.numerator)
        ]
        denominator_table = region_table[
            region_table["condition"].isin(contrast.denominator)
        ]
        row = {"region_id": region_id}
        effect_sizes = []
        for feature_name in feature_columns:
            numerator_mean = float(numerator_table[feature_name].mean())
            denominator_mean = float(denominator_table[feature_name].mean())
            delta_mean = numerator_mean - denominator_mean
            row[f"{feature_name}_numerator_mean"] = numerator_mean
            row[f"{feature_name}_denominator_mean"] = denominator_mean
            row[f"{feature_name}_delta_mean"] = delta_mean
            effect_sizes.append(abs(delta_mean))
        row["effect_size"] = max(effect_sizes, default=0.0)
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    if summary.empty:
        summary = pd.DataFrame(columns=["region_id", "effect_size"])
    else:
        summary = summary.sort_values(
            by=["effect_size", "region_id"],
            ascending=[False, True],
            kind="mergesort",
        ).reset_index(drop=True)
    summary["rank"] = range(1, len(summary) + 1)
    return evidence, summary.drop(columns="effect_size", errors="ignore")


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
    # STATISTICS CHANGE (review fix #1): only rows with >= 2 replicates per side
    # are genuinely testable by Welch's t-test; single-replicate rows return a
    # degenerate p=1.0 and must not be counted in the BH total.
    testable = scored.apply(
        lambda row: len(list(row["numerator_sample_values"])) >= 2
        and len(list(row["denominator_sample_values"])) >= 2,
        axis=1,
    )
    scored["adjusted_p_value"] = _adjust_p_values_bh(
        scored["p_value"], testable=testable
    )
    return scored


def _dominant_label(values: pd.Series) -> str | None:
    non_null = values.dropna()
    if non_null.empty:
        return None
    counts = non_null.value_counts()
    max_count = counts.max()
    return sorted(counts[counts == max_count].index.tolist())[0]


def _dominant_clusters_differ(
    dominant_cluster: str | None,
    reference_dominant_cluster: str | None,
) -> bool:
    dominant_missing = pd.isna(dominant_cluster)
    reference_missing = pd.isna(reference_dominant_cluster)
    if dominant_missing and reference_missing:
        return False
    if dominant_missing or reference_missing:
        return True
    return bool(dominant_cluster != reference_dominant_cluster)


def _score_cluster_occupancy(
    *,
    evidence: pd.DataFrame,
    contrast: ContrastSpec,
    representation: str,
    test: str,
    multiple_testing: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pseudocount = 1e-6
    _require_supported_cluster_occupancy_mode(contrast)

    if representation == "cluster_fraction":
        _require_occupancy_columns(
            evidence,
            required_columns={
                "region_id",
                "sample_id",
                "condition",
                "cluster",
                "fraction",
            },
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
        regions_table["reference_fraction"] = regions_table[
            "reference_fraction"
        ].fillna(0.0)
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
                by=[
                    "adjusted_p_value",
                    "p_value",
                    "effect_size",
                    "region_id",
                    "cluster",
                ],
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
            required_columns={
                "region_id",
                "sample_id",
                "condition",
                "dominant_cluster",
            },
        )
        side_specs = {
            "numerator": contrast.numerator or [],
            "denominator": contrast.denominator or [],
        }
        _require_requested_conditions_present(
            evidence,
            side_specs=side_specs,
        )
        sample_level = evidence.loc[
            :, ["region_id", "sample_id", "condition", "dominant_cluster"]
        ].drop_duplicates()
        numerator = (
            sample_level.loc[sample_level["condition"].isin(side_specs["numerator"])]
            .groupby("region_id", dropna=False, sort=False)
            .agg(
                dominant_cluster=("dominant_cluster", _dominant_label),
                numerator_replicate_n=("sample_id", "nunique"),
            )
            .reset_index()
        )
        denominator = (
            sample_level.loc[sample_level["condition"].isin(side_specs["denominator"])]
            .groupby("region_id", dropna=False, sort=False)
            .agg(
                reference_dominant_cluster=("dominant_cluster", _dominant_label),
                denominator_replicate_n=("sample_id", "nunique"),
            )
            .reset_index()
        )
        regions_table = numerator.merge(
            denominator,
            on="region_id",
            how="outer",
            sort=False,
        )
        missing_region_side = (
            regions_table["numerator_replicate_n"].isna()
            ^ regions_table["denominator_replicate_n"].isna()
        )
        if missing_region_side.any():
            raise ValueError(
                "Missing dominant_cluster evidence for one or more requested contrast "
                "sides at some regions."
            )
        regions_table["numerator_replicate_n"] = (
            regions_table["numerator_replicate_n"].fillna(0).astype(int)
        )
        regions_table["denominator_replicate_n"] = (
            regions_table["denominator_replicate_n"].fillna(0).astype(int)
        )
        regions_table["dominant_cluster_changed"] = regions_table.apply(
            lambda row: _dominant_clusters_differ(
                row["dominant_cluster"],
                row["reference_dominant_cluster"],
            ),
            axis=1,
        )
        regions_table["effect_size"] = regions_table["dominant_cluster_changed"].astype(
            int
        )
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
        evidence.loc[
            :, ["region_id", "sample_id", "condition", "cluster_entropy"]
        ].drop_duplicates(),
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
    read_table: pd.DataFrame | None = None,
    feature_table: pd.DataFrame | None = None,
) -> RegionContrastResult:
    validate_region_contrast_request(
        analysis_unit=analysis_unit,
        representation=representation,
        signal_source=signal_source,
        test=test,
    )
    contrast_metadata = contrast.metadata or {}

    if analysis_unit == "single_read":
        if representation == "read_mod_fraction":
            if read_table is None:
                raise ValueError(
                    "single_read read_mod_fraction scoring currently requires read_table."
                )
            evidence = build_single_read_mod_fraction_evidence_table(
                extract_table=read_table,
                pairing_key=contrast.pairing_key
                if contrast.mode == "matched_pairwise"
                else None,
            )
            regions_table, summary = _score_single_read_mod_fraction(
                evidence=evidence,
                contrast=contrast,
            )
            metadata = {
                "contrast_mode": contrast.mode,
                "analysis_unit": analysis_unit,
                "representation": representation,
                "signal_source": signal_source,
                "test": test,
                "multiple_testing": multiple_testing,
                "normalization_mode": contrast_metadata.get(
                    "normalization_mode", "none"
                ),
                "biological_interpretation": contrast_metadata.get(
                    "biological_interpretation",
                    "region-level difference in sample-level read modification fraction",
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

        if representation == "read_window_features":
            if feature_table is None:
                feature_table = _load_builtin_single_read_feature_table(
                    samples=samples,
                    regions=regions,
                    motifs=motifs,
                )
            evidence = build_single_read_feature_evidence_table(
                feature_table=feature_table,
                pairing_key=contrast.pairing_key
                if contrast.mode == "matched_pairwise"
                else None,
            )
            regions_table, summary = _score_single_read_features(
                evidence=evidence,
                contrast=contrast,
            )
            metadata = {
                "contrast_mode": contrast.mode,
                "analysis_unit": analysis_unit,
                "representation": representation,
                "signal_source": signal_source,
                "test": test,
                "multiple_testing": multiple_testing,
                "normalization_mode": contrast_metadata.get(
                    "normalization_mode", "none"
                ),
                "biological_interpretation": contrast_metadata.get(
                    "biological_interpretation",
                    "region-level difference in sample-level read feature summaries",
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

        raise NotImplementedError(
            f"Single-read representation '{representation}' is not implemented yet."
        )

    if analysis_unit == "cluster_occupancy":
        if occupancy_table is None:
            raise ValueError(
                "cluster_occupancy scoring currently requires occupancy_table."
            )
        _validate_occupancy_table_numeric_columns(occupancy_table)
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

    # STATISTICS CHANGE (review fix #3): record which side actually contributed a
    # row before we fillna(0). An outer merge fabricates a zero reference (or a
    # zero numerator) for regions present on only one side; without these flags
    # such rows report a huge log2fc (~+/-20) and can be ranked as top hits even
    # though the "effect" is an artifact of missing coverage on the other side.
    merged["numerator_present"] = merged["numerator_valid_count"].notna()
    merged["denominator_present"] = merged["denominator_valid_count"].notna()

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
        (merged["fraction"] + pseudocount)
        / (merged["reference_fraction"] + pseudocount)
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

    # STATISTICS CHANGE (review fix #3): a row backed by both sides is a genuine
    # contrast; a row present on only one side is a fabricated-zero comparison.
    merged["both_sides_present"] = (
        merged["numerator_present"] & merged["denominator_present"]
    )

    if test == "beta_binomial":
        regions_table = _add_beta_binomial_scores(
            merged,
            multiple_testing=multiple_testing,
            replicate_counts=_replicate_counts_by_region(evidence, contrast),
        )
        regions_table = regions_table.sort_values(
            by=["adjusted_p_value", "p_value"],
            ascending=[True, True],
            kind="mergesort",
        ).reset_index(drop=True)
    else:
        # Rank genuine (both-sides-present) rows ahead of fabricated-zero rows so a
        # zero-coverage row is never silently surfaced as the top hit.
        regions_table = merged.sort_values(
            by=["both_sides_present", "effect_size"],
            ascending=[False, False],
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
        "numerator_present",
        "denominator_present",
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
        summary_columns.extend(
            ["count", "reference_count", "delta_count", "log2_fc_count"]
        )
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
