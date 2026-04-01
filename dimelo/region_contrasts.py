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
    if analysis_unit != "ensemble_region":
        raise ValueError(
            "V1 region_contrasts inference requires analysis_unit='ensemble_region'."
        )
    if signal_source != "pileup_counts":
        raise ValueError(
            "V1 region_contrasts inference requires signal_source='pileup_counts'."
        )
    if representation not in {"modified_fraction", "modified_count"}:
        raise ValueError(
            "V1 region_contrasts inference requires representation to be "
            "'modified_fraction' or 'modified_count'."
        )
    if test not in {"beta_binomial", "effect_size_only"}:
        raise ValueError(
            "V1 region_contrasts inference requires test to be "
            "'beta_binomial' or 'effect_size_only'."
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


def _zero_safe_divide(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    return numerator.div(denominator.where(denominator != 0), fill_value=0).fillna(0.0)


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
    test: str = "beta_binomial",
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

    if test != "effect_size_only":
        raise NotImplementedError(f"score_regions does not yet implement test='{test}'.")

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
    merged["abs_delta_fraction"] = merged["delta_fraction"].abs()

    integer_columns = [
        "numerator_modified_count",
        "numerator_valid_count",
        "numerator_replicate_n",
        "denominator_modified_count",
        "denominator_valid_count",
        "denominator_replicate_n",
    ]
    merged[integer_columns] = merged[integer_columns].astype(int)

    regions_table = merged.sort_values(
        by="abs_delta_fraction",
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
