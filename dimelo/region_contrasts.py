from __future__ import annotations

from functools import partial
from typing import Iterable

import pandas as pd

from . import load_processed, utils


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
