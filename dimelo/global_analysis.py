from __future__ import annotations

from pathlib import Path
from typing import Sequence

import pandas as pd
import pysam

from . import load_processed, utils
from .models import SampleSpec


def _global_counts_from_bedmethyl(
    bedmethyl_file: str | Path,
    motif: str,
    quiet: bool = True,
) -> tuple[int, int]:
    """
    Sum modified and valid pileup counts across all contigs for one motif.
    """

    del quiet

    parsed_motif = utils.ParsedMotif(motif)
    modified_count = 0
    valid_count = 0

    with pysam.TabixFile(str(bedmethyl_file)) as tabix_file:
        for contig in tabix_file.contigs:
            for row in tabix_file.fetch(contig):
                keep_basemod, _, modified_in_row, valid_in_row = (
                    load_processed.process_pileup_row(
                        row=row,
                        parsed_motif=parsed_motif,
                        region_strand=".",
                    )
                )
                if keep_basemod:
                    modified_count += modified_in_row
                    valid_count += valid_in_row

    return modified_count, valid_count


def summarize_global_samples(
    samples: Sequence[SampleSpec],
    motifs: Sequence[str],
    quiet: bool = True,
) -> pd.DataFrame:
    """
    Build a long-form table of global motif fractions per sample.
    """

    rows: list[dict[str, object]] = []

    for sample in samples:
        if not sample.metadata or "pileup_path" not in sample.metadata:
            raise ValueError(
                f"Sample {sample.sample_id!r} is missing metadata['pileup_path'] for global_analysis mode."
            )

        pileup_path = sample.metadata["pileup_path"]
        for motif in motifs:
            modified_count, valid_count = _global_counts_from_bedmethyl(
                bedmethyl_file=pileup_path,
                motif=motif,
                quiet=quiet,
            )
            rows.append(
                {
                    "sample_id": sample.sample_id,
                    "condition": sample.condition,
                    "replicate": sample.replicate,
                    "motif": motif,
                    "modified_count": modified_count,
                    "valid_count": valid_count,
                    "global_fraction": (
                        modified_count / valid_count if valid_count else 0.0
                    ),
                }
            )

    return pd.DataFrame(
        rows,
        columns=[
            "sample_id",
            "condition",
            "replicate",
            "motif",
            "modified_count",
            "valid_count",
            "global_fraction",
        ],
    )
