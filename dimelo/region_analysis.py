from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np

from . import cluster
from .models import SampleSpec


def build_region_feature_table(
    *,
    samples: Sequence[SampleSpec],
    motifs: Sequence[str],
    matched_regions: str | Path | list[str | Path] | None,
    pileup_paths: dict[str, str | Path] | None = None,
) -> tuple[np.ndarray, list[dict[str, object]]]:
    """
    Build one region-level feature row per sample and matched region.
    """

    if not motifs:
        raise ValueError("motifs must contain at least one motif.")
    if len(motifs) != 1:
        raise ValueError(
            "build_region_feature_table currently supports exactly one motif for "
            "region_anchored mode."
        )

    matrices: list[np.ndarray] = []
    metadata_rows: list[dict[str, object]] = []

    for sample in samples:
        pileup_path = None
        if pileup_paths is not None:
            pileup_path = pileup_paths.get(sample.sample_id)
        if pileup_path is None:
            if not sample.metadata or "pileup_path" not in sample.metadata:
                raise ValueError(
                    f"Sample {sample.sample_id!r} is missing a pileup_path for region_anchored mode."
                )
            pileup_path = sample.metadata["pileup_path"]

        regions = matched_regions or sample.regions_bed
        if regions is None:
            raise ValueError(
                f"Sample {sample.sample_id!r} requires matched_regions or regions_bed."
            )

        matrix, region_metadata = cluster.region_feature_matrix_from_pileup(
            bedmethyl_file=pileup_path,
            motif=motifs[0],
            regions=regions,
        )
        matrices.append(np.asarray(matrix, dtype=np.float32))
        for chrom, start, end, strand in region_metadata:
            metadata_rows.append(
                {
                    "sample_id": sample.sample_id,
                    "condition": sample.condition,
                    "replicate": sample.replicate,
                    "region_id": f"{chrom}:{start}-{end}:{strand}",
                    "chromosome": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                }
            )

    if not matrices:
        raise ValueError("No region feature rows were generated.")
    return np.vstack(matrices).astype(np.float32, copy=False), metadata_rows
