from __future__ import annotations

import concurrent.futures
from pathlib import Path
from typing import Sequence

import numpy as np
from tqdm.auto import tqdm

from . import cluster, utils
from .models import SampleSpec


def _build_sample_region_features(
    *,
    sample: SampleSpec,
    motif: str,
    matched_regions: str | Path | list[str | Path] | None,
    pileup_paths: dict[str, str | Path] | None,
    cores: int | None,
    regions_executor: concurrent.futures.Executor | None = None,
) -> tuple[np.ndarray, list[dict[str, object]]]:
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
        motif=motif,
        regions=regions,
        cores=cores,
        quiet=True,
        regions_executor=regions_executor,
    )

    matrix_float32 = np.asarray(matrix, dtype=np.float32)
    metadata_rows = [
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
        for chrom, start, end, strand in region_metadata
    ]
    return matrix_float32, metadata_rows


def build_region_feature_table(
    *,
    samples: Sequence[SampleSpec],
    motifs: Sequence[str],
    matched_regions: str | Path | list[str | Path] | None,
    pileup_paths: dict[str, str | Path] | None = None,
    cores: int | None = None,
    quiet: bool = True,
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

    sample_list = list(samples)
    if not sample_list:
        raise ValueError("No region feature rows were generated.")

    if cores is None:
        sample_workers = 1
        per_sample_cores: int | None = None
    else:
        cores_to_run = max(1, utils.cores_to_run(cores))
        sample_workers = min(len(sample_list), cores_to_run)
        per_sample_cores = max(1, cores_to_run // sample_workers)

    if sample_workers > 1:
        ordered_results: list[tuple[np.ndarray, list[dict[str, object]]] | None] = [
            None
        ] * len(sample_list)
        with concurrent.futures.ThreadPoolExecutor(max_workers=sample_workers) as executor:
            future_by_index = {
                index: executor.submit(
                    _build_sample_region_features,
                    sample=sample,
                    motif=motifs[0],
                    matched_regions=matched_regions,
                    pileup_paths=pileup_paths,
                    cores=per_sample_cores,
                    regions_executor=None,
                )
                for index, sample in enumerate(sample_list)
            }
            for index, future in future_by_index.items():
                ordered_results[index] = future.result()
        results = [result for result in ordered_results if result is not None]
    else:
        shared_regions_executor = None
        if len(sample_list) > 1 and (cores is None or cores > 1):
            try:
                shared_regions_executor = concurrent.futures.ProcessPoolExecutor(
                    max_workers=max(1, utils.cores_to_run(cores)),
                )
            except Exception:
                # Some environments (for example restricted CI sandboxes) cannot
                # create process pools; fall back to the legacy per-sample path.
                shared_regions_executor = None
        try:
            results = []
            for sample in tqdm(sample_list, desc="Building region features", disable=quiet):
                result = _build_sample_region_features(
                    sample=sample,
                    motif=motifs[0],
                    matched_regions=matched_regions,
                    pileup_paths=pileup_paths,
                    cores=per_sample_cores,
                    regions_executor=shared_regions_executor,
                )
                results.append(result)
        finally:
            if shared_regions_executor is not None:
                shared_regions_executor.shutdown(wait=True, cancel_futures=False)

    matrices = [matrix for matrix, _ in results]
    metadata_rows = [row for _, rows in results for row in rows]

    if not matrices:
        raise ValueError("No region feature rows were generated.")
    return np.vstack(matrices).astype(np.float32, copy=False), metadata_rows
