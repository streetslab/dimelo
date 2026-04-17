from __future__ import annotations

from functools import partial
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd
import pysam
from tqdm.auto import tqdm

from . import load_processed, utils
from .models import GlobalAnalysisResult, SampleSpec


def _global_counts_for_motifs_from_bedmethyl(
    bedmethyl_file: str | Path,
    motifs: Sequence[str],
    quiet: bool = True,
) -> dict[str, tuple[int, int]]:
    """
    Sum modified and valid pileup counts across all contigs for all motifs in one pass.
    """

    del quiet

    unique_motifs = list(dict.fromkeys(motifs))
    if len(unique_motifs) == 0:
        return {}

    parsed_motifs = {
        motif: utils.ParsedMotif(motif) for motif in unique_motifs
    }
    counts_by_motif = {motif: [0, 0] for motif in unique_motifs}

    with pysam.TabixFile(str(bedmethyl_file)) as tabix_file:
        for contig in tabix_file.contigs:
            for row in tabix_file.fetch(contig):
                for motif, parsed_motif in parsed_motifs.items():
                    keep_basemod, _, modified_in_row, valid_in_row = (
                        load_processed.process_pileup_row(
                            row=row,
                            parsed_motif=parsed_motif,
                            region_strand=".",
                        )
                    )
                    if keep_basemod:
                        counts_by_motif[motif][0] += modified_in_row
                        counts_by_motif[motif][1] += valid_in_row

    return {
        motif: (modified_count, valid_count)
        for motif, (modified_count, valid_count) in counts_by_motif.items()
    }


def _global_counts_from_bedmethyl(
    bedmethyl_file: str | Path,
    motif: str,
    quiet: bool = True,
) -> tuple[int, int]:
    """
    Sum modified and valid pileup counts across all contigs for one motif.
    """

    counts_by_motif = _global_counts_for_motifs_from_bedmethyl(
        bedmethyl_file=bedmethyl_file,
        motifs=[motif],
        quiet=quiet,
    )
    return counts_by_motif[motif]


def summarize_global_samples(
    samples: Sequence[SampleSpec],
    motifs: Sequence[str],
    quiet: bool = True,
) -> pd.DataFrame:
    """
    Build a long-form table of global motif fractions per sample.
    """

    rows: list[dict[str, object]] = []

    for sample in tqdm(samples, desc="Summarizing global samples", disable=quiet):
        if not sample.metadata or "pileup_path" not in sample.metadata:
            raise ValueError(
                f"Sample {sample.sample_id!r} is missing metadata['pileup_path'] for global_analysis mode."
            )

        pileup_path = sample.metadata["pileup_path"]
        counts_by_motif = _global_counts_for_motifs_from_bedmethyl(
            bedmethyl_file=pileup_path,
            motifs=motifs,
            quiet=quiet,
        )
        for motif in motifs:
            modified_count, valid_count = counts_by_motif[motif]
            rows.append(
                {
                    "sample_id": sample.sample_id,
                    "condition": sample.condition,
                    "replicate": sample.replicate,
                    "motif": motif,
                    "modified_count": modified_count,
                    "valid_count": valid_count,
                    "global_fraction": (
                        modified_count / valid_count if valid_count else float("nan")
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


def tile_genome_windows(
    *,
    genome_sizes: dict[str, int],
    window_size: int,
    step_size: int,
    include_contigs: Iterable[str] | None = None,
    exclude_contigs: Iterable[str] | None = None,
) -> pd.DataFrame:
    if window_size <= 0:
        raise ValueError("window_size must be positive.")
    if step_size <= 0:
        raise ValueError("step_size must be positive.")

    include_set = set(include_contigs) if include_contigs is not None else None
    exclude_set = set(exclude_contigs or [])

    rows: list[dict[str, object]] = []
    for chromosome, contig_length in genome_sizes.items():
        if include_set is not None and chromosome not in include_set:
            continue
        if chromosome in exclude_set:
            continue
        if contig_length <= 0:
            continue

        start = 0
        while start < contig_length:
            end = min(start + window_size, contig_length)
            rows.append(
                {
                    "window_id": f"{chromosome}:{start}-{end}",
                    "chromosome": chromosome,
                    "start": start,
                    "end": end,
                    "strand": ".",
                }
            )
            if end >= contig_length:
                break
            start += step_size

    return pd.DataFrame(
        rows,
        columns=["window_id", "chromosome", "start", "end", "strand"],
    )


def build_window_summary(
    *,
    samples: Sequence[SampleSpec],
    motifs: Sequence[str],
    genome_sizes: dict[str, int],
    window_size: int,
    step_size: int,
    include_contigs: Iterable[str] | None = None,
    exclude_contigs: Iterable[str] | None = None,
    quiet: bool = True,
    cores: int | None = None,
) -> pd.DataFrame:
    motifs = list(motifs)
    if len(motifs) != 1:
        raise ValueError("build_window_summary currently supports exactly one motif.")

    motif = motifs[0]
    windows = tile_genome_windows(
        genome_sizes=genome_sizes,
        window_size=window_size,
        step_size=step_size,
        include_contigs=include_contigs,
        exclude_contigs=exclude_contigs,
    )

    region_strings = [
        f"{row.chromosome}:{row.start}-{row.end},{row.strand}"
        for row in windows.itertuples(index=False)
    ]

    rows: list[dict[str, object]] = []
    for sample in tqdm(samples, desc="Processing tiled windows", disable=quiet):
        metadata = sample.metadata or {}
        if "pileup_path" not in metadata:
            raise ValueError(
                f"Sample {sample.sample_id!r} is missing metadata['pileup_path'] for global_analysis mode."
            )

        pileup_loader = partial(
            load_processed.pileup_counts_from_bedmethyl,
            bedmethyl_file=metadata["pileup_path"],
            motif=motif,
        )
        counts_by_window = load_processed.regions_to_list(
            function_handle=pileup_loader,
            regions=region_strings,
            window_size=None,
            quiet=quiet,
            cores=cores,
        )

        if len(counts_by_window) != len(windows):
            raise RuntimeError(
                "Pileup window count length did not match the number of tiled windows."
            )

        for window_row, (modified_count, valid_count) in zip(
            windows.itertuples(index=False),
            counts_by_window,
        ):
            window_fraction = (
                0.0 if valid_count == 0 else modified_count / valid_count
            )
            rows.append(
                {
                    "sample_id": sample.sample_id,
                    "condition": sample.condition,
                    "replicate": sample.replicate,
                    "motif": motif,
                    "window_id": window_row.window_id,
                    "chromosome": window_row.chromosome,
                    "start": window_row.start,
                    "end": window_row.end,
                    "strand": window_row.strand,
                    "modified_count": modified_count,
                    "valid_count": valid_count,
                    "window_fraction": window_fraction,
                }
            )

    return pd.DataFrame(
        rows,
        columns=[
            "sample_id",
            "condition",
            "replicate",
            "motif",
            "window_id",
            "chromosome",
            "start",
            "end",
            "strand",
            "modified_count",
            "valid_count",
            "window_fraction",
        ],
    )


def compute_global_normalization_factors(summary: pd.DataFrame) -> pd.DataFrame:
    factors = summary.copy()
    factors["reference_fraction"] = factors.groupby("motif", sort=True)[
        "global_fraction"
    ].transform("mean")
    factors["global_offset"] = (
        factors["global_fraction"] - factors["reference_fraction"]
    )
    return factors.loc[
        :,
        [
            "sample_id",
            "condition",
            "replicate",
            "motif",
            "global_fraction",
            "reference_fraction",
            "global_offset",
        ],
    ].copy()


def run_global_analysis(
    *,
    samples: Sequence[SampleSpec],
    motifs: Iterable[str],
    genome_sizes: dict[str, int],
    window_size: int,
    step_size: int,
    include_contigs: Iterable[str] | None = None,
    exclude_contigs: Iterable[str] | None = None,
    quiet: bool = True,
    cores: int | None = None,
) -> GlobalAnalysisResult:
    motifs = list(motifs)
    if len(motifs) == 0:
        raise ValueError("run_global_analysis requires at least one motif.")

    summary = summarize_global_samples(samples=samples, motifs=motifs, quiet=quiet)
    per_motif_windows = [
        build_window_summary(
            samples=samples,
            motifs=[motif],
            genome_sizes=genome_sizes,
            window_size=window_size,
            step_size=step_size,
            include_contigs=include_contigs,
            exclude_contigs=exclude_contigs,
            quiet=quiet,
            cores=cores,
        )
        for motif in motifs
    ]
    windows = pd.concat(per_motif_windows, ignore_index=True)
    normalization_factors = compute_global_normalization_factors(summary)
    plot_data = {
        "global_fraction_bar": summary.loc[
            :,
            ["sample_id", "condition", "motif", "global_fraction"],
        ].copy(),
        "window_fraction_table": windows.copy(),
    }

    return GlobalAnalysisResult(
        summary=summary,
        windows=windows,
        normalization_factors=normalization_factors,
        plot_data=plot_data,
        metadata={
            "normalization_mode": "per_sample_global",
            "window_size": window_size,
            "step_size": step_size,
            "motifs": motifs,
        },
        figures={},
    )
