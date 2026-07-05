"""Integrative population / genome-architecture views (Q10).

The capstone: combine the per-site and per-read outputs of Q2-Q8 into a fuller picture of
the cell population and genome architecture. This is synthesis (joins, cross-tabs), not a
new statistic.

- ``build_integrated_site_table`` joins the per-site metrics from the other workflows
  (occupancy/S2, binding strength/S3, a per-site footprint summary/S4, a per-site
  joint-occupancy summary/S5, external-track & Hi-C signal/S6, trans-correction/S7) into
  one denormalized per-site table — the backbone for ranking, filtering, and figures.
- ``single_molecule_state_composition`` cross-tabulates per-read annotations
  (bound/unbound, footprint present, co-occupancy state) into the fraction of single
  molecules in each combinatorial chromatin state per site — the cell-population
  heterogeneity view long reads uniquely enable.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from functools import reduce

import pandas as pd

_SITE_SOURCES = (
    "occupancy",
    "binding_strength",
    "footprint",
    "joint_occupancy",
    "tracks",
    "trans",
)


def build_integrated_site_table(
    *,
    site_column: str = "region_id",
    occupancy: pd.DataFrame | None = None,
    binding_strength: pd.DataFrame | None = None,
    footprint: pd.DataFrame | None = None,
    joint_occupancy: pd.DataFrame | None = None,
    tracks: pd.DataFrame | None = None,
    trans: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Outer-join per-site metric frames into one denormalized table keyed on ``site_column``.

    Each source is an optional per-site ``DataFrame`` that must contain ``site_column``.
    Columns that appear in more than one source are prefixed with the source name (e.g.
    ``occupancy_n_reads`` vs ``binding_strength_n_reads``); columns unique to one source keep
    their name. Requires at least one non-``None`` source.
    """
    sources = {
        "occupancy": occupancy,
        "binding_strength": binding_strength,
        "footprint": footprint,
        "joint_occupancy": joint_occupancy,
        "tracks": tracks,
        "trans": trans,
    }
    provided = {name: frame for name, frame in sources.items() if frame is not None}
    if not provided:
        raise ValueError(
            "build_integrated_site_table requires at least one non-None source."
        )
    for name, frame in provided.items():
        if site_column not in frame.columns:
            raise ValueError(
                f"integrated source {name!r} requires the site column {site_column!r}."
            )

    # prefix only columns that collide across sources
    column_counts: Counter[str] = Counter()
    for frame in provided.values():
        for column in frame.columns:
            if column != site_column:
                column_counts[column] += 1
    colliding = {column for column, count in column_counts.items() if count > 1}

    frames: list[pd.DataFrame] = []
    for name, frame in provided.items():
        renamed = frame.rename(
            columns={
                column: f"{name}_{column}"
                for column in frame.columns
                if column in colliding
            }
        )
        frames.append(renamed)

    return reduce(
        lambda left, right: left.merge(right, on=site_column, how="outer"), frames
    ).reset_index(drop=True)


def single_molecule_state_composition(
    read_states: pd.DataFrame,
    *,
    state_columns: Sequence[str],
    site_column: str = "region_id",
) -> pd.DataFrame:
    """Per-site composition of combinatorial single-molecule states.

    ``read_states`` is a per-read frame with ``site_column`` and the given ``state_columns``
    (e.g. ``is_true_signal``, ``footprint_present``, ``joint_state`` — assemble by joining
    the per-read outputs on ``read_id``). The state columns are combined into a
    ``combined_state`` label and the fraction of reads in each state is computed per site.
    Returns ``site_column, combined_state, n_reads, fraction`` (fractions sum to 1 per site).
    """
    state_columns = list(state_columns)
    if not state_columns:
        raise ValueError("state_columns must be non-empty.")
    missing = ({site_column, *state_columns}) - set(read_states.columns)
    if missing:
        raise ValueError(
            "single_molecule_state_composition requires columns: "
            f"{', '.join(sorted(missing))}."
        )
    columns = [site_column, "combined_state", "n_reads", "fraction"]
    if read_states.empty:
        return pd.DataFrame(columns=columns)

    frame = read_states.copy()
    frame["combined_state"] = (
        frame[state_columns].astype(str).agg("|".join, axis=1)
    )
    grouped = (
        frame.groupby([site_column, "combined_state"], sort=False)
        .size()
        .reset_index(name="n_reads")
    )
    totals = grouped.groupby(site_column)["n_reads"].transform("sum")
    grouped["fraction"] = grouped["n_reads"] / totals.where(totals != 0)
    return grouped.loc[:, columns]
