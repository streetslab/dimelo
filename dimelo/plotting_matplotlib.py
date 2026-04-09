from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import pandas as pd


def save_figure(
    fig,
    path,
    *,
    dpi: int = 300,
    bbox_inches: str = "tight",
    transparent: bool = False,
) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output_path,
        dpi=dpi,
        bbox_inches=bbox_inches,
        transparent=transparent,
    )
    return output_path


def _require_payload_table(payload: Mapping[str, object], key: str) -> pd.DataFrame:
    if key not in payload:
        raise ValueError(f"plot payload is missing required key: {key}")
    table = payload[key]
    if not isinstance(table, pd.DataFrame):
        raise TypeError(f"plot payload key {key!r} must be a pandas DataFrame.")
    return table


def _import_pyplot():
    from matplotlib import pyplot as plt

    return plt


def plot_region_contrast_profile_matplotlib(payload: Mapping[str, object]):
    plot_table = _require_payload_table(payload, "plot_table")
    metadata = payload.get("metadata", {}) if isinstance(payload.get("metadata", {}), Mapping) else {}

    plt = _import_pyplot()
    fig, ax = plt.subplots(figsize=(8, 4))

    if plot_table.empty:
        ax.set_title("Region contrast profile")
        ax.set_xlabel("position")
        ax.set_ylabel("value")
        return fig, ax

    x_column = "plot_x" if "plot_x" in plot_table.columns else "position"
    y_column = "value"

    if "value_mode" in plot_table.columns:
        for value_mode, subset in plot_table.groupby("value_mode", sort=False):
            subset = subset.sort_values([x_column, "region_id"], kind="stable")
            ax.plot(subset[x_column], subset[y_column], marker="o", linewidth=1.5, label=str(value_mode))
        ax.legend(title="value_mode")
    else:
        plot_table = plot_table.sort_values([x_column, "region_id"], kind="stable")
        ax.plot(plot_table[x_column], plot_table[y_column], marker="o", linewidth=1.5)

    axis_table = payload.get("axis_table")
    if isinstance(axis_table, pd.DataFrame) and not axis_table.empty and {"axis_min", "axis_max"}.issubset(axis_table.columns):
        ax.set_xlim(axis_table["axis_min"].min(), axis_table["axis_max"].max())

    ax.set_title("Region contrast profile")
    ax.set_xlabel("position" if x_column == "position" else x_column)
    ax.set_ylabel("value")

    if metadata.get("orientation") == "region_5to3":
        ax.set_xlabel("oriented position")

    return fig, ax


def plot_region_contrast_heatmap_matplotlib(payload: Mapping[str, object]):
    plot_table = _require_payload_table(payload, "plot_table")
    summary_table = _require_payload_table(payload, "summary_table")

    plt = _import_pyplot()
    fig, ax = plt.subplots(figsize=(8, 4))

    if plot_table.empty:
        ax.set_title("Region contrast heatmap")
        ax.set_xlabel("position")
        ax.set_ylabel("region")
        return fig, ax

    x_column = "plot_x" if "plot_x" in plot_table.columns else "position"
    y_column = "row_order" if "row_order" in plot_table.columns else "region_id"

    heatmap_table = plot_table.loc[:, [x_column, y_column, "value"]].copy()
    heatmap = heatmap_table.pivot_table(
        index=y_column,
        columns=x_column,
        values="value",
        aggfunc="mean",
    ).sort_index(axis=0).sort_index(axis=1)

    image = ax.imshow(
        heatmap.to_numpy(),
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    ax.figure.colorbar(image, ax=ax, label="value")

    if y_column == "row_order" and "region_id" in summary_table.columns:
        summary_lookup = (
            summary_table.loc[:, ["region_id", "row_order"]]
            .drop_duplicates()
            .sort_values("row_order", kind="stable")
        )
        ax.set_yticks(range(len(summary_lookup)))
        ax.set_yticklabels(summary_lookup["region_id"].astype(str).tolist())
    else:
        ax.set_ylabel("region")

    ax.set_xticks(range(len(heatmap.columns)))
    ax.set_xticklabels([str(value) for value in heatmap.columns.tolist()], rotation=45, ha="right")
    ax.set_title("Region contrast heatmap")
    ax.set_xlabel("position" if x_column == "position" else x_column)

    return fig, ax
