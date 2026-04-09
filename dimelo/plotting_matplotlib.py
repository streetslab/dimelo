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


def _require_payload_keys(payload: Mapping[str, object], keys: tuple[str, ...], *, owner: str) -> None:
    missing = [key for key in keys if key not in payload]
    if missing:
        raise ValueError(f"{owner} payload is missing required keys: {', '.join(missing)}")


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


def _make_axis(*, ax=None, figsize=(8, 4)):
    plt = _import_pyplot()
    if ax is not None:
        return ax.figure, ax
    fig, created_ax = plt.subplots(figsize=figsize)
    return fig, created_ax


def _prepare_region_contrast_value_mode_table(
    plot_table: pd.DataFrame,
    *,
    value_mode: str,
) -> pd.DataFrame:
    if value_mode not in {"numerator", "denominator", "delta"}:
        raise ValueError("Unsupported region contrast value_mode.")

    if "value_mode" not in plot_table.columns:
        return plot_table.copy()

    filtered = plot_table.loc[plot_table["value_mode"] == value_mode].copy()
    if filtered.empty:
        raise ValueError(f"plot payload does not contain rows for value_mode={value_mode!r}.")
    return filtered


def plot_region_contrast_profile_matplotlib(
    payload: Mapping[str, object],
    *,
    value_mode: str = "delta",
    ax=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload,
        ("plot_table", "metadata"),
        owner="plot_region_contrast_profile_matplotlib",
    )
    plot_table = _require_payload_table(payload, "plot_table")
    metadata = payload.get("metadata", {}) if isinstance(payload.get("metadata", {}), Mapping) else {}
    plot_table = _prepare_region_contrast_value_mode_table(plot_table, value_mode=value_mode)

    fig, ax = _make_axis(ax=ax, figsize=(8, 4))

    if plot_table.empty:
        ax.set_title(title or "Region contrast profile")
        ax.set_xlabel("position")
        ax.set_ylabel(value_mode.replace("_", " ").title())
        return fig, ax

    x_column = "plot_x" if "plot_x" in plot_table.columns else "position"
    grouped = (
        plot_table.loc[:, [x_column, "value"]]
        .groupby(x_column, as_index=False, sort=True)
        .mean(numeric_only=True)
    )
    ax.plot(grouped[x_column], grouped["value"], marker="o", linewidth=1.5)

    axis_table = payload.get("axis_table")
    if isinstance(axis_table, pd.DataFrame) and not axis_table.empty and {"axis_min", "axis_max"}.issubset(axis_table.columns):
        ax.set_xlim(axis_table["axis_min"].min(), axis_table["axis_max"].max())

    ax.set_title(title or "Region contrast profile")
    ax.set_xlabel("position" if x_column == "position" else x_column)
    ax.set_ylabel(value_mode.replace("_", " ").title())

    if metadata.get("orientation") == "region_5to3":
        ax.set_xlabel("oriented position")

    return fig, ax


def plot_region_contrast_heatmap_matplotlib(
    payload: Mapping[str, object],
    *,
    value_mode: str = "delta",
    ax=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload,
        ("plot_table", "metadata"),
        owner="plot_region_contrast_heatmap_matplotlib",
    )
    plot_table = _require_payload_table(payload, "plot_table")
    summary_table = payload.get("summary_table")
    if summary_table is not None and not isinstance(summary_table, pd.DataFrame):
        raise TypeError("plot payload key 'summary_table' must be a pandas DataFrame.")
    plot_table = _prepare_region_contrast_value_mode_table(plot_table, value_mode=value_mode)

    fig, ax = _make_axis(ax=ax, figsize=(8, 4))

    if plot_table.empty:
        ax.set_title(title or "Region contrast heatmap")
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

    if y_column == "row_order":
        if isinstance(summary_table, pd.DataFrame) and {"region_id", "row_order"}.issubset(summary_table.columns):
            summary_lookup = (
                summary_table.loc[:, ["region_id", "row_order"]]
                .drop_duplicates()
                .sort_values("row_order", kind="stable")
            )
        elif "region_id" in plot_table.columns:
            summary_lookup = (
                plot_table.loc[:, ["region_id", "row_order"]]
                .drop_duplicates()
                .sort_values("row_order", kind="stable")
            )
        else:
            summary_lookup = None
        if summary_lookup is not None:
            ax.set_yticks(range(len(summary_lookup)))
            ax.set_yticklabels(summary_lookup["region_id"].astype(str).tolist())
        else:
            ax.set_ylabel("region")
    else:
        ax.set_ylabel("region")

    ax.set_xticks(range(len(heatmap.columns)))
    ax.set_xticklabels([str(value) for value in heatmap.columns.tolist()], rotation=45, ha="right")
    ax.set_title(title or "Region contrast heatmap")
    ax.set_xlabel("position" if x_column == "position" else x_column)

    return fig, ax
