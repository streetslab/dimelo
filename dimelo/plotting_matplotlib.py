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


def _make_axes(*, axes=None, n_axes: int, figsize=(8, 4)):
    plt = _import_pyplot()
    n_axes = max(int(n_axes), 1)
    if axes is not None:
        if isinstance(axes, (list, tuple)):
            normalized_axes = list(axes)
        else:
            try:
                normalized_axes = list(axes.ravel())
            except AttributeError:
                normalized_axes = [axes]
        if not normalized_axes:
            raise ValueError("axes must contain at least one Matplotlib axis.")
        return normalized_axes[0].figure, normalized_axes

    fig, created_axes = plt.subplots(
        n_axes,
        1,
        figsize=(figsize[0], max(figsize[1], 2.5 * n_axes)),
        squeeze=False,
    )
    return fig, list(created_axes.ravel())


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


def plot_region_discovery_scan_matplotlib(
    payload: Mapping[str, object],
    *,
    axes=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload,
        ("scan_table", "hit_table", "metadata"),
        owner="plot_region_discovery_scan_matplotlib",
    )
    scan_table = _require_payload_table(payload, "scan_table")
    hit_table = _require_payload_table(payload, "hit_table")
    metadata = payload.get("metadata", {}) if isinstance(payload.get("metadata", {}), Mapping) else {}
    score_column = str(metadata.get("score_column") or "score_value")
    contigs = list(metadata.get("contig_order") or scan_table.get("contig", pd.Series(dtype="object")).dropna().unique())

    fig, axes = _make_axes(axes=axes, n_axes=len(contigs) or 1, figsize=(8, 3))

    if not contigs:
        axes[0].set_title(title or "Region discovery scan")
        axes[0].set_xlabel("Window midpoint")
        axes[0].set_ylabel(score_column.replace("_", " ").title())
        return fig, axes

    for index, contig in enumerate(contigs):
        ax = axes[index]
        contig_table = scan_table.loc[scan_table["contig"] == contig]
        if not contig_table.empty:
            ax.plot(
                contig_table["window_midpoint"],
                contig_table[score_column],
                marker="o",
                linewidth=1.5,
            )
        contig_hits = hit_table.loc[hit_table["contig"] == contig] if not hit_table.empty else hit_table
        if not contig_hits.empty and {"start", "end"}.issubset(contig_hits.columns):
            for _, hit in contig_hits.iterrows():
                ax.axvspan(float(hit["start"]), float(hit["end"]), color="tab:red", alpha=0.15)
        ax.set_title(str(contig))
        ax.set_xlabel("Window midpoint")
        ax.set_ylabel(score_column.replace("_", " ").title())

    for ax in axes[len(contigs):]:
        ax.set_visible(False)

    if title:
        fig.suptitle(title)
    return fig, axes


def plot_region_discovery_hit_context_matplotlib(
    payload: Mapping[str, object],
    *,
    axes=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload,
        ("context_table", "selected_hits", "metadata"),
        owner="plot_region_discovery_hit_context_matplotlib",
    )
    context_table = _require_payload_table(payload, "context_table")
    selected_hits = _require_payload_table(payload, "selected_hits")
    metadata = payload.get("metadata", {}) if isinstance(payload.get("metadata", {}), Mapping) else {}
    score_column = str(metadata.get("score_column") or "score_value")
    hit_ids = list(
        selected_hits.get("window_id", pd.Series(dtype="object")).dropna().astype(str).tolist()
    )

    fig, axes = _make_axes(axes=axes, n_axes=len(hit_ids) or 1, figsize=(8, 3))

    if not hit_ids:
        axes[0].set_title(title or "Region discovery hit context")
        axes[0].set_xlabel("Window midpoint")
        axes[0].set_ylabel(score_column.replace("_", " ").title())
        return fig, axes

    for index, hit_id in enumerate(hit_ids):
        ax = axes[index]
        hit_context = context_table.loc[context_table["selected_hit_id"].astype(str) == hit_id]
        if not hit_context.empty:
            ax.plot(
                hit_context["window_midpoint"],
                hit_context[score_column],
                marker="o",
                linewidth=1.5,
            )
            selected_rows = hit_context.loc[hit_context["is_selected_hit"]]
            if not selected_rows.empty:
                ax.scatter(
                    selected_rows["window_midpoint"],
                    selected_rows[score_column],
                    color="tab:red",
                    zorder=3,
                )
        ax.set_title(hit_id)
        ax.set_xlabel("Window midpoint")
        ax.set_ylabel(score_column.replace("_", " ").title())

    for ax in axes[len(hit_ids):]:
        ax.set_visible(False)

    if title:
        fig.suptitle(title)
    return fig, axes


def plot_shared_cluster_distribution_matplotlib(
    payload: Mapping[str, object],
    *,
    level: str = "condition",
    ax=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload,
        ("sample_distribution", "condition_distribution", "metadata"),
        owner="plot_shared_cluster_distribution_matplotlib",
    )
    if level not in {"sample", "condition"}:
        raise ValueError("level must be 'sample' or 'condition'.")

    table = _require_payload_table(
        payload,
        "sample_distribution" if level == "sample" else "condition_distribution",
    )
    fig, ax = _make_axis(ax=ax, figsize=(8, 4))

    if not table.empty:
        x_column = "sample_id" if level == "sample" else "condition"
        pivot = (
            table.pivot_table(
                index=x_column,
                columns="cluster",
                values="fraction",
                aggfunc="mean",
                fill_value=0.0,
            )
            .sort_index(axis=0)
            .sort_index(axis=1)
        )
        pivot.plot(kind="bar", stacked=True, ax=ax)
        ax.tick_params(axis="x", rotation=45)

    ax.set_xlabel("Sample" if level == "sample" else "Condition")
    ax.set_ylabel("Fraction")
    ax.set_title(title or "Shared cluster distribution")
    return fig, ax


def plot_shared_cluster_change_matplotlib(
    payload: Mapping[str, object],
    *,
    ax=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload,
        ("distribution_change", "metadata"),
        owner="plot_shared_cluster_change_matplotlib",
    )
    change_table = _require_payload_table(payload, "distribution_change")
    fig, ax = _make_axis(ax=ax, figsize=(8, 4))

    if not change_table.empty:
        matrix = (
            change_table.pivot_table(
                index="condition",
                columns="cluster",
                values="delta_fraction",
                aggfunc="mean",
                fill_value=0.0,
            )
            .sort_index(axis=0)
            .sort_index(axis=1)
        )
        max_abs = float(matrix.abs().to_numpy().max()) if not matrix.empty else 0.0
        if max_abs == 0.0:
            max_abs = 1.0
        image = ax.imshow(
            matrix.to_numpy(),
            aspect="auto",
            origin="lower",
            interpolation="nearest",
            cmap="coolwarm",
            vmin=-max_abs,
            vmax=max_abs,
        )
        ax.figure.colorbar(image, ax=ax, label="Delta Fraction")
        ax.set_xticks(range(len(matrix.columns)))
        ax.set_xticklabels([str(value) for value in matrix.columns.tolist()], rotation=45, ha="right")
        ax.set_yticks(range(len(matrix.index)))
        ax.set_yticklabels([str(value) for value in matrix.index.tolist()])

    ax.set_xlabel("Cluster")
    ax.set_ylabel("Condition")
    ax.set_title(title or "Shared cluster change")
    return fig, ax


def plot_global_analysis_summary_matplotlib(
    payload: Mapping[str, object],
    *,
    level: str = "condition",
    ax=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload,
        ("sample_summary", "condition_summary", "normalization_table", "metadata"),
        owner="plot_global_analysis_summary_matplotlib",
    )
    if level not in {"sample", "condition"}:
        raise ValueError("level must be 'sample' or 'condition'.")

    table = _require_payload_table(payload, "sample_summary" if level == "sample" else "condition_summary")
    fig, ax = _make_axis(ax=ax, figsize=(8, 4))

    if not table.empty:
        if level == "sample":
            x_values = table["sample_id"].astype(str)
            y_values = table["global_fraction"]
        else:
            x_values = table["condition"].astype(str)
            y_values = table["global_fraction_mean"]
        ax.bar(x_values, y_values)
        ax.tick_params(axis="x", rotation=45)

    ax.set_ylabel("Global Fraction")
    ax.set_xlabel("Sample" if level == "sample" else "Condition")
    ax.set_title(title or "Global analysis summary")
    return fig, ax


def plot_global_analysis_window_matplotlib(
    payload: Mapping[str, object],
    *,
    axes=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload,
        ("window_table", "metadata"),
        owner="plot_global_analysis_window_matplotlib",
    )
    window_table = _require_payload_table(payload, "window_table")
    metadata = payload.get("metadata", {}) if isinstance(payload.get("metadata", {}), Mapping) else {}
    contigs = list(metadata.get("contig_order") or window_table.get("contig", pd.Series(dtype="object")).dropna().unique())

    fig, axes = _make_axes(axes=axes, n_axes=len(contigs) or 1, figsize=(8, 3))

    if not contigs:
        axes[0].set_title(title or "Global analysis window")
        axes[0].set_xlabel("Window midpoint")
        axes[0].set_ylabel("Window Fraction")
        return fig, axes

    for index, contig in enumerate(contigs):
        ax = axes[index]
        contig_table = window_table.loc[window_table["contig"] == contig]
        if not contig_table.empty:
            grouped = (
                contig_table.loc[:, ["window_midpoint", "window_fraction"]]
                .groupby("window_midpoint", as_index=False, sort=True)
                .mean(numeric_only=True)
            )
            ax.plot(grouped["window_midpoint"], grouped["window_fraction"], marker="o", linewidth=1.5)
        ax.set_title(str(contig))
        ax.set_xlabel("Window midpoint")
        ax.set_ylabel("Window Fraction")

    for ax in axes[len(contigs):]:
        ax.set_visible(False)

    if title:
        fig.suptitle(title)
    return fig, axes
