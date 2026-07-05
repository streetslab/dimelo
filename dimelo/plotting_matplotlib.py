from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
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


def _require_payload_keys(
    payload: Mapping[str, object], keys: tuple[str, ...], *, owner: str
) -> None:
    missing = [key for key in keys if key not in payload]
    if missing:
        raise ValueError(
            f"{owner} payload is missing required keys: {', '.join(missing)}"
        )


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


def _ordered_cluster_labels(
    metadata: Mapping[str, object], observed_clusters: list[object]
) -> list[object]:
    cluster_labels = list(metadata.get("cluster_labels") or [])
    if not cluster_labels:
        return list(observed_clusters)

    observed_cluster_set = set(observed_clusters)
    ordered_clusters = [
        cluster for cluster in cluster_labels if cluster in observed_cluster_set
    ]
    ordered_clusters.extend(
        [cluster for cluster in observed_clusters if cluster not in cluster_labels]
    )
    return ordered_clusters


def _ordered_unique_values(table: pd.DataFrame, column: str) -> list[object]:
    return table.loc[:, column].drop_duplicates().tolist()


def _format_sort_value(value: object) -> str:
    if isinstance(value, (list, tuple)):
        return " > ".join(str(item) for item in value)
    return str(value)


def _read_cluster_association_title(
    metadata: Mapping[str, object],
    *,
    row_annotation_columns: Sequence[str],
    group_region_labels: bool | None,
    region_label_mode: str,
) -> str:
    value_mode = str(metadata.get("value_mode") or "fraction").replace("_", " ")
    region_sort = _format_sort_value(metadata.get("region_sort", "input"))
    cluster_sort = str(metadata.get("cluster_sort", "input"))
    parts = [
        f"Read-cluster association ({value_mode})",
        f"regions: {region_sort}",
        f"clusters: {cluster_sort}",
    ]
    metadata_region_sort = metadata.get("region_sort")
    if metadata_region_sort == "association_strength" or (
        isinstance(metadata_region_sort, (list, tuple))
        and "association_strength" in metadata_region_sort
    ):
        parts.append(
            f"strength: {metadata.get('association_strength_aggregate', 'max')}"
        )
    if row_annotation_columns:
        parts.append(f"annotated by {', '.join(row_annotation_columns)}")
    if group_region_labels:
        parts.append("grouped labels")
    elif region_label_mode != "auto":
        parts.append(f"labels: {region_label_mode}")
    return " | ".join(parts)


def _normalize_row_annotation_columns(
    *,
    row_annotation_column: str | None,
    row_annotation_columns: Sequence[str] | None,
) -> list[str]:
    columns: list[str] = []
    if row_annotation_columns is not None:
        if isinstance(row_annotation_columns, str):
            raise TypeError(
                "row_annotation_columns must be a sequence of column names, not a string."
            )
        columns.extend(str(column) for column in row_annotation_columns)
    if row_annotation_column is not None:
        legacy_column = str(row_annotation_column)
        if legacy_column not in columns:
            columns.insert(0, legacy_column)
    return list(dict.fromkeys(columns))


def _row_annotation_title(
    column: str,
    *,
    row_annotation_title: str | None,
    row_annotation_titles: Mapping[str, str] | None,
    n_columns: int,
) -> str:
    if row_annotation_titles is not None and column in row_annotation_titles:
        return str(row_annotation_titles[column])
    if n_columns == 1 and row_annotation_title is not None:
        return row_annotation_title
    return column


def _row_annotation_palette(
    column: str,
    *,
    row_annotation_palette: Mapping[str, str] | None,
    row_annotation_palettes: Mapping[str, Mapping[str, str]] | None,
    n_columns: int,
) -> Mapping[str, str] | None:
    if row_annotation_palettes is not None and column in row_annotation_palettes:
        return row_annotation_palettes[column]
    if n_columns == 1:
        return row_annotation_palette
    return None


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
        raise ValueError(
            f"plot payload does not contain rows for value_mode={value_mode!r}."
        )
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
    metadata = (
        payload.get("metadata", {})
        if isinstance(payload.get("metadata", {}), Mapping)
        else {}
    )
    plot_table = _prepare_region_contrast_value_mode_table(
        plot_table, value_mode=value_mode
    )

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
    if (
        isinstance(axis_table, pd.DataFrame)
        and not axis_table.empty
        and {"axis_min", "axis_max"}.issubset(axis_table.columns)
    ):
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
    plot_table = _prepare_region_contrast_value_mode_table(
        plot_table, value_mode=value_mode
    )

    fig, ax = _make_axis(ax=ax, figsize=(8, 4))

    if plot_table.empty:
        ax.set_title(title or "Region contrast heatmap")
        ax.set_xlabel("position")
        ax.set_ylabel("region")
        return fig, ax

    x_column = "plot_x" if "plot_x" in plot_table.columns else "position"
    y_column = "row_order" if "row_order" in plot_table.columns else "region_id"

    heatmap_table = plot_table.loc[:, [x_column, y_column, "value"]].copy()
    heatmap = (
        heatmap_table.pivot_table(
            index=y_column,
            columns=x_column,
            values="value",
            aggfunc="mean",
        )
        .sort_index(axis=0)
        .sort_index(axis=1)
    )

    image = ax.imshow(
        heatmap.to_numpy(),
        aspect="auto",
        # Review fix #10: use origin="upper" to match every sibling heatmap in
        # this module; origin="lower" rendered the region-contrast panel flipped
        # vertically relative to its row-order/summary annotations.
        origin="upper",
        interpolation="nearest",
    )
    ax.figure.colorbar(image, ax=ax, label="value")

    if y_column == "row_order":
        if isinstance(summary_table, pd.DataFrame) and {
            "region_id",
            "row_order",
        }.issubset(summary_table.columns):
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
    ax.set_xticklabels(
        [str(value) for value in heatmap.columns.tolist()], rotation=45, ha="right"
    )
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
    metadata = (
        payload.get("metadata", {})
        if isinstance(payload.get("metadata", {}), Mapping)
        else {}
    )
    score_column = str(metadata.get("score_column") or "score_value")
    contigs = list(
        metadata.get("contig_order")
        or scan_table.get("contig", pd.Series(dtype="object")).dropna().unique()
    )

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
        contig_hits = (
            hit_table.loc[hit_table["contig"] == contig]
            if not hit_table.empty
            else hit_table
        )
        if not contig_hits.empty and {"start", "end"}.issubset(contig_hits.columns):
            for _, hit in contig_hits.iterrows():
                ax.axvspan(
                    float(hit["start"]), float(hit["end"]), color="tab:red", alpha=0.15
                )
        ax.set_title(str(contig))
        ax.set_xlabel("Window midpoint")
        ax.set_ylabel(score_column.replace("_", " ").title())

    for ax in axes[len(contigs) :]:
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
    metadata = (
        payload.get("metadata", {})
        if isinstance(payload.get("metadata", {}), Mapping)
        else {}
    )
    score_column = str(metadata.get("score_column") or "score_value")
    hit_ids = list(
        selected_hits.get("window_id", pd.Series(dtype="object"))
        .dropna()
        .astype(str)
        .tolist()
    )

    fig, axes = _make_axes(axes=axes, n_axes=len(hit_ids) or 1, figsize=(8, 3))

    if not hit_ids:
        axes[0].set_title(title or "Region discovery hit context")
        axes[0].set_xlabel("Window midpoint")
        axes[0].set_ylabel(score_column.replace("_", " ").title())
        return fig, axes

    for index, hit_id in enumerate(hit_ids):
        ax = axes[index]
        hit_context = context_table.loc[
            context_table["selected_hit_id"].astype(str) == hit_id
        ]
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

    for ax in axes[len(hit_ids) :]:
        ax.set_visible(False)

    if title:
        fig.suptitle(title)
    return fig, axes


def plot_joint_state_distribution_matplotlib(
    payload: Mapping[str, object],
    *,
    ax=None,
    title: str | None = None,
):
    """Stacked bar of per-pair joint-state fractions (both / a_only / b_only / neither)."""
    _require_payload_keys(
        payload,
        ("distribution_table", "metadata"),
        owner="plot_joint_state_distribution_matplotlib",
    )
    table = _require_payload_table(payload, "distribution_table")
    metadata = (
        payload.get("metadata", {})
        if isinstance(payload.get("metadata", {}), Mapping)
        else {}
    )
    states = list(metadata.get("states") or ["both", "a_only", "b_only", "neither"])
    fig, ax = _make_axis(ax=ax, figsize=(max(6.0, len(table) * 0.8), 4))
    if table.empty:
        ax.set_title(title or "Joint occupancy states")
        ax.set_ylabel("fraction of spanning reads")
        return fig, ax
    positions = np.arange(len(table))
    colors = {
        "both": "tab:red",
        "a_only": "tab:orange",
        "b_only": "tab:blue",
        "neither": "tab:gray",
    }
    bottom = np.zeros(len(table))
    for state in states:
        heights = table[f"frac_{state}"].astype(float).fillna(0.0).to_numpy()
        ax.bar(positions, heights, bottom=bottom, label=state, color=colors.get(state))
        bottom += heights
    ax.set_xticks(positions)
    ax.set_xticklabels(table["pair_id"].astype(str).tolist(), rotation=45, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("fraction of spanning reads")
    ax.legend(frameon=False, ncol=len(states), fontsize=8)
    ax.set_title(title or "Joint occupancy states")
    return fig, ax


def plot_joint_occupancy_distance_matplotlib(
    payload: Mapping[str, object],
    *,
    ax=None,
    title: str | None = None,
):
    """Scatter of log2 observed/expected co-occupancy vs 1D anchor distance; significant
    pairs highlighted. Positive coupling at short distance flags a possible trans-contact
    artifact (Q4)."""
    _require_payload_keys(
        payload,
        ("distance_table", "metadata"),
        owner="plot_joint_occupancy_distance_matplotlib",
    )
    table = _require_payload_table(payload, "distance_table")
    fig, ax = _make_axis(ax=ax, figsize=(7, 4))
    if table.empty:
        ax.set_title(title or "Joint occupancy vs anchor distance")
        ax.set_xlabel("anchor distance (bp)")
        ax.set_ylabel("log2 observed / expected co-occupancy")
        return fig, ax
    significant = table["significant"].astype(bool).to_numpy()
    ax.scatter(
        table.loc[~significant, "distance_bp"], table.loc[~significant, "log2_obs_exp"],
        color="tab:gray", alpha=0.7, label="n.s.",
    )
    ax.scatter(
        table.loc[significant, "distance_bp"], table.loc[significant, "log2_obs_exp"],
        color="tab:red", label="significant",
    )
    ax.axhline(0.0, color="black", linewidth=0.6)
    ax.set_xlabel("anchor distance (bp)")
    ax.set_ylabel("log2 observed / expected co-occupancy")
    ax.legend(frameon=False)
    ax.set_title(title or "Joint occupancy vs anchor distance")
    return fig, ax


def plot_footprint_profile_matplotlib(
    payload: Mapping[str, object],
    *,
    ax=None,
    title: str | None = None,
):
    """Aggregate footprint profile: methylation density (left axis) and protected-state
    posterior (right axis) vs position; the footprint is a methylation dip / posterior
    rise at the anchor."""
    _require_payload_keys(
        payload, ("profile_table", "metadata"), owner="plot_footprint_profile_matplotlib"
    )
    table = _require_payload_table(payload, "profile_table")
    metadata = (
        payload.get("metadata", {})
        if isinstance(payload.get("metadata", {}), Mapping)
        else {}
    )
    fig, ax = _make_axis(ax=ax, figsize=(8, 4))
    if table.empty:
        ax.set_title(title or "Aggregate footprint profile")
        ax.set_xlabel("position")
        return fig, ax
    ax.plot(
        table["position"], table["mean_methylation"],
        color="tab:blue", marker=".", label="methylation",
    )
    ax.set_xlabel("position")
    ax.set_ylabel("mean methylation", color="tab:blue")
    ax.set_ylim(0, 1)
    ax.tick_params(axis="y", labelcolor="tab:blue")

    ax_right = ax.twinx()
    ax_right.plot(
        table["position"], table["mean_protected_posterior"],
        color="tab:red", marker=".", linestyle="--", label="protected posterior",
    )
    ax_right.set_ylabel("mean protected posterior", color="tab:red")
    ax_right.set_ylim(0, 1)
    ax_right.tick_params(axis="y", labelcolor="tab:red")

    anchor_index = metadata.get("anchor_index")
    if anchor_index is not None:
        ax.axvline(float(anchor_index), color="black", linewidth=0.8, linestyle=":")
    ax.set_title(title or "Aggregate footprint profile")
    return fig, ax


def plot_binding_strength_matplotlib(
    payload: Mapping[str, object],
    *,
    ax=None,
    title: str | None = None,
):
    """Ranked per-site binding strength: target occupancy (Beta-posterior mean + credible
    interval) with control occupancy overlaid; significant sites highlighted."""
    _require_payload_keys(
        payload, ("binding_table", "metadata"), owner="plot_binding_strength_matplotlib"
    )
    table = _require_payload_table(payload, "binding_table")
    fig, ax = _make_axis(ax=ax, figsize=(max(6.0, len(table) * 0.7), 4))
    if table.empty:
        ax.set_title(title or "Binding strength (control-calibrated occupancy)")
        ax.set_ylabel("occupancy")
        return fig, ax

    positions = np.arange(len(table))
    heights = table["occupancy_posterior_mean"].astype(float).to_numpy()
    lower = np.clip(heights - table["occupancy_ci_low"].astype(float).to_numpy(), 0, None)
    upper = np.clip(table["occupancy_ci_high"].astype(float).to_numpy() - heights, 0, None)
    significant = table["significant"].astype(bool).to_numpy()
    colors = np.where(significant, "tab:red", "tab:gray")
    ax.bar(
        positions, heights, color=colors,
        yerr=np.vstack([lower, upper]), capsize=3, ecolor="black",
    )
    # control occupancy overlaid as markers
    ax.scatter(
        positions, table["control_occupancy"].astype(float),
        color="black", marker="_", s=120, zorder=3, label="control occupancy",
    )
    ax.set_xticks(positions)
    ax.set_xticklabels(table["region_id"].astype(str).tolist(), rotation=45, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("occupancy (posterior mean ± CI)")
    ax.legend(frameon=False, loc="upper right")
    ax.set_title(title or "Binding strength (control-calibrated occupancy)")
    return fig, ax


def plot_true_signal_read_raster_matplotlib(
    payload: Mapping[str, object],
    *,
    axes=None,
    title: str | None = None,
):
    """Single-molecule raster: reads (rows) vs read_mod_fraction, colored by call."""
    _require_payload_keys(
        payload, ("read_table", "metadata"), owner="plot_true_signal_read_raster_matplotlib"
    )
    read_table = _require_payload_table(payload, "read_table")
    metadata = (
        payload.get("metadata", {})
        if isinstance(payload.get("metadata", {}), Mapping)
        else {}
    )
    color_by = str(metadata.get("color_by") or "is_true_signal")
    region_ids = list(
        metadata.get("region_ids")
        or read_table.get("region_id", pd.Series(dtype="object")).drop_duplicates()
    )

    fig, axes = _make_axes(axes=axes, n_axes=len(region_ids) or 1, figsize=(7, 3))
    if not region_ids:
        axes[0].set_title(title or "Single-molecule true-signal raster")
        axes[0].set_xlabel("read modification fraction")
        axes[0].set_ylabel("read")
        return fig, axes

    for index, region in enumerate(region_ids):
        ax = axes[index]
        sub = read_table.loc[read_table["region_id"] == region]
        if not sub.empty:
            if color_by == "occupancy_posterior":
                sc = ax.scatter(
                    sub["read_mod_fraction"], sub["read_order"],
                    c=sub["occupancy_posterior"], cmap="viridis", vmin=0.0, vmax=1.0, s=18,
                )
                fig.colorbar(sc, ax=ax, label="occupancy posterior")
            else:
                is_sig = sub["is_true_signal"].astype(bool)
                ax.scatter(
                    sub.loc[~is_sig, "read_mod_fraction"], sub.loc[~is_sig, "read_order"],
                    color="tab:gray", s=18, alpha=0.7, label="background",
                )
                ax.scatter(
                    sub.loc[is_sig, "read_mod_fraction"], sub.loc[is_sig, "read_order"],
                    color="tab:red", s=18, label="true signal",
                )
                ax.legend(frameon=False, loc="lower right", fontsize=8)
            background_rate = sub["background_rate"].iloc[0]
            ax.axvline(float(background_rate), color="black", linewidth=0.8, linestyle="--")
        ax.set_xlim(-0.02, 1.02)
        ax.set_title(str(region))
        ax.set_xlabel("read modification fraction")
        ax.set_ylabel("read")

    for ax in axes[len(region_ids) :]:
        ax.set_visible(False)
    if title:
        fig.suptitle(title)
    return fig, axes


def plot_occupancy_rate_track_matplotlib(
    payload: Mapping[str, object],
    *,
    ax=None,
    title: str | None = None,
):
    """Bar of per-site control-calibrated occupancy rate."""
    _require_payload_keys(
        payload, ("track_table", "metadata"), owner="plot_occupancy_rate_track_matplotlib"
    )
    track = _require_payload_table(payload, "track_table")
    fig, ax = _make_axis(ax=ax, figsize=(max(6.0, len(track) * 0.6), 4))
    if track.empty:
        ax.set_title(title or "Control-calibrated occupancy rate")
        ax.set_ylabel("occupancy rate")
        return fig, ax
    positions = np.arange(len(track))
    has_posterior = {
        "occupancy_posterior_mean",
        "occupancy_ci_low",
        "occupancy_ci_high",
    }.issubset(track.columns)
    if has_posterior:
        # Bayesian Beta-posterior mean with a credible interval
        heights = track["occupancy_posterior_mean"].astype(float).to_numpy()
        lower = heights - track["occupancy_ci_low"].astype(float).to_numpy()
        upper = track["occupancy_ci_high"].astype(float).to_numpy() - heights
        ax.bar(
            positions,
            heights,
            color="tab:red",
            yerr=np.vstack([np.clip(lower, 0, None), np.clip(upper, 0, None)]),
            capsize=3,
            ecolor="black",
        )
        ylabel = "occupancy (posterior mean ± CI)"
    else:
        ax.bar(positions, track["occupancy_rate"].astype(float), color="tab:red")
        ylabel = "occupancy rate"
    ax.set_xticks(positions)
    ax.set_xticklabels(track["region_id"].astype(str).tolist(), rotation=45, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel(ylabel)
    ax.set_title(title or "Control-calibrated occupancy rate")
    return fig, ax


def plot_background_removed_pileup_overlay_matplotlib(
    payload: Mapping[str, object],
    *,
    ax=None,
    title: str | None = None,
):
    """Grouped bars of raw vs background-removed modification fraction per site."""
    _require_payload_keys(
        payload,
        ("overlay_table", "metadata"),
        owner="plot_background_removed_pileup_overlay_matplotlib",
    )
    overlay = _require_payload_table(payload, "overlay_table")
    metadata = (
        payload.get("metadata", {})
        if isinstance(payload.get("metadata", {}), Mapping)
        else {}
    )
    region_order = list(
        metadata.get("region_order")
        or overlay.get("region_id", pd.Series(dtype="object")).drop_duplicates()
    )
    fig, ax = _make_axis(ax=ax, figsize=(max(6.0, len(region_order) * 0.7), 4))
    if overlay.empty or not region_order:
        ax.set_title(title or "Raw vs background-removed pileup")
        ax.set_ylabel("modification fraction")
        return fig, ax
    raw = overlay.loc[overlay["stage"] == "raw"].set_index("region_id")[
        "mod_fraction"
    ].to_dict()
    removed = overlay.loc[overlay["stage"] == "background_removed"].set_index(
        "region_id"
    )["mod_fraction"].to_dict()
    positions = np.arange(len(region_order))
    width = 0.4
    ax.bar(positions - width / 2, [float(raw.get(r, 0.0)) for r in region_order],
           width, label="raw", color="tab:gray")
    ax.bar(positions + width / 2, [float(removed.get(r, 0.0)) for r in region_order],
           width, label="background-removed", color="tab:red")
    ax.set_xticks(positions)
    ax.set_xticklabels([str(r) for r in region_order], rotation=45, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("modification fraction")
    ax.legend(frameon=False)
    ax.set_title(title or "Raw vs background-removed pileup")
    return fig, ax


def plot_region_contrast_correction_overlay_matplotlib(
    payload: Mapping[str, object],
    *,
    ax=None,
    title: str | None = None,
):
    """Grouped bars of raw vs background-corrected ``delta_fraction`` per region."""
    _require_payload_keys(
        payload,
        ("overlay_table", "metadata"),
        owner="plot_region_contrast_correction_overlay_matplotlib",
    )
    overlay = _require_payload_table(payload, "overlay_table")
    metadata = (
        payload.get("metadata", {})
        if isinstance(payload.get("metadata", {}), Mapping)
        else {}
    )
    region_order = list(
        metadata.get("region_order")
        or overlay.get("region_id", pd.Series(dtype="object")).drop_duplicates()
    )

    fig, ax = _make_axis(ax=ax, figsize=(max(6.0, len(region_order) * 0.7), 4))

    if overlay.empty or not region_order:
        ax.set_title(title or "Background correction: raw vs corrected")
        ax.set_xlabel("Region")
        ax.set_ylabel("delta fraction (numerator - denominator)")
        return fig, ax

    raw = (
        overlay.loc[overlay["stage"] == "raw"]
        .set_index("region_id")["delta_fraction"]
        .to_dict()
    )
    corrected = (
        overlay.loc[overlay["stage"] == "corrected"]
        .set_index("region_id")["delta_fraction"]
        .to_dict()
    )
    positions = np.arange(len(region_order))
    width = 0.4
    ax.bar(
        positions - width / 2,
        [float(raw.get(region, 0.0)) for region in region_order],
        width,
        label="raw",
        color="tab:gray",
    )
    ax.bar(
        positions + width / 2,
        [float(corrected.get(region, 0.0)) for region in region_order],
        width,
        label="corrected",
        color="tab:blue",
    )
    ax.axhline(0.0, color="black", linewidth=0.6)
    ax.set_xticks(positions)
    ax.set_xticklabels([str(region) for region in region_order], rotation=45, ha="right")
    ax.set_ylabel("delta fraction (numerator - denominator)")
    ax.set_xlabel("Region")
    ax.legend(frameon=False)
    ax.set_title(title or "Background correction: raw vs corrected")
    return fig, ax


def plot_dmr_site_matplotlib(
    payload: Mapping[str, object],
    *,
    axes=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload,
        ("site_table", "metadata"),
        owner="plot_dmr_site_matplotlib",
    )
    site_table = _require_payload_table(payload, "site_table")
    metadata = (
        payload.get("metadata", {})
        if isinstance(payload.get("metadata", {}), Mapping)
        else {}
    )
    effect_column = str(metadata.get("effect_size_column") or "effect_size")
    y_column = effect_column if effect_column in site_table.columns else "neg_log10_pvalue"
    contigs = list(
        metadata.get("contig_order")
        or site_table.get("contig", pd.Series(dtype="object")).dropna().unique()
    )

    fig, axes = _make_axes(axes=axes, n_axes=len(contigs) or 1, figsize=(8, 3))

    y_label = y_column.replace("_", " ").title()
    if not contigs:
        axes[0].set_title(title or "DMR sites")
        axes[0].set_xlabel("Genomic midpoint")
        axes[0].set_ylabel(y_label)
        return fig, axes

    for index, contig in enumerate(contigs):
        ax = axes[index]
        contig_table = site_table.loc[site_table["contig"] == contig]
        if not contig_table.empty and y_column in contig_table.columns:
            is_hc = (
                contig_table["is_high_confidence"].astype(bool)
                if "is_high_confidence" in contig_table.columns
                else pd.Series(False, index=contig_table.index)
            )
            background = contig_table.loc[~is_hc]
            highlighted = contig_table.loc[is_hc]
            if not background.empty:
                ax.scatter(
                    background["midpoint"],
                    background[y_column],
                    s=12,
                    color="tab:gray",
                    alpha=0.6,
                )
            if not highlighted.empty:
                ax.scatter(
                    highlighted["midpoint"],
                    highlighted[y_column],
                    s=20,
                    color="tab:red",
                    zorder=3,
                    label="high confidence",
                )
        ax.axhline(0.0, color="black", linewidth=0.6, alpha=0.4)
        ax.set_title(str(contig))
        ax.set_xlabel("Genomic midpoint")
        ax.set_ylabel(y_label)

    for ax in axes[len(contigs) :]:
        ax.set_visible(False)

    if title:
        fig.suptitle(title)
    return fig, axes


def plot_dmr_segment_matplotlib(
    payload: Mapping[str, object],
    *,
    axes=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload,
        ("segment_table", "metadata"),
        owner="plot_dmr_segment_matplotlib",
    )
    segment_table = _require_payload_table(payload, "segment_table")
    metadata = (
        payload.get("metadata", {})
        if isinstance(payload.get("metadata", {}), Mapping)
        else {}
    )
    score_column = str(metadata.get("score_column") or "score")
    contigs = list(
        metadata.get("contig_order")
        or segment_table.get("contig", pd.Series(dtype="object")).dropna().unique()
    )

    fig, axes = _make_axes(axes=axes, n_axes=len(contigs) or 1, figsize=(8, 3))

    y_label = score_column.replace("_", " ").title()
    if not contigs:
        axes[0].set_title(title or "DMR segments")
        axes[0].set_xlabel("Genomic position")
        axes[0].set_ylabel(y_label)
        return fig, axes

    for index, contig in enumerate(contigs):
        ax = axes[index]
        contig_table = segment_table.loc[segment_table["contig"] == contig]
        has_score = score_column in contig_table.columns
        for _, segment in contig_table.iterrows():
            height = float(segment[score_column]) if has_score else 1.0
            ax.plot(
                [float(segment["start"]), float(segment["end"])],
                [height, height],
                color="tab:blue",
                linewidth=3.0,
                solid_capstyle="butt",
            )
        ax.set_title(str(contig))
        ax.set_xlabel("Genomic position")
        ax.set_ylabel(y_label)

    for ax in axes[len(contigs) :]:
        ax.set_visible(False)

    if title:
        fig.suptitle(title)
    return fig, axes


def plot_dmr_multi_summary_matplotlib(
    payload: Mapping[str, object],
    *,
    ax=None,
    title: str | None = None,
    value_column: str = "n_significant",
):
    _require_payload_keys(
        payload,
        ("summary_table", "metadata"),
        owner="plot_dmr_multi_summary_matplotlib",
    )
    summary_table = _require_payload_table(payload, "summary_table")
    fig, ax = _make_axis(ax=ax, figsize=(8, 4))

    if summary_table.empty or value_column not in summary_table.columns:
        ax.set_title(title or "DMR multi-sample summary")
        ax.set_xlabel("Sample pair")
        ax.set_ylabel(value_column.replace("_", " ").title())
        return fig, ax

    positions = range(len(summary_table))
    ax.bar(list(positions), summary_table[value_column].astype(float), color="tab:blue")
    ax.set_xticks(list(positions))
    ax.set_xticklabels(
        summary_table["pair_name"].astype(str).tolist(),
        rotation=45,
        ha="right",
    )
    ax.set_xlabel("Sample pair")
    ax.set_ylabel(value_column.replace("_", " ").title())
    ax.set_title(title or "DMR multi-sample summary")
    return fig, ax


def plot_shared_cluster_distribution_matplotlib(
    payload: Mapping[str, object],
    *,
    level: str = "condition",
    ax=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload, ("metadata",), owner="plot_shared_cluster_distribution_matplotlib"
    )
    if level not in {"sample", "condition"}:
        raise ValueError("level must be 'sample' or 'condition'.")

    table = _require_payload_table(
        payload,
        "sample_distribution" if level == "sample" else "condition_distribution",
    )
    metadata = payload.get("metadata")
    if not isinstance(metadata, Mapping):
        raise TypeError("plot payload key 'metadata' must be a mapping.")

    fig, ax = _make_axis(ax=ax, figsize=(8, 4))

    if not table.empty:
        x_column = "sample_id" if level == "sample" else "condition"
        value_table = table.loc[:, [x_column, "cluster", "fraction"]].copy()
        if value_table.duplicated([x_column, "cluster"]).any():
            raise ValueError(
                "plot payload contains duplicate cluster fractions for the same x value."
            )
        x_order = _ordered_unique_values(value_table, x_column)
        cluster_order = _ordered_cluster_labels(
            metadata, _ordered_unique_values(value_table, "cluster")
        )
        pivot = value_table.set_index([x_column, "cluster"])["fraction"].unstack(
            "cluster", fill_value=0.0
        )
        pivot = pivot.reindex(index=x_order, columns=cluster_order, fill_value=0.0)
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
        payload, ("metadata",), owner="plot_shared_cluster_change_matplotlib"
    )
    change_table = _require_payload_table(payload, "distribution_change")
    metadata = payload.get("metadata")
    if not isinstance(metadata, Mapping):
        raise TypeError("plot payload key 'metadata' must be a mapping.")
    fig, ax = _make_axis(ax=ax, figsize=(8, 4))

    if not change_table.empty:
        value_table = change_table.loc[
            :, ["condition", "cluster", "delta_fraction"]
        ].copy()
        if value_table.duplicated(["condition", "cluster"]).any():
            raise ValueError(
                "plot payload contains duplicate delta fractions for the same condition."
            )
        condition_order = list(metadata.get("change_condition_order") or [])
        if not condition_order:
            condition_order = _ordered_unique_values(value_table, "condition")
        observed_clusters = _ordered_unique_values(value_table, "cluster")
        cluster_order = list(metadata.get("cluster_labels") or [])
        if cluster_order:
            cluster_order.extend(
                [
                    cluster
                    for cluster in observed_clusters
                    if cluster not in cluster_order
                ]
            )
        else:
            cluster_order = observed_clusters
        matrix = value_table.set_index(["condition", "cluster"])[
            "delta_fraction"
        ].unstack("cluster")
        matrix = matrix.reindex(index=condition_order, columns=cluster_order)
        finite_values = pd.Series(matrix.abs().to_numpy().ravel()).dropna()
        max_abs = float(finite_values.max()) if not finite_values.empty else 0.0
        if max_abs == 0.0:
            max_abs = 1.0
        image = ax.imshow(
            matrix.to_numpy(),
            aspect="auto",
            origin="upper",
            interpolation="nearest",
            cmap="coolwarm",
            vmin=-max_abs,
            vmax=max_abs,
        )
        ax.figure.colorbar(image, ax=ax, label="Delta Fraction")
        ax.set_xticks(range(len(matrix.columns)))
        ax.set_xticklabels(
            [str(value) for value in matrix.columns.tolist()], rotation=45, ha="right"
        )
        ax.set_yticks(range(len(matrix.index)))
        ax.set_yticklabels([str(value) for value in matrix.index.tolist()])

    ax.set_xlabel("Cluster")
    ax.set_ylabel("Condition")
    ax.set_title(title or "Shared cluster change")
    return fig, ax


def plot_shared_cluster_profile_heatmap_matplotlib(
    payload: Mapping[str, object],
    *,
    ax=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload,
        ("profile_table", "metadata"),
        owner="plot_shared_cluster_profile_heatmap_matplotlib",
    )
    profile_table = _require_payload_table(payload, "profile_table")
    metadata = payload.get("metadata")
    if not isinstance(metadata, Mapping):
        raise TypeError("plot payload key 'metadata' must be a mapping.")
    fig, ax = _make_axis(ax=ax, figsize=(8, 4))

    value_table = profile_table.loc[:, ["cluster", "feature", "value"]].copy()
    if value_table.duplicated(["cluster", "feature"]).any():
        raise ValueError(
            "plot payload contains duplicate profile values for the same cluster and feature."
        )

    feature_order = list(metadata.get("feature_names") or [])
    observed_features = _ordered_unique_values(value_table, "feature")
    if feature_order:
        feature_order.extend(
            [feature for feature in observed_features if feature not in feature_order]
        )
    else:
        feature_order = observed_features

    observed_clusters = _ordered_unique_values(value_table, "cluster")
    cluster_order = list(metadata.get("cluster_labels") or [])
    if cluster_order:
        cluster_order.extend(
            [cluster for cluster in observed_clusters if cluster not in cluster_order]
        )
    else:
        cluster_order = observed_clusters

    if feature_order:
        ax.set_xticks(range(len(feature_order)))
        ax.set_xticklabels(
            [str(value) for value in feature_order], rotation=45, ha="right"
        )
    if cluster_order:
        ax.set_yticks(range(len(cluster_order)))
        ax.set_yticklabels([str(value) for value in cluster_order])

    if feature_order and cluster_order:
        if value_table.empty:
            matrix = pd.DataFrame(
                index=cluster_order, columns=feature_order, dtype=float
            )
        else:
            matrix = value_table.set_index(["cluster", "feature"])["value"].unstack(
                "feature"
            )
            matrix = matrix.reindex(index=cluster_order, columns=feature_order)
        if matrix.notna().any().any():
            image = ax.imshow(
                matrix.to_numpy(),
                aspect="auto",
                origin="upper",
                interpolation="nearest",
            )
            ax.figure.colorbar(image, ax=ax, label="Value")
            ax.set_xticks(range(len(matrix.columns)))
            ax.set_xticklabels(
                [str(value) for value in matrix.columns.tolist()],
                rotation=45,
                ha="right",
            )
            ax.set_yticks(range(len(matrix.index)))
            ax.set_yticklabels([str(value) for value in matrix.index.tolist()])

    ax.set_xlabel("Feature")
    ax.set_ylabel("Cluster")
    ax.set_title(title or "Shared cluster profiles")
    return fig, ax


def plot_shared_cluster_profile_series_matplotlib(
    payload: Mapping[str, object],
    *,
    ax=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload,
        ("profile_table", "metadata"),
        owner="plot_shared_cluster_profile_series_matplotlib",
    )
    profile_table = _require_payload_table(payload, "profile_table")
    metadata = payload.get("metadata")
    if not isinstance(metadata, Mapping):
        raise TypeError("plot payload key 'metadata' must be a mapping.")
    fig, ax = _make_axis(ax=ax, figsize=(8, 4))

    value_table = profile_table.loc[:, ["cluster", "feature", "value"]].copy()
    if value_table.duplicated(["cluster", "feature"]).any():
        raise ValueError(
            "plot payload contains duplicate profile values for the same cluster and feature."
        )

    feature_order = list(metadata.get("feature_names") or [])
    observed_features = _ordered_unique_values(value_table, "feature")
    if feature_order:
        feature_order.extend(
            [feature for feature in observed_features if feature not in feature_order]
        )
    else:
        feature_order = observed_features
    feature_lookup = {feature: index for index, feature in enumerate(feature_order)}

    observed_clusters = _ordered_unique_values(value_table, "cluster")
    cluster_order = list(metadata.get("cluster_labels") or [])
    if cluster_order:
        cluster_order.extend(
            [cluster for cluster in observed_clusters if cluster not in cluster_order]
        )
    else:
        cluster_order = observed_clusters

    for cluster in cluster_order:
        cluster_table = value_table.loc[
            value_table["cluster"] == cluster, ["feature", "value"]
        ].copy()
        if feature_order:
            cluster_table["feature_order"] = cluster_table["feature"].map(
                feature_lookup
            )
            cluster_table = cluster_table.sort_values("feature_order", kind="stable")
            cluster_series = cluster_table.set_index("feature_order")["value"].reindex(
                range(len(feature_order))
            )
            cluster_series = cluster_series.astype(float)
            x_values = cluster_series.index.tolist()
            y_values = cluster_series.to_numpy()
        else:
            x_values = []
            y_values = []
        ax.plot(
            x_values,
            y_values,
            marker="o",
            linewidth=1.5,
            label=str(cluster),
        )

    if feature_order:
        ax.set_xticks(range(len(feature_order)))
        ax.set_xticklabels(feature_order, rotation=45, ha="right")
    if ax.lines:
        ax.legend(title="Cluster")

    ax.set_xlabel("Feature")
    ax.set_ylabel("Value")
    ax.set_title(title or "Shared cluster profile series")
    return fig, ax


def plot_shared_cluster_region_matplotlib(
    payload: Mapping[str, object],
    *,
    level: str = "condition",
    ax=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload,
        ("region_table", "condition_region_table", "metadata"),
        owner="plot_shared_cluster_region_matplotlib",
    )
    if level not in {"sample", "condition"}:
        raise ValueError("level must be 'sample' or 'condition'.")

    table = _require_payload_table(
        payload,
        "condition_region_table" if level == "condition" else "region_table",
    )
    metadata = payload.get("metadata")
    if not isinstance(metadata, Mapping):
        raise TypeError("plot payload key 'metadata' must be a mapping.")
    fig, ax = _make_axis(ax=ax, figsize=(8, 4))

    value_column = "fraction_mean" if level == "condition" else "fraction"
    label_column = "condition" if level == "condition" else "sample_id"
    value_table = table.loc[
        :, ["region_id", label_column, "cluster", value_column]
    ].copy()
    value_table["row_id"] = (
        value_table.loc[:, ["region_id", label_column]]
        .astype(str)
        .agg(" | ".join, axis=1)
    )
    if value_table.duplicated(["row_id", "cluster"]).any():
        raise ValueError(
            "plot payload contains duplicate region occupancy values for the same row and cluster."
        )

    row_order = _ordered_unique_values(value_table, "row_id")
    observed_clusters = _ordered_unique_values(value_table, "cluster")
    cluster_order = list(metadata.get("cluster_labels") or [])
    if cluster_order:
        cluster_order.extend(
            [cluster for cluster in observed_clusters if cluster not in cluster_order]
        )
    else:
        cluster_order = observed_clusters

    if cluster_order:
        ax.set_xticks(range(len(cluster_order)))
        ax.set_xticklabels(
            [str(value) for value in cluster_order], rotation=45, ha="right"
        )
    if row_order:
        ax.set_yticks(range(len(row_order)))
        ax.set_yticklabels([str(value) for value in row_order])

    if row_order and cluster_order:
        matrix = value_table.set_index(["row_id", "cluster"])[value_column].unstack(
            "cluster"
        )
        matrix = matrix.reindex(index=row_order, columns=cluster_order)
        if matrix.notna().any().any():
            image = ax.imshow(
                matrix.to_numpy(),
                aspect="auto",
                origin="upper",
                interpolation="nearest",
            )
            ax.figure.colorbar(
                image, ax=ax, label=value_column.replace("_", " ").title()
            )

    ax.set_xlabel("Cluster")
    ax.set_ylabel("Region / Condition" if level == "condition" else "Region / Sample")
    ax.set_title(title or "Shared cluster region occupancy")
    return fig, ax


def plot_read_cluster_region_association_heatmap_matplotlib(
    payload: Mapping[str, object],
    *,
    ax=None,
    title: str | None = None,
    region_sort: str | Sequence[str] | None = None,
    association_strength_aggregate: str | None = None,
    cluster_sort: str | None = None,
    cluster_order: Sequence[object] | None = None,
    region_label_mode: str = "auto",
    max_region_labels: int = 50,
    row_annotation_column: str | None = None,
    row_annotation_columns: Sequence[str] | None = None,
    row_annotation_title: str | None = None,
    row_annotation_titles: Mapping[str, str] | None = None,
    row_annotation_palette: Mapping[str, str] | None = None,
    row_annotation_palettes: Mapping[str, Mapping[str, str]] | None = None,
    group_region_labels: bool | None = None,
    group_label_columns: Sequence[str] | None = None,
):
    _require_payload_keys(
        payload,
        ("matrix_table", "metadata"),
        owner="plot_read_cluster_region_association_heatmap_matplotlib",
    )
    matrix_table = _require_payload_table(payload, "matrix_table")
    metadata = payload.get("metadata")
    if not isinstance(metadata, Mapping):
        raise TypeError("plot payload key 'metadata' must be a mapping.")

    if region_sort is not None:
        association_table = payload.get("association_table")
        if not isinstance(association_table, pd.DataFrame):
            raise ValueError(
                "region_sort override requires payload['association_table']; "
                "build payload with prepare_read_cluster_region_association_data(...)."
            )
        from dimelo import plotting as plotting_api

        value_mode = str(metadata.get("value_mode") or "fraction")
        top_n_regions_per_cluster = metadata.get("top_n_regions_per_cluster")
        aggregate_mode = (
            association_strength_aggregate
            if association_strength_aggregate is not None
            else metadata.get("association_strength_aggregate", "max")
        )
        cluster_sort_mode = (
            cluster_sort if cluster_sort is not None else metadata.get("cluster_sort", "input")
        )
        override_payload = plotting_api.prepare_read_cluster_region_association_data(
            association_table=association_table,
            value_mode=value_mode,
            top_n_regions_per_cluster=top_n_regions_per_cluster,
            region_sort=region_sort,
            association_strength_aggregate=str(aggregate_mode),
            cluster_sort=str(cluster_sort_mode),
            cluster_order=cluster_order,
        )

        matrix_table = _require_payload_table(override_payload, "matrix_table")
        metadata = override_payload.get("metadata")
        if not isinstance(metadata, Mapping):
            raise TypeError(
                "region_sort override produced payload metadata that is not a mapping."
            )

        original_axis_table = payload.get("region_axis_table")
        region_axis_table = override_payload.get("region_axis_table")
        if (
            isinstance(original_axis_table, pd.DataFrame)
            and isinstance(region_axis_table, pd.DataFrame)
            and "region_id" in original_axis_table.columns
            and "region_id" in region_axis_table.columns
        ):
            extra_cols = [
                c
                for c in original_axis_table.columns
                if c not in region_axis_table.columns
            ]
            if extra_cols:
                extras = original_axis_table.loc[
                    :, ["region_id", *extra_cols]
                ].drop_duplicates(subset=["region_id"])
                region_axis_table = region_axis_table.merge(
                    extras, on="region_id", how="left"
                )
    elif cluster_sort is not None or cluster_order is not None:
        association_table = payload.get("association_table")
        if not isinstance(association_table, pd.DataFrame):
            raise ValueError(
                "cluster_sort or cluster_order override requires payload['association_table']; "
                "build payload with prepare_read_cluster_region_association_data(...)."
            )
        from dimelo import plotting as plotting_api

        override_payload = plotting_api.prepare_read_cluster_region_association_data(
            association_table=association_table,
            value_mode=str(metadata.get("value_mode") or "fraction"),
            top_n_regions_per_cluster=metadata.get("top_n_regions_per_cluster"),
            region_sort=metadata.get("region_sort") or "cluster_fraction",
            association_strength_aggregate=str(
                metadata.get("association_strength_aggregate", "max")
            ),
            cluster_sort=str(cluster_sort or metadata.get("cluster_sort", "input")),
            cluster_order=cluster_order,
        )
        matrix_table = _require_payload_table(override_payload, "matrix_table")
        metadata = override_payload.get("metadata")
        original_axis_table = payload.get("region_axis_table")
        region_axis_table = override_payload.get("region_axis_table")
        if (
            isinstance(original_axis_table, pd.DataFrame)
            and isinstance(region_axis_table, pd.DataFrame)
            and "region_id" in original_axis_table.columns
            and "region_id" in region_axis_table.columns
        ):
            extra_cols = [
                c
                for c in original_axis_table.columns
                if c not in region_axis_table.columns
            ]
            if extra_cols:
                extras = original_axis_table.loc[
                    :, ["region_id", *extra_cols]
                ].drop_duplicates(subset=["region_id"])
                region_axis_table = region_axis_table.merge(
                    extras, on="region_id", how="left"
                )
    else:
        region_axis_table = payload.get("region_axis_table")

    if matrix_table.empty:
        fig, ax = _make_axis(ax=ax, figsize=(8, 4))
        ax.set_title(title or "Read-cluster association heatmap")
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Region")
        return fig, ax

    annotation_columns = _normalize_row_annotation_columns(
        row_annotation_column=row_annotation_column,
        row_annotation_columns=row_annotation_columns,
    )
    region_column = (
        "region_id" if "region_id" in matrix_table.columns else matrix_table.columns[0]
    )
    cluster_columns = [
        column for column in matrix_table.columns if column != region_column
    ]
    fig_width = max(8.0, 4.5 + (0.75 * len(cluster_columns)))
    if annotation_columns:
        fig_width += 1.1 + (0.45 * len(annotation_columns))
    fig_height = max(4.5, min(16.0, 2.0 + (0.10 * len(matrix_table))))
    fig, ax = _make_axis(ax=ax, figsize=(fig_width, fig_height))
    value_mode = str(metadata.get("value_mode") or "fraction")
    if region_axis_table is not None and not isinstance(
        region_axis_table, pd.DataFrame
    ):
        raise TypeError(
            "plot payload key 'region_axis_table' must be a pandas DataFrame when provided."
        )
    if region_label_mode not in {"auto", "region_id", "genomic", "chromosome"}:
        raise ValueError(
            "region_label_mode must be 'auto', 'region_id', 'genomic', or 'chromosome'."
        )
    max_region_labels = max(1, int(max_region_labels))

    heatmap = matrix_table.loc[:, cluster_columns].copy()
    image = ax.imshow(
        heatmap.to_numpy(),
        aspect="auto",
        origin="upper",
        interpolation="nearest",
    )
    colorbar = ax.figure.colorbar(
        image,
        ax=ax,
        label=value_mode.replace("_", " ").title(),
        fraction=0.045,
        pad=0.05,
    )
    ax.set_xticks(range(len(cluster_columns)))
    ax.set_xticklabels(
        [str(value) for value in cluster_columns], rotation=45, ha="right"
    )

    n_regions = len(matrix_table)
    default_ticks = np.arange(n_regions)
    default_labels = matrix_table[region_column].astype(str).tolist()
    ticks = default_ticks
    labels = default_labels
    step = max(1, int(np.ceil(n_regions / max_region_labels)))
    if step > 1:
        ticks = default_ticks[::step]
        labels = [default_labels[idx] for idx in ticks]
    annotation_values_by_column: dict[str, list[str]] = {}

    axis_table = None
    if isinstance(region_axis_table, pd.DataFrame) and not region_axis_table.empty:
        axis_table = region_axis_table.copy()
        if {"chrom", "start", "end"}.issubset(axis_table.columns):
            genomic_labels = [
                f"{chrom}:{int(start):,}-{int(end):,}"
                for chrom, start, end in axis_table[
                    ["chrom", "start", "end"]
                ].itertuples(index=False)
            ]
            chrom_labels = axis_table["chrom"].astype(str).tolist()
            effective_mode = region_label_mode
            if effective_mode == "auto":
                default_region_sort = str(metadata.get("region_sort") or "input")
                if default_region_sort == "genomic":
                    effective_mode = (
                        "chromosome" if n_regions > max_region_labels else "genomic"
                    )
                else:
                    effective_mode = "region_id"
            if effective_mode == "genomic":
                labels = [genomic_labels[idx] for idx in ticks]
            elif effective_mode == "chromosome":
                labels = [chrom_labels[idx] for idx in ticks]
                # Mark chromosome boundaries to visually group regions.
                chrom_codes = pd.Categorical(axis_table["chrom"].astype(str)).codes
                boundaries = np.flatnonzero(np.diff(chrom_codes)) + 1
                for boundary in boundaries:
                    ax.axhline(boundary - 0.5, color="white", linewidth=0.6, alpha=0.6)
            elif effective_mode == "region_id":
                labels = [default_labels[idx] for idx in ticks]
        for annotation_column in annotation_columns:
            if annotation_column in axis_table.columns:
                annotation_values_by_column[annotation_column] = (
                    axis_table[annotation_column].astype(str).tolist()
                )

    if group_region_labels is None:
        available_annotation_columns = [
            column for column in annotation_columns if column in annotation_values_by_column
        ]
        group_region_labels = (
            axis_table is not None
            and region_label_mode in {"genomic", "chromosome", "auto"}
            and n_regions > max_region_labels
            and (
                ("chrom" in axis_table.columns)
                or bool(available_annotation_columns)
            )
        )

    if axis_table is not None and group_region_labels:
        requested_group_cols = list(group_label_columns or [])
        if not requested_group_cols:
            available_annotation_columns = [
                column for column in annotation_columns if column in axis_table.columns
            ]
            if available_annotation_columns:
                requested_group_cols = [available_annotation_columns[0], "chrom"]
            elif "chrom" in axis_table.columns:
                requested_group_cols = ["chrom"]
        requested_group_cols = [
            col for col in requested_group_cols if col in axis_table.columns
        ]
        if requested_group_cols:
            group_parts = axis_table.loc[:, requested_group_cols].astype(str)
            group_key = group_parts.agg(" | ".join, axis=1).tolist()
            boundaries = (
                np.flatnonzero(np.array(group_key[1:]) != np.array(group_key[:-1])) + 1
            )
            group_starts = np.concatenate(([0], boundaries))
            group_ends = np.concatenate((boundaries, [len(group_key)]))
            centers = ((group_starts + group_ends - 1) / 2.0).astype(float)
            group_labels = [group_key[int(start)] for start in group_starts]
            # Downsample grouped labels when there are still too many.
            if len(group_labels) > max_region_labels:
                gstep = max(1, int(np.ceil(len(group_labels) / max_region_labels)))
                centers = centers[::gstep]
                group_labels = group_labels[::gstep]
            ticks = centers
            labels = group_labels
            for boundary in boundaries:
                ax.axhline(boundary - 0.5, color="white", linewidth=0.6, alpha=0.6)

    ax.set_yticks(np.asarray(ticks).tolist())
    use_legacy_annotation_labels = (
        bool(annotation_values_by_column)
        and row_annotation_columns is None
        and row_annotation_column is not None
        and len(annotation_values_by_column) == 1
        and not group_region_labels
    )
    if use_legacy_annotation_labels:
        legacy_values = next(iter(annotation_values_by_column.values()))
        labels = [f"{labels[pos]} | {legacy_values[idx]}" for pos, idx in enumerate(ticks)]
    ytick_fontsize = 8 if len(labels) <= 30 else 7
    ax.set_yticklabels(labels, fontsize=ytick_fontsize)
    if labels:
        max_label_len = max(len(str(label)) for label in labels)
        left_margin = min(0.42, max(0.14, 0.08 + (0.0035 * max_label_len)))
        fig.subplots_adjust(left=left_margin)
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Region")
    ax.set_title(
        title
        or _read_cluster_association_title(
            metadata,
            row_annotation_columns=annotation_columns,
            group_region_labels=group_region_labels,
            region_label_mode=region_label_mode,
        ),
        pad=14,
    )

    if annotation_values_by_column and n_regions > 0:
        from matplotlib.colors import ListedColormap
        from matplotlib.patches import Patch

        ax.set_ylabel("")
        # Place the annotation strip in figure coordinates so it never overlaps
        # the main heatmap region, regardless of y-label length.
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        ax_pos = ax.get_position()
        tick_labels = [tick for tick in ax.get_yticklabels() if tick.get_text()]
        if tick_labels:
            tick_left = min(
                tick.get_window_extent(renderer)
                .transformed(fig.transFigure.inverted())
                .x0
                for tick in tick_labels
            )
        else:
            tick_left = ax_pos.x0

        strip_width = 0.012
        strip_gap = 0.006
        label_gap = 0.018
        n_annotation_strips = len(annotation_values_by_column)
        total_strip_width = (
            n_annotation_strips * strip_width
            + max(0, n_annotation_strips - 1) * strip_gap
        )
        strip_x1 = tick_left - label_gap
        first_strip_x0 = strip_x1 - total_strip_width
        min_left = 0.012
        if first_strip_x0 < min_left:
            needed_left = min(0.56, ax_pos.x0 + (min_left - first_strip_x0) + 0.02)
            fig.subplots_adjust(left=needed_left)
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            ax_pos = ax.get_position()
            tick_labels = [tick for tick in ax.get_yticklabels() if tick.get_text()]
            if tick_labels:
                tick_left = min(
                    tick.get_window_extent(renderer)
                    .transformed(fig.transFigure.inverted())
                    .x0
                    for tick in tick_labels
                )
            else:
                tick_left = ax_pos.x0
            strip_x1 = tick_left - label_gap
            first_strip_x0 = max(min_left, strip_x1 - total_strip_width)

        plt = _import_pyplot()
        cmap = plt.get_cmap("tab10")
        for strip_idx, (annotation_column, annotation_values) in enumerate(
            annotation_values_by_column.items()
        ):
            ordered_annotations = _ordered_unique_values(
                pd.DataFrame({"annotation": annotation_values}),
                "annotation",
            )
            annotation_palette = _row_annotation_palette(
                annotation_column,
                row_annotation_palette=row_annotation_palette,
                row_annotation_palettes=row_annotation_palettes,
                n_columns=len(annotation_values_by_column),
            )
            if annotation_palette is None:
                annotation_palette = {
                    value: cmap(i % 10) for i, value in enumerate(ordered_annotations)
                }
            annotation_colors = [
                annotation_palette.get(value, "0.6") for value in ordered_annotations
            ]
            color_lookup = {value: i for i, value in enumerate(ordered_annotations)}
            color_codes = np.array(
                [color_lookup.get(value, 0) for value in annotation_values], dtype=int
            )
            strip_x0 = first_strip_x0 + strip_idx * (strip_width + strip_gap)
            strip_ax = fig.add_axes(
                [strip_x0, ax_pos.y0, strip_width, ax_pos.height]
            )
            strip_ax.imshow(
                color_codes.reshape(-1, 1),
                aspect="auto",
                origin="upper",
                interpolation="nearest",
                cmap=ListedColormap(annotation_colors),
                vmin=0,
                vmax=max(1, len(annotation_colors) - 1),
            )
            strip_ax.set_xticks([])
            strip_ax.set_yticks([])
            strip_title = _row_annotation_title(
                annotation_column,
                row_annotation_title=row_annotation_title,
                row_annotation_titles=row_annotation_titles,
                n_columns=len(annotation_values_by_column),
            )
            strip_ax.set_xlabel(strip_title, fontsize=7, labelpad=8)
            strip_ax.xaxis.label.set_rotation(90)
            strip_ax.xaxis.label.set_verticalalignment("top")
            strip_ax.xaxis.label.set_horizontalalignment("center")

            legend_handles = [
                Patch(
                    facecolor=annotation_palette.get(value, "0.6"),
                    edgecolor="none",
                    label=str(value),
                )
                for value in ordered_annotations
            ]
            if legend_handles:
                fig.canvas.draw()
                legend_x = max(
                    colorbar.ax.get_position().x1 + 0.055,
                    ax.get_position().x1 + 0.12,
                )
                legend = fig.legend(
                    handles=legend_handles,
                    title=strip_title,
                    loc="upper left",
                    bbox_to_anchor=(legend_x, ax.get_position().y1 - (0.16 * strip_idx)),
                    bbox_transform=fig.transFigure,
                    frameon=False,
                )
                fig.add_artist(legend)
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

    table = _require_payload_table(
        payload, "sample_summary" if level == "sample" else "condition_summary"
    )
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
    metadata = (
        payload.get("metadata", {})
        if isinstance(payload.get("metadata", {}), Mapping)
        else {}
    )
    contigs = list(
        metadata.get("contig_order")
        or window_table.get("contig", pd.Series(dtype="object")).dropna().unique()
    )

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
            ax.plot(
                grouped["window_midpoint"],
                grouped["window_fraction"],
                marker="o",
                linewidth=1.5,
            )
        ax.set_title(str(contig))
        ax.set_xlabel("Window midpoint")
        ax.set_ylabel("Window Fraction")

    for ax in axes[len(contigs) :]:
        ax.set_visible(False)

    if title:
        fig.suptitle(title)
    return fig, axes
