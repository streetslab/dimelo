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


def _ordered_cluster_labels(metadata: Mapping[str, object], observed_clusters: list[object]) -> list[object]:
    cluster_labels = list(metadata.get("cluster_labels") or [])
    if not cluster_labels:
        return list(observed_clusters)

    observed_cluster_set = set(observed_clusters)
    ordered_clusters = [cluster for cluster in cluster_labels if cluster in observed_cluster_set]
    ordered_clusters.extend([cluster for cluster in observed_clusters if cluster not in cluster_labels])
    return ordered_clusters


def _ordered_unique_values(table: pd.DataFrame, column: str) -> list[object]:
    return table.loc[:, column].drop_duplicates().tolist()


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
    _require_payload_keys(payload, ("metadata",), owner="plot_shared_cluster_distribution_matplotlib")
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
            raise ValueError("plot payload contains duplicate cluster fractions for the same x value.")
        x_order = _ordered_unique_values(value_table, x_column)
        cluster_order = _ordered_cluster_labels(metadata, _ordered_unique_values(value_table, "cluster"))
        pivot = value_table.set_index([x_column, "cluster"])["fraction"].unstack("cluster", fill_value=0.0)
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
    _require_payload_keys(payload, ("metadata",), owner="plot_shared_cluster_change_matplotlib")
    change_table = _require_payload_table(payload, "distribution_change")
    metadata = payload.get("metadata")
    if not isinstance(metadata, Mapping):
        raise TypeError("plot payload key 'metadata' must be a mapping.")
    fig, ax = _make_axis(ax=ax, figsize=(8, 4))

    if not change_table.empty:
        value_table = change_table.loc[:, ["condition", "cluster", "delta_fraction"]].copy()
        if value_table.duplicated(["condition", "cluster"]).any():
            raise ValueError("plot payload contains duplicate delta fractions for the same condition.")
        condition_order = list(metadata.get("change_condition_order") or [])
        if not condition_order:
            condition_order = _ordered_unique_values(value_table, "condition")
        observed_clusters = _ordered_unique_values(value_table, "cluster")
        cluster_order = list(metadata.get("cluster_labels") or [])
        if cluster_order:
            cluster_order.extend(
                [cluster for cluster in observed_clusters if cluster not in cluster_order]
            )
        else:
            cluster_order = observed_clusters
        matrix = value_table.set_index(["condition", "cluster"])["delta_fraction"].unstack("cluster")
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
        ax.set_xticklabels([str(value) for value in matrix.columns.tolist()], rotation=45, ha="right")
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
        raise ValueError("plot payload contains duplicate profile values for the same cluster and feature.")

    feature_order = list(metadata.get("feature_names") or [])
    observed_features = _ordered_unique_values(value_table, "feature")
    if feature_order:
        feature_order.extend([feature for feature in observed_features if feature not in feature_order])
    else:
        feature_order = observed_features

    observed_clusters = _ordered_unique_values(value_table, "cluster")
    cluster_order = list(metadata.get("cluster_labels") or [])
    if cluster_order:
        cluster_order.extend([cluster for cluster in observed_clusters if cluster not in cluster_order])
    else:
        cluster_order = observed_clusters

    if feature_order:
        ax.set_xticks(range(len(feature_order)))
        ax.set_xticklabels([str(value) for value in feature_order], rotation=45, ha="right")
    if cluster_order:
        ax.set_yticks(range(len(cluster_order)))
        ax.set_yticklabels([str(value) for value in cluster_order])

    if feature_order and cluster_order:
        if value_table.empty:
            matrix = pd.DataFrame(index=cluster_order, columns=feature_order, dtype=float)
        else:
            matrix = value_table.set_index(["cluster", "feature"])["value"].unstack("feature")
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
            ax.set_xticklabels([str(value) for value in matrix.columns.tolist()], rotation=45, ha="right")
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
        raise ValueError("plot payload contains duplicate profile values for the same cluster and feature.")

    feature_order = list(metadata.get("feature_names") or [])
    observed_features = _ordered_unique_values(value_table, "feature")
    if feature_order:
        feature_order.extend([feature for feature in observed_features if feature not in feature_order])
    else:
        feature_order = observed_features
    feature_lookup = {feature: index for index, feature in enumerate(feature_order)}

    observed_clusters = _ordered_unique_values(value_table, "cluster")
    cluster_order = list(metadata.get("cluster_labels") or [])
    if cluster_order:
        cluster_order.extend([cluster for cluster in observed_clusters if cluster not in cluster_order])
    else:
        cluster_order = observed_clusters

    for cluster in cluster_order:
        cluster_table = value_table.loc[value_table["cluster"] == cluster, ["feature", "value"]].copy()
        if feature_order:
            cluster_table["feature_order"] = cluster_table["feature"].map(feature_lookup)
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
    value_table = table.loc[:, ["region_id", label_column, "cluster", value_column]].copy()
    value_table["row_id"] = value_table.loc[:, ["region_id", label_column]].astype(str).agg(" | ".join, axis=1)
    if value_table.duplicated(["row_id", "cluster"]).any():
        raise ValueError("plot payload contains duplicate region occupancy values for the same row and cluster.")

    row_order = _ordered_unique_values(value_table, "row_id")
    observed_clusters = _ordered_unique_values(value_table, "cluster")
    cluster_order = list(metadata.get("cluster_labels") or [])
    if cluster_order:
        cluster_order.extend([cluster for cluster in observed_clusters if cluster not in cluster_order])
    else:
        cluster_order = observed_clusters

    if cluster_order:
        ax.set_xticks(range(len(cluster_order)))
        ax.set_xticklabels([str(value) for value in cluster_order], rotation=45, ha="right")
    if row_order:
        ax.set_yticks(range(len(row_order)))
        ax.set_yticklabels([str(value) for value in row_order])

    if row_order and cluster_order:
        matrix = value_table.set_index(["row_id", "cluster"])[value_column].unstack("cluster")
        matrix = matrix.reindex(index=row_order, columns=cluster_order)
        if matrix.notna().any().any():
            image = ax.imshow(
                matrix.to_numpy(),
                aspect="auto",
                origin="upper",
                interpolation="nearest",
            )
            ax.figure.colorbar(image, ax=ax, label=value_column.replace("_", " ").title())

    ax.set_xlabel("Cluster")
    ax.set_ylabel("Region / Condition" if level == "condition" else "Region / Sample")
    ax.set_title(title or "Shared cluster region occupancy")
    return fig, ax


def plot_read_cluster_region_association_heatmap_matplotlib(
    payload: Mapping[str, object],
    *,
    ax=None,
    title: str | None = None,
    region_label_mode: str = "auto",
    max_region_labels: int = 50,
    row_annotation_column: str | None = None,
    row_annotation_title: str | None = None,
    row_annotation_palette: Mapping[str, str] | None = None,
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

    if matrix_table.empty:
        fig, ax = _make_axis(ax=ax, figsize=(8, 4))
        ax.set_title(title or "Read-cluster association heatmap")
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Region")
        return fig, ax

    region_column = "region_id" if "region_id" in matrix_table.columns else matrix_table.columns[0]
    cluster_columns = [column for column in matrix_table.columns if column != region_column]
    fig_width = max(8.0, 4.5 + (0.75 * len(cluster_columns)))
    fig_height = max(4.5, min(16.0, 2.0 + (0.10 * len(matrix_table))))
    fig, ax = _make_axis(ax=ax, figsize=(fig_width, fig_height))
    value_mode = str(metadata.get("value_mode") or "fraction")
    region_axis_table = payload.get("region_axis_table")
    if region_axis_table is not None and not isinstance(region_axis_table, pd.DataFrame):
        raise TypeError("plot payload key 'region_axis_table' must be a pandas DataFrame when provided.")
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
    ax.figure.colorbar(image, ax=ax, label=value_mode.replace("_", " ").title())
    ax.set_xticks(range(len(cluster_columns)))
    ax.set_xticklabels([str(value) for value in cluster_columns], rotation=45, ha="right")

    n_regions = len(matrix_table)
    default_ticks = np.arange(n_regions)
    default_labels = matrix_table[region_column].astype(str).tolist()
    ticks = default_ticks
    labels = default_labels
    step = max(1, int(np.ceil(n_regions / max_region_labels)))
    if step > 1:
        ticks = default_ticks[::step]
        labels = [default_labels[idx] for idx in ticks]
    annotation_values: list[str] | None = None

    axis_table = None
    if isinstance(region_axis_table, pd.DataFrame) and not region_axis_table.empty:
        axis_table = region_axis_table.copy()
        if {"chrom", "start", "end"}.issubset(axis_table.columns):
            genomic_labels = [
                f"{chrom}:{int(start):,}-{int(end):,}"
                for chrom, start, end in axis_table[["chrom", "start", "end"]].itertuples(index=False)
            ]
            chrom_labels = axis_table["chrom"].astype(str).tolist()
            effective_mode = region_label_mode
            if effective_mode == "auto":
                default_region_sort = str(metadata.get("region_sort") or "input")
                if default_region_sort == "genomic":
                    effective_mode = "chromosome" if n_regions > max_region_labels else "genomic"
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
        if row_annotation_column is not None and row_annotation_column in axis_table.columns:
            annotation_values = axis_table[row_annotation_column].astype(str).tolist()

    if group_region_labels is None:
        group_region_labels = (
            axis_table is not None
            and region_label_mode in {"genomic", "chromosome", "auto"}
            and n_regions > max_region_labels
            and (
                ("chrom" in axis_table.columns)
                or (row_annotation_column is not None and row_annotation_column in axis_table.columns)
            )
        )

    if axis_table is not None and group_region_labels:
        requested_group_cols = list(group_label_columns or [])
        if not requested_group_cols:
            if row_annotation_column is not None and row_annotation_column in axis_table.columns:
                requested_group_cols = [row_annotation_column, "chrom"]
            elif "chrom" in axis_table.columns:
                requested_group_cols = ["chrom"]
        requested_group_cols = [col for col in requested_group_cols if col in axis_table.columns]
        if requested_group_cols:
            group_parts = axis_table.loc[:, requested_group_cols].astype(str)
            group_key = group_parts.agg(" | ".join, axis=1).tolist()
            boundaries = np.flatnonzero(np.array(group_key[1:]) != np.array(group_key[:-1])) + 1
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
    if annotation_values is not None and not group_region_labels:
        labels = [f"{labels[pos]} | {annotation_values[idx]}" for pos, idx in enumerate(ticks)]
    ytick_fontsize = 8 if len(labels) <= 30 else 7
    ax.set_yticklabels(labels, fontsize=ytick_fontsize)
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Region")
    ax.set_title(title or "Read-cluster association heatmap")

    if annotation_values is not None and n_regions > 0:
        from matplotlib.colors import ListedColormap
        from matplotlib.patches import Patch

        ordered_annotations = _ordered_unique_values(
            pd.DataFrame({"annotation": annotation_values}),
            "annotation",
        )
        if row_annotation_palette is None:
            plt = _import_pyplot()
            cmap = plt.get_cmap("tab10")
            row_annotation_palette = {
                value: cmap(i % 10) for i, value in enumerate(ordered_annotations)
            }
        annotation_colors = [row_annotation_palette.get(value, "0.6") for value in ordered_annotations]
        color_lookup = {value: i for i, value in enumerate(ordered_annotations)}
        color_codes = np.array([color_lookup.get(value, 0) for value in annotation_values], dtype=int)
        strip_ax = ax.inset_axes([-0.07, 0.0, 0.022, 1.0], transform=ax.transAxes)
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
        strip_ax.set_title(
            row_annotation_title or row_annotation_column,
            fontsize=8,
            pad=2,
        )

        legend_handles = [
            Patch(facecolor=row_annotation_palette.get(value, "0.6"), edgecolor="none", label=str(value))
            for value in ordered_annotations
        ]
        if legend_handles:
            ax.legend(
                handles=legend_handles,
                title=row_annotation_title or row_annotation_column,
                loc="upper left",
                bbox_to_anchor=(1.01, 1.0),
                frameon=False,
            )
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
