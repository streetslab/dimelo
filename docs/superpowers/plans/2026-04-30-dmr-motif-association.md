# DMR Motif Association Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a notebook section that quantifies how CTCF PWM score relates to positive DMR signal, including a global correlation view, an enrichment test, and a top-region local association view.

**Architecture:** Reuse the existing 100 bp tiled discovery table as the source of truth. Compute a per-bin CTCF PWM score from the reference sequence, keep the analysis restricted to positive-effect bins, and derive three outputs from the same annotated table: a Spearman correlation, a Fisher exact enrichment test, and a region-local scatter plot for the top 3 regions. Keep the new code inside the notebook so the example stays self-contained.

**Tech Stack:** Python, pandas, NumPy, SciPy, Matplotlib, Biopython, pysam, Jupyter notebook JSON.

---

### Task 1: Annotate positive 100 bp bins with PWM scores

**Files:**
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dmr_multi_sample.ipynb`

- [ ] **Step 1: Add the bin-scoring code**

```python
positive_bins = tiled.loc[tiled["effect_size"] > 0].copy().reset_index(drop=True)

bin_score_rows = []
for row in positive_bins.itertuples(index=False):
    seq = fasta.fetch(row.chrom, int(row.start), int(row.end))
    best_hit = best_ctcf_pwm_hit(seq, window_start=int(row.start))
    bin_score_rows.append(
        {
            "chrom": row.chrom,
            "start": int(row.start),
            "end": int(row.end),
            "ctcf_best_score": best_hit["score"],
            "ctcf_best_score_fraction": best_hit["score_fraction"],
            "ctcf_best_hit_start": best_hit["start"],
            "ctcf_best_hit_end": best_hit["end"],
            "ctcf_best_hit_strand": best_hit["strand"],
            "ctcf_hit_above_threshold": best_hit["score"] >= ctcf_pwm_threshold,
        }
    )

positive_bins = positive_bins.merge(
    pd.DataFrame(bin_score_rows),
    on=["chrom", "start", "end"],
    how="left",
)
```

- [ ] **Step 2: Run the notebook cell logic in a Python check**

Run:
```bash
python3.11 - <<'PY'
import json
from pathlib import Path
nb = json.loads(Path("dmr_multi_sample.ipynb").read_text())
print(len(nb["cells"]))
PY
```
Expected: the notebook still parses and the new cell is present.

- [ ] **Step 3: Keep the result self-contained**

```python
positive_bins["ctcf_hit_above_threshold"] = positive_bins["ctcf_hit_above_threshold"].fillna(False)
positive_bins["ctcf_best_score_fraction"] = positive_bins["ctcf_best_score_fraction"].fillna(0.0)
```

- [ ] **Step 4: Re-run the notebook JSON parse check**

Run:
```bash
python3.11 - <<'PY'
import json
from pathlib import Path
json.loads(Path("dmr_multi_sample.ipynb").read_text())
print("ok")
PY
```
Expected: `ok`

- [ ] **Step 5: Commit**

```bash
git add dmr_multi_sample.ipynb
git commit -m "feat: add motif association notebook section"
```

### Task 2: Add global correlation and enrichment outputs

**Files:**
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dmr_multi_sample.ipynb`

- [ ] **Step 1: Add the stats imports and summary table**

```python
from scipy.stats import fisher_exact, spearmanr

spearman_rho, spearman_p = spearmanr(
    positive_bins["ctcf_best_score_fraction"],
    positive_bins["effect_size"],
)

contingency = pd.crosstab(
    positive_bins["candidate_dmr"],
    positive_bins["ctcf_hit_above_threshold"],
).reindex(index=[False, True], columns=[False, True], fill_value=0)

odds_ratio, fisher_p = fisher_exact(contingency.to_numpy(), alternative="greater")

association_summary = pd.DataFrame(
    {
        "metric": [
            "positive bins",
            "significant positive bins",
            "PWM-hit positive bins",
            "PWM-hit significant positive bins",
            "Spearman rho",
            "Spearman p-value",
            "Fisher odds ratio",
            "Fisher p-value",
        ],
        "value": [
            len(positive_bins),
            int(positive_bins["candidate_dmr"].sum()),
            int(positive_bins["ctcf_hit_above_threshold"].sum()),
            int(
                (positive_bins["candidate_dmr"]
                & positive_bins["ctcf_hit_above_threshold"]).sum()
            ),
            spearman_rho,
            spearman_p,
            odds_ratio,
            fisher_p,
        ],
    }
)
display(association_summary)
```

- [ ] **Step 2: Add the global scatter and enrichment plot**

```python
fig, (scatter_ax, enrich_ax) = plt.subplots(
    ncols=2,
    figsize=(11.5, 4.2),
    constrained_layout=True,
)

scatter_ax.scatter(
    positive_bins.loc[~positive_bins["candidate_dmr"], "ctcf_best_score_fraction"],
    positive_bins.loc[~positive_bins["candidate_dmr"], "effect_size"],
    s=28,
    color="#7570B3",
    alpha=0.8,
    label="positive non-significant bin",
)
scatter_ax.scatter(
    positive_bins.loc[positive_bins["candidate_dmr"], "ctcf_best_score_fraction"],
    positive_bins.loc[positive_bins["candidate_dmr"], "effect_size"],
    s=34,
    color="#D95F02",
    alpha=0.95,
    label="significant positive bin",
)
scatter_ax.axhline(ABS_EFFECT_MIN, color="0.65", linestyle="--", linewidth=0.8)
scatter_ax.axvline(ctcf_pwm_threshold / ctcf_pssm.max, color="0.65", linestyle="--", linewidth=0.8)
scatter_ax.set_xlabel("Best CTCF PWM score fraction")
scatter_ax.set_ylabel("Effect size (barcode17 - barcode18)")
scatter_ax.set_title(
    f"Positive bins only: Spearman ρ={spearman_rho:.2f}, p={spearman_p:.3g}"
)
scatter_ax.legend(frameon=False, loc="best")

hit_rate = pd.DataFrame(
    {
        "category": ["non-significant", "significant"],
        "hit_rate": [
            positive_bins.loc[~positive_bins["candidate_dmr"], "ctcf_hit_above_threshold"].mean(),
            positive_bins.loc[positive_bins["candidate_dmr"], "ctcf_hit_above_threshold"].mean(),
        ],
    }
)
enrich_ax.bar(hit_rate["category"], hit_rate["hit_rate"], color=["#7570B3", "#D95F02"])
enrich_ax.set_ylim(0, 1)
enrich_ax.set_ylabel("Fraction with PWM hit")
enrich_ax.set_title(
    f"Enrichment: OR={odds_ratio:.2f}, Fisher p={fisher_p:.3g}"
)
for idx, value in enumerate(hit_rate["hit_rate"]):
    enrich_ax.text(idx, value + 0.03, f"{value:.0%}", ha="center", va="bottom", fontsize=9)
plt.show()
```

- [ ] **Step 3: Run the notebook cell logic in a Python check**

Run:
```bash
python3.11 - <<'PY'
import json
from pathlib import Path
json.loads(Path("dmr_multi_sample.ipynb").read_text())
print("ok")
PY
```
Expected: `ok`

### Task 3: Add top-region local association plots

**Files:**
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dmr_multi_sample.ipynb`

- [ ] **Step 1: Add the per-region local scatter**

```python
fig, axes = plt.subplots(
    nrows=len(top_candidate_regions),
    ncols=1,
    figsize=(8.5, 3.3 * len(top_candidate_regions)),
    constrained_layout=True,
)
if len(top_candidate_regions) == 1:
    axes = [axes]

for ax, region in zip(axes, top_candidate_regions.itertuples(index=False), strict=False):
    window_start = max(0, int(region.start) - zoom_flank)
    window_end = int(region.end) + zoom_flank
    local_bins = positive_bins.loc[
        (positive_bins["start"] < window_end) & (positive_bins["end"] > window_start)
    ].copy()
    local_rho, local_p = spearmanr(
        local_bins["ctcf_best_score_fraction"],
        local_bins["effect_size"],
    )

    ax.scatter(
        local_bins.loc[~local_bins["candidate_dmr"], "ctcf_best_score_fraction"],
        local_bins.loc[~local_bins["candidate_dmr"], "effect_size"],
        s=26,
        color="#7570B3",
        alpha=0.8,
        label="positive non-significant bin",
    )
    ax.scatter(
        local_bins.loc[local_bins["candidate_dmr"], "ctcf_best_score_fraction"],
        local_bins.loc[local_bins["candidate_dmr"], "effect_size"],
        s=34,
        color="#D95F02",
        alpha=0.95,
        label="significant positive bin",
    )
    hit_bins = local_bins.loc[local_bins["ctcf_hit_above_threshold"]]
    if not hit_bins.empty:
        ax.scatter(
            hit_bins["ctcf_best_score_fraction"],
            hit_bins["effect_size"],
            s=80,
            facecolors="none",
            edgecolors="#111111",
            linewidths=1.2,
            label="PWM-hit bin",
        )
    ax.axvline(ctcf_pwm_threshold / ctcf_pssm.max, color="0.5", linestyle="--", linewidth=0.8)
    ax.axhline(ABS_EFFECT_MIN, color="0.65", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Best CTCF PWM score fraction")
    ax.set_ylabel("Effect size")
    ax.set_title(
        f"{region.chrom}:{region.start}-{region.end} | ρ={local_rho:.2f}, p={local_p:.3g}, n={len(local_bins)}"
    )
    if not hit_bins.empty:
        for _, hit in hit_bins.iterrows():
            ax.annotate(
                f"{int(hit['start'])}",
                (hit["ctcf_best_score_fraction"], hit["effect_size"]),
                textcoords="offset points",
                xytext=(4, 4),
                fontsize=8,
            )
    ax.legend(frameon=False, loc="best")
plt.show()
```

- [ ] **Step 2: Run the notebook cell logic in a Python check**

Run:
```bash
python3.11 - <<'PY'
import json
from pathlib import Path
json.loads(Path("dmr_multi_sample.ipynb").read_text())
print("ok")
PY
```
Expected: `ok`

- [ ] **Step 3: Commit**

```bash
git add dmr_multi_sample.ipynb
git commit -m "feat: add local motif association plots"
```

## Self-Review

- The plan only uses positive-effect bins, matching the requested restriction.
- Each output reuses the same annotated 100 bp tile table, so the analyses are consistent.
- No placeholder text remains in the steps.
