# Region Discovery

`dimelo.region_discovery` is the de novo locus discovery layer. Use it when you do not yet know which genomic intervals matter and you want a deterministic tiled scan over pileup-backed inputs. It also supports paired discovery when you have matched samples or ordered time courses and want to keep the output in the same data-first shape as pooled discovery.

## When To Use Discovery vs Contrasts

- Use `region_discovery` when you need to find candidate loci first.
- Use `region_contrasts` when you already know the regions and want to test them formally.
- If you already have a BED file or a matched region set, skip discovery and go straight to `region_contrasts`.

## Paired Discovery

Paired discovery is still discovery-first. The scan resolves complete matched units from `SampleSpec.metadata[pairing_key]`, filters out incomplete pairs, and returns the same `RegionDiscoveryResult` container as pooled discovery.

- Use `ContrastSpec(mode="matched_pairwise", numerator=[...], denominator=[...], pairing_key="...")` when you have matched two-condition samples.
- Use `ContrastSpec(mode="time_course", time_order=[...], pairing_key="...")` when you want an ordered trajectory across paired time points.
- Inspect `result.windows` for the full tiled table and `result.hits` for ranked discovery hits.
- Use `result.metadata["pairing_key"]`, `result.metadata["pairing_policy"]`, `result.metadata["n_pairs_used"]`, and `result.metadata["n_pairs_dropped"]` to audit the pairing resolution that was applied.
- `result.plot_data` stays optional; the canonical handoff artifact is still the table data.

## Minimal Example

```python
from dimelo import region_discovery
from dimelo.models import ContrastSpec, SampleSpec

result = region_discovery.scan_genome(
    samples=[
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="output/s1.extract.h5",
            metadata={"pileup_path": "output/s1.pileup.sorted.bed.gz"},
        ),
        SampleSpec(
            sample_id="s2",
            condition="treated",
            extract_h5="output/s2.extract.h5",
            metadata={"pileup_path": "output/s2.pileup.sorted.bed.gz"},
        ),
    ],
    motifs=["A,0"],
    genome_sizes={"chr1": 248956422},
    window_size=2000,
    step_size=500,
    contrast=ContrastSpec(
        mode="pairwise",
        numerator=["treated"],
        denominator=["NS"],
        reference_condition="NS",
    ),
    score="beta_binomial",
    min_coverage=10,
)
```

The canonical outputs are data-first:

- `result.windows`
- `result.hits`
- `result.plot_data["window_score_table"]`
- `result.plot_data["top_hits_table"]`

As with the newer clustering and contrast workflows, these plotting payloads are renderer-neutral tables first. Legacy orientation flags such as `regions_5to3prime` map onto the shared region-alignment model for downstream follow-up plots, while any scaled segment normalization is kept to aggregate views rather than single-read displays.

## Plotting Prep Helpers

Use the shared plotting helpers when you want plot-ready discovery tables without committing to a specific renderer:

```python
from dimelo import plotting

scan_payload = plotting.prepare_region_discovery_scan_data(
    result=result,
    top_n_hits=50,
)

context_payload = plotting.prepare_region_discovery_hit_context_data(
    result=result,
    top_n=12,
    padding_windows=5,
)
```

If you want a built-in Matplotlib figure after preparing the payload, pass it through `dimelo.plotting_matplotlib`:

```python
from dimelo import plotting_matplotlib

fig, axes = plotting_matplotlib.plot_region_discovery_scan_matplotlib(scan_payload)
plotting_matplotlib.save_figure(fig, "region-discovery-scan.png")
```

- `prepare_region_discovery_scan_data(...)` consumes a `RegionDiscoveryResult` and returns renderer-neutral `scan_table`, `hit_table`, and `metadata` payloads. Scan payloads stay per-contig by default rather than projecting windows onto a cumulative genome axis.
- `prepare_region_discovery_hit_context_data(...)` consumes the same `RegionDiscoveryResult` and returns renderer-neutral `context_table`, `selected_hits`, and `metadata` payloads for small-multiple or local hit-context views.

For broad whole-sample summaries instead of locus discovery, use [global analysis](global-analysis.md). For defined-region testing after discovery, use [region contrasts](region-contrasts.md).

## Handoff Guidance

- For formal region testing, convert discovered hits into BED and write them to disk before passing the BED path into `region_contrasts.score_regions(...)`.
- This handoff works the same for paired discovery outputs; the paired score columns stay in `result.hits`, while `hits_to_bed()` projects the ranked loci into a plain BED table for downstream region testing.
- For downstream clustering, the recommended next step is read-window clustering on the discovered loci, followed by region-level occupancy follow-up from the resulting read-cluster labels.
- If you specifically want shape-based exploratory analysis, you can still pass the discovered loci into the optional `mode="region_anchored"` workflow.
- Keep discovery and contrast roles separate: discovery finds candidates, contrasts tests known regions, and clustering explains state mixtures.

The combined discovery-to-clustering workflow uses the same handoff, but keeps it in memory. It is most useful when you explicitly want the optional region-anchored exploratory path on discovery hits:

- discovery produces ranked hit rows
- `hits_to_bed()` turns those rows into BED-style `selected_regions`
- `discovery_cluster_workflow()` derives a serializable region-spec from those rows and passes it to `shared_cluster_distribution(..., matched_regions=...)`
- the workflow returns the original discovery result, the BED-style `selected_regions` table, and the clustering result together

```python
from dimelo import workflows

result = workflows.discovery_cluster_workflow(
    samples=samples,
    motifs=["A,0"],
    genome_sizes={"chr1": 248956422},
    discovery={
        "window_size": 2000,
        "step_size": 500,
        "score": "beta_binomial",
        "contrast": contrast,
    },
    clustering={"mode": "region_anchored", "n_clusters": 6},
    selection={"mode": "top_n", "top_n": 250},
)
```

When you want discovery, clustering, and defined-region contrasts in one pass, use `workflows.discovery_cluster_contrast_workflow()`:

```python
result = workflows.discovery_cluster_contrast_workflow(
    samples=samples,
    motifs=["A,0"],
    genome_sizes={"chr1": 248956422},
    discovery={
        "window_size": 2000,
        "step_size": 500,
        "score": "beta_binomial",
        "contrast": discovery_contrast,
    },
    clustering={
        "mode": "region_anchored",
        "n_clusters": 6,
    },
    contrasts={
        "contrast": region_contrast,
        "test": "beta_binomial",
    },
    selection={"mode": "top_n", "top_n": 250},
)

selected = result.selected_regions
clustered = result.clustering.assignments
scored = result.contrasts.regions
windows = result.metadata["full_scan_windows"]
```

The contract is intentionally simple:

- `result.selected_regions` is the BED-style selected follow-up set
- clustering receives a serializable region-spec derived from those rows
- clustering-side `region_id` values in the combined result are normalized to the same `chr:start-end,strand` key used by contrasts
- `result.contrasts` scores the same selected loci by default, or `contrasts["regions"]` if you provide an explicit override
- `result.metadata["full_scan_windows"]` carries the full discovery scan for context

```python
bed_df = region_discovery.hits_to_bed(result.discovery.hits)
bed_df.to_csv("discovered_hits.bed", sep="\t", header=False, index=False)
```

Paired discovery example:

```python
paired_result = region_discovery.scan_genome(
    samples=paired_samples,
    motifs=["A,0"],
    genome_sizes=genome_sizes,
    window_size=2000,
    step_size=500,
    contrast=ContrastSpec(
        mode="matched_pairwise",
        numerator=["targeting"],
        denominator=["nontargeting"],
        pairing_key="pair_id",
    ),
    score="effect_size_only",
)

paired_bed = region_discovery.hits_to_bed(paired_result.hits)
paired_bed.to_csv("paired_discovered_hits.bed", sep="\t", header=False, index=False)
```
