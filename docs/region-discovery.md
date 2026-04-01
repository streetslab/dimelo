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

## Handoff Guidance

- For formal region testing, convert discovered hits into BED and write them to disk before passing the BED path into `region_contrasts.score_regions(...)`.
- This handoff works the same for paired discovery outputs; the paired score columns stay in `result.hits`, while `hits_to_bed()` projects the ranked loci into a plain BED table for downstream region testing.
- For downstream clustering, use the discovered loci to focus the region set before running the clustering workflow.
- Keep discovery and contrast roles separate: discovery finds candidates, contrasts tests known regions, and clustering explains state mixtures.

The combined discovery-to-clustering workflow uses the same handoff, but keeps it in memory:

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

```python
bed_df = region_discovery.hits_to_bed(result.hits)
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
