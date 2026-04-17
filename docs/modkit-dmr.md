# modkit DMR Workflows

`dimelo` now exposes native wrappers for `modkit dmr` so region discovery and differential methylation can be driven directly by ONT's maintained DMR implementation.

## Pairwise DMR (+ optional HMM segmentation)

Use `workflows.modkit_dmr_pair_workflow(...)` (or `dimelo.dmr.run_dmr_pair(...)`) for two-sample comparison.

- Site-level output includes `map_pvalue` and `effect_size` columns from `modkit dmr pair`.
- Optional HMM segmentation is enabled by passing `segment_path=...`.
- `regions_bed` and `segment_path` are mutually exclusive in a single run because `modkit dmr pair` enforces that constraint.

Example:

```python
from dimelo import workflows

result = workflows.modkit_dmr_pair_workflow(
    control_bed_methyl="sample_a.pileup.sorted.bed.gz",
    experiment_bed_methyl="sample_b.pileup.sorted.bed.gz",
    ref_genome="ref.fa",
    out_path="dmr_sites.bed",
    segment_path="dmr_segments.bed",
    bases=["A"],
    threads=8,
)

print(result.sites.head())
print(result.segments.head())
print(result.high_confidence_sites.head())
```

## Multi-sample regional DMR

Use `workflows.modkit_dmr_multi_workflow(...)` when comparing all sample pairs over a fixed region BED.

Example:

```python
multi = workflows.modkit_dmr_multi_workflow(
    samples={
        "s1": "s1.pileup.sorted.bed.gz",
        "s2": "s2.pileup.sorted.bed.gz",
        "s3": "s3.pileup.sorted.bed.gz",
    },
    regions_bed="regions.bed",
    ref_genome="ref.fa",
    out_dir="dmr_multi",
    bases=["A"],
    threads=8,
)

print(multi.pair_files)
```

If you already use `SampleSpec`, you can route through:

`workflows.modkit_dmr_multi_from_samples_workflow(...)` (reads `sample.metadata["pileup_path"]`).
