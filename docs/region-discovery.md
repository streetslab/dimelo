# Region Discovery

`dimelo.region_discovery` is the de novo locus discovery layer. Use it when you do not yet know which genomic intervals matter and you want a deterministic tiled scan over pileup-backed inputs.

## When To Use Discovery vs Contrasts

- Use `region_discovery` when you need to find candidate loci first.
- Use `region_contrasts` when you already know the regions and want to test them formally.
- If you already have a BED file or a matched region set, skip discovery and go straight to `region_contrasts`.

## Minimal Example

```python
from dimelo import region_discovery
from dimelo.models import SampleSpec

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

- For formal region testing, convert discovered hits into BED and pass them into `region_contrasts.score_regions(...)`.
- For downstream clustering, use the discovered loci to focus the region set before running the clustering workflow.
- Keep discovery and contrast roles separate: discovery finds candidates, contrasts tests known regions, and clustering explains state mixtures.

