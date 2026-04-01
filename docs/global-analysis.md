# Global Analysis

`dimelo.global_analysis` provides pileup-backed whole-sample summaries and broad tiled-window summaries on top of `parse_bam.pileup()` outputs.

## What It Covers

- global motif fractions per sample and condition
- tiled window summaries across genome-size inputs
- reusable normalization factors derived from the global fractions

## Example

```python
from dimelo import global_analysis
from dimelo.models import SampleSpec

result = global_analysis.run_global_analysis(
    samples=[
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="output/s1.extract.h5",
            metadata={"pileup_path": "output/s1.pileup.sorted.bed.gz"},
        )
    ],
    motifs=["A,0"],
    genome_sizes={"chr1": 248956422},
    window_size=100000,
    step_size=50000,
)
```

## Canonical Outputs

`run_global_analysis(...)` returns a `GlobalAnalysisResult` with data-first tables:

- `result.summary`
- `result.windows`
- `result.normalization_factors`
- `result.plot_data["global_fraction_bar"]`
- `result.plot_data["window_fraction_table"]`

The wrapper does not couple the workflow to a plotting backend. If you want to render the results yourself, use the returned tables directly.
