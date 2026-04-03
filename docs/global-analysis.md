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

## Plotting Prep Helpers

Use the shared plotting helpers when you want plot-ready global-analysis tables without committing to a specific renderer:

```python
from dimelo import plotting

summary_payload = plotting.prepare_global_analysis_summary_data(
    result=result,
    aggregate_conditions=True,
)

window_payload = plotting.prepare_global_analysis_window_data(
    result=result,
    aggregate_conditions=True,
)
```

- `prepare_global_analysis_summary_data(...)` consumes a `GlobalAnalysisResult` and returns renderer-neutral `sample_summary`, `condition_summary`, `normalization_table`, and `metadata` payloads. The summary payload exposes both sample-level rows and optional condition-level views, while normalization values stay available through `normalization_table`.
- `prepare_global_analysis_window_data(...)` consumes the same `GlobalAnalysisResult` and returns renderer-neutral `window_table`, `condition_window_table`, and `metadata` payloads. Broad-window payloads stay per-contig by default rather than flattening onto a cumulative genome axis.

For defined-region follow-up, pair these global summaries with [region contrasts](region-contrasts.md). For de novo locus finding first, use [region discovery](region-discovery.md).
