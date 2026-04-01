# dimelo-toolkit
## Introduction

`dimelo-toolkit` package provides an integrated pipeline for the analysis of multimodal single molecule epigenetic measurements. It is designed with long-read sequencing measurements in mind, but is compatible with any sequencing processing pipeline that generates specification-compliant modbam files. If you need to cite `dimelo-toolkit`, please reference https://www.biorxiv.org/content/10.1101/2025.11.09.687458v1.full. 

`dimelo-toolkit v0.2.x` supersedes the previous tool for processing DiMeLo-seq data, `dimelo v0.1.0`, while supporting numerous new use cases.

<ol>
    <li>Support multicolor data / any base modification context (GpC, CpC, etc) for any properly formatted .bam file, including PacBio sequencing</li>
    <li>Vector extraction for all data types</li>
    <li>Enhanced speed and reliability, sometimes 100-1000x faster, enabling e.g. whole genome processing</li>
    <li>Maintainability -> using a small number of standard dependencies, outsourcing as much as possible to well-maintained third-party packages (e.g. modkit, pysam, h5py, and a few others)</li>
    <li>Modularity in both architecture and operation</li>
    <li>Ease of use, especially for multiplatform installation</li>
    <li>More powerful plotting e.g. bam files from different basecallers, single read sorting, rapid iteration</li>
</ol>

This README document contains installation instructions and documentation for various use cases. There is a [tutorial](#basic-use) jupyter notebook that will take you through the core functionality of the package step-by-step. For Google Colab, the notebook already contains the necessary code to set up and run `dimelo-toolkit`, whereas for local operation you will first need to follow the [local install instructions](#local-install-via-conda). Be sure to check that your system meets our [specifications](#system-requirements). The software is still in early release, and as such, user feedback and requests are greatly appreciated.

## Contents
[1.0 Install Instructions](#Install-instructions)

-[1.1 Local Install via Conda](#Local-Install-via-Conda)
    
-[1.2 Google Colab Installation](#Google-Colab-Installation)

-[1.3 Alternative Installations](#Alternative-Installations)

-[1.4 Developer Installation](#Developer-Installation)

[2.0 Basic Use](#Basic-Use)

-[2.1 Parameters and what they mean](#Parameters-and-what-they-mean)
    
-[2.2 Parsing and processing](#Parsing-and-processing)

-[2.3 Plotting](#Plotting)

-[2.4 Load values from processed files](#load-values-from-processed-files)

-[2.5 Clustering scaffolding](#clustering-scaffolding)

-[2.6 Region contrasts](#region-contrasts)

[3.0 Known Issues](#known-issues)

-[3.1 No progress bars](#no-progress-bars)

# Install Instructions
## Local Install via Conda

### System Requirements

**System Specs:** You will need at least 10GB of disk space and 10GB of RAM for the tutorial (largely due to the reference genome). More disk space may be required if processing large datasets. Additionally, if you want to run on many cores, you should have at least 4GB of RAM per core. If you have less per-core memory than this, consider specifying a subset of cores when calling parsing methods. See the [parameters](#parameters-and-what-they-mean) section for more information. `cores=1` will use the least memory and thus is the least likely to be terminated by your OS.

**Platforms:** Mac and Linux operating systems, ARM (e.g. M1/M2 mac) and x86 (e.g. Intel mac) architectures. The package has been tested on HPC clusters, but there may be additional complexities depending on how these systems are set up.

*For Windows,  we recommend using [Google Colab](https://colab.research.google.com/). We have not tested on [Windows Linux Subsystem](https://learn.microsoft.com/en-us/windows/wsl/install) but in principle that should work too. Windows support is possible in future, but blocked by [conda availability for modkit executables](https://anaconda.org/nanoporetech/modkit) and the [current implementation](dimelo-toolkit/run_modkit.py) of live error/progress tracking during modkit execution, which relies on a unix-only library as of Python 3.11. The urgency of a Windows implementation will depend on user need, so please let us know if this is important for you.*

**Conda and Python:** The default installation requires conda, or alternatives like mamba. See [here](https://www.anaconda.com/download) for conda installation. The installation instructions below will install Python 3.11 for you within a conda virtual environment, but depending on your system configuration you may need to ensure that you are not also loading a different version of Python on your path. If you encounter unexpected errors when importing `dimelo`, e.g. complaining about syntax, consider checking your Python version.

### Load source code for dimelo-toolkit

Open your terminal or command line and navigate to wherever you want to keep the `dimelo-toolkit` source code (e.g. your Documents folder, `cd Documents`) and clone the repo

```
git clone https://github.com/streetslab/dimelo-toolkit
```

### Set up virtual environment

Navigate into the dimelo-toolkit directory

```
cd dimelo-toolkit
```

Create a conda environment using environment.yml. This will make a new conda environment with the name `dimelo-toolkit`. 

```
conda env create -f environment.yml
```

*If you want to handle environment creation yourself, see [the alternative installation instructions](#alternative-installations).*

### Install pip dependencies and core dimelo-toolkit package

Activate your conda environment, which should now contain python 3.11 and a modkit executable on the path and executable on your system.

```
conda activate dimelo-toolkit
```

Ensure that you are still in the top-level dimelo-toolkit directory. Install the dimelo-toolkit package and its dependencies from source.

```
pip install .
```

## Google Colab Installation

Run the following code in the first cell of your notebook to grab `modkit v0.2.4` from conda and install the `dimelo-toolkit main` branch. This will have to be run whenever you make a new Colab instance, unless you have a better way of managing this, in which case please reach out. The tutorial notebook runs equivalent code blocks to set up your environment, so if you are trying to run the tutorial you can skip to [Basic Use](#basic-use).

```
from google.colab import drive
drive.mount('/content/drive')
!pip install -q condacolab
import condacolab
condacolab.install()
!conda install nanoporetech::modkit==0.2.4
!git clone https://github.com/streetslab/dimelo-toolkit
!cd dimelo-toolkit && pip install ipywidgets==7.7.1 .
import dimelo
```

## Alternative Installations

Alternatively, you can install modkit into any conda environment you like. If you want to, you can install modkit some other way, and then add it to the path of your notebook or script. *NOTE: if you are creating the environment yourself, be sure to use python 3.10 or greater. Some dimelo-toolkit features require relatively new python releases.*

```
conda install nanoporetech::modkit==0.2.4
```
OR
```
# install modkit some other way
# add to path in python before importing dimelo
import sys
sys.path.append('path_to_modkit_executable_directory')
```

## Developer Installation
If you are planning on developing for the `dimelo-toolkit` package, change the `pip install` command to install the package in "editable" mode, so that your code changes are reflected in your environment:
```
pip install . -e
```

Additionally, be aware that this package uses [ruff](https://docs.astral.sh/ruff/) to enforce code standards. To make it easy to check that your changes meet standards, we provide [pre-commit](https://pre-commit.com/) hooks that run the checks automatically when you commit.

After installing [pre-commit](https://pre-commit.com/) on your system, run the following from the top level of the repository to set up the hooks:
```
pre-commit install
```

If you need to manually trigger a formatting check, the following command will forcibly run all checks on the entire repository:
```
pre-commit run --all-files
```

To run functionality tests on your machine, ensure that `pytest` is installed in your conda environment:
```
conda install pytest
```

The tests can be run from the top level of the repository using the following command:
```
pytest
```

# Basic Use

See the [tutorial](tutorial.ipynb) as a starting point.

For local operation on Mac or Linux, you will already have cloned the repo to disk in the installation step. Activate your conda environment, make sure you have jupyter installed, and then launch a jupyter notebook server and navigate to `tutorial.ipynb`. You can also use other tools to open the jupyter notebook or you can simply reference it as an example.

```
conda activate dimelo-toolkit
jupyter notebook
```

If you want to run the tutorial on Google Colab, you can download [tutorial.ipynb](tutorial.ipynb), upload it to [Google Colab](https://colab.research.google.com/), and follow the instructions in the cells.

## Parsing and processing

The general workflow of this package is as follows:
```
Parsing: aligned modbam file (latest .bam spec) --> processed file
Loading: processed file --> python objects
Plotting: python objects --> visualizations
```

Both pileup and extract are typically run with a .bed file of regions, which can then be also passed to the plotting functions. All regions are processed into a file called `regions.processed.bed` which follows the format required by `modkit`:
```
chr14	44123158	44123308	+  .  .
```

`parse_bam.pileup` creates a bedmethyl genome-position-wise pileup for profiles and enrichment plotting between regions/modifications or for pulling out a genomic track at one or more regions.

```
def pileup(
    input_file: str | Path,
    output_name: str,
    ref_genome: str | Path,
    output_directory: str | Path = None,
    regions: str | Path | list[str | Path] = None,
    motifs: list = ['A,0','CG,0'],
    thresh: float = None,
    window_size: int = None,
    cores: int = None,
    log: bool = False,
    cleanup: bool = True,
    quiet: bool = False,
    override_checks: bool = False,) -> Path, Path:
```

`parse_bam.extract` creates an hdf5 file with datasets for different aspects of single read data, which can then be passed to plot single reads.
```
def extract(
    input_file: str | Path,
    output_name: str,
    ref_genome: str | Path,
    output_directory: str | Path = None,
    regions: str | Path | list[str | Path] = None,
    motifs: list = ['A,0','CG,0','GCH,1'],
    thresh: float = None,
    window_size: int = None,
    cores: int = None,
    log: bool = False,
    cleanup: bool = True,
    quiet: bool = False,
    override_checks: bool = False,) -> Path, Path:
```

## Analysis Guides

Three higher-level analysis guides now sit on top of the existing parsing layer:

- [docs/global-analysis.md](docs/global-analysis.md) for pileup-backed global summaries and normalization factors
- [docs/shared-clustering.md](docs/shared-clustering.md) for shared-boundary clustering workflows
- [docs/region-contrasts.md](docs/region-contrasts.md) for known-region motif-abundance contrasts

Use `parse_bam.pileup()` when you want defined-region abundance testing, and use `parse_bam.extract()` when you want read-level clustering or single-read follow-up.

For human-readable pileups (bedmethyl files, .bed) and extracted reads (.txt tab-separated values), run with `cleanup=False`. `cleanup=True` will clear these outputs because they can take up a lot of space.

### Parsing outputs
You should expect to see some text outputs and a series of progress bars. Progress bars tell you an estimated time remaining (typically an overestimate by 2-3x at the beginning of contig/chromosome). If you do not see progress bars, go to the [known issues: no progress bars](#no-progress-bars) section for possible fixes. 

There should not be such issues for command line operation. See below an example of command line progress outputs: you should expect relatively fast pre-processing, 10-90 seconds, and then contig processing times depending heavily on the size of your `.bam` file and the extent of your `regions`.

```
(dimelo-toolkit) % python dimelo_cmd.py
modkit found with expected version 0.2.4
No specified number of cores requested. 8 available on machine, allocating all.
Modification threshold of 0.9 will be treated as coming from range 0-1.
████████████████████| Preprocessing complete for motifs ['A,0'] in chm13.draft_v1.1.fasta:  100% | 00:30
███████████████████| All regions complete in mod_mappings.01.retagged.ma.sorted.bam:  100% | 02:23<00:00
████████████████| processed 218324 reads, 13323144 rows, skipped ~184465 reads, failed ~0 reads:  100%
```

You should see no outputs at all if `quiet=True`.

## Plotting

`plot_enrichment_profile` module for pileup line plot profiles across one or more region
```
def plot_enrichment_profile(
    mod_file_names: list[str | Path],
    regions_list: list[str | Path | list[str | Path]],
    motifs: list[str],
    sample_names: list[str],
    window_size: int,
    relative: bool = True,
    single_strand: bool = False,
    regions_5to3prime: bool = False,
    smooth_window: int | None = None,
    quiet: bool = False,
    cores: int | None = None,
    **kwargs,
) -> Axes:
    """
    Plot enrichment profiles, overlaying the resulting traces on top of each other.

    Each input list is expected to be parallel and the same length. Each index represents one analysis condition across the lists.
    Using the same file for multiple conditions requires adding the same file multiple times, in the appropriate indices.

    This is the most flexible method for enrichment profile plotting. For most use cases, consider
    using one of the plot_enrichment_profile.by_* methods.

    Args:
        mod_file_names: list of paths to modified base data files
        bed_file_names: list of paths to bed files specifying centered equal-length regions
        mod_names: list of modifications to extract; expected to match mods available in the relevant mod_files
        sample_names: list of names to use for labeling traces in the output; legend entries
        window_size: half-size of the desired window to plot; how far the window stretches on either side of the center point
        relative: True means x-axis is centered around region centers, False means x-axis is absolute genome positions. Must be True when plotting more than one region.
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        regions_5to3prime: True means negative strand regions get flipped, False means no flipping
        smooth_window: size of the moving window to use for smoothing. If set to None, no smoothing is performed
        quiet: disables progress bars
        cores: CPU cores across which to parallelize processing. Default to None, which means all available.
        kwargs: other keyword parameters passed through to utils.line_plot

    Returns:
        Axes object containing the plot
    """
def by_modification(
    mod_file_name: str | Path,
    regions: str | Path,
    motifs: list[str],
    **kwargs,
) -> Axes:
    """
    Plot enrichment profile, holding modification file and regions constant, varying modification types

    See plot_enrichment_profile for details.
    """
def by_regions(
    mod_file_name: str | Path,
    regions_list: list[str | Path | list[str | Path]],
    motif: str,
    sample_names: list[str] | None = None,
    **kwargs,
) -> Axes:
    """
    Plot enrichment profile, holding modification file and modification types constant, varying regions

    Note: Sample names default to the names of the bed files.

    See plot_enrichment_profile for details.
    """
def by_dataset(
    mod_file_names: list[str | Path],
    regions: str | Path | list[str | Path],
    motif: str,
    sample_names: list[str] | None = None,
    **kwargs,
) -> Axes:
    """
    Plot enrichment profile, holding modification types and regions constant, varying modification files

    Note: Sample names default to the names of the modification files.

    See plot_enrichment_profile for details.
    """
```
`plot_enrichment` module for enrichment (e.g. mA/A) bar plot comparisons
```
def plot_enrichment(
    mod_file_names: list[str | Path],
    regions_list: list[str | Path | list[str | Path]],
    motifs: list[str],
    sample_names: list[str],
    window_size: int | None = None,
    single_strand: bool = False,
    quiet: bool = False,
    cores: int | None = None,
    **kwargs,
) -> Axes:
    """
    Plot enrichment comparison barplots using the given list of pre-processed input files.

    Each input list is expected to be parallel and the same length. Each index represents one analysis condition across the lists.
    Using the same file for multiple conditions requires adding the same file multiple times, in the appropriate indices.

    This is the most flexible method for enrichment plotting. For most use cases, consider
    using one of the plot_enrichment.by_* methods.

    Args:
        mod_file_names: list of paths to modified base pileup data files
        bed_file_names: list of paths to bed files specifying regions to extract
        mod_names: list of modifications to extract; expected to match mods available in the relevant mod_files
        sample_names: list of names to use for labeling bars in the output; x-axis labels
        window_size: (currently disabled) window around center of region, +-window_size//2
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        quiet: disables progress bars
        cores: CPU cores across which to parallelize processing. Default to None, which means all available.
        kwargs: other keyword parameters passed through to utils.bar_plot

    Returns:
        Axes object containing the plot
    """
def by_modification(
    mod_file_name: str | Path,
    regions: str | Path | list[str | Path],
    motifs: list[str],
    **kwargs,
) -> Axes:
    """
    Plot enrichment bar plots, holding modification file and regions constant, varying modification types

    See plot_enrichment for details.
    """
def by_regions(
    mod_file_name: str | Path,
    regions_list: list[str | Path | list[str | Path]],
    motif: str,
    sample_names: list[str] | None = None,
    **kwargs,
) -> Axes:
    """
    Plot enrichment bar plots, holding modification file and modification types constant, varying regions

    Note: Sample names default to the names of the bed files.

    See plot_enrichment for details.
    """
def by_dataset(
    mod_file_names: list[str | Path],
    regions: str | Path | list[str | Path],
    motif: str,
    sample_names: list[str] | None = None,
    **kwargs,
) -> Axes:
    """
    Plot enrichment bar plots, holding modification types and regions constant, varying modification files

    Note: Sample names default to the names of the modification files.

    See plot_enrichment for details.
    """
```
`plot_reads` module for single read plots
```
def plot_reads(
    mod_file_name: str | Path,
    regions: str | Path | list[str | Path],
    motifs: list[str],
    window_size: int | None = None,
    single_strand: bool = False,
    regions_5to3prime: bool = False,
    sort_by: str | list[str] = "shuffle",
    thresh: float | None = None,
    relative: bool = True,
    **kwargs,
) -> Axes:
    """
    Plots centered single reads as a scatterplot, cut off at the boundaries of the requested regions

    Args:
        mod_file_name: path to file containing modification data for single reads
        regions: path to bed file specifying regions to extract
        motifs: list of modifications to extract; expected to match mods available in the relevant mod_files
        window_size: we plot +-window_size//2 from the center of the region(s)
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        regions_5to3prime: True means negative strand regions get flipped, False means no flipping. Only works if relative=True
        sort_by: ordered list for hierarchical sort. Currently only smallest to biggest.
        thresh: if no threshold has been applied already, this will threshold the mod calls for plotting (method is only boolean)
        relative: if True, all regions are centered

    Returns:
        Axes object containing the plot
    """
```
`plot_read_browser` module for single read browser plots
```
def plot_read_browser(
    mod_file_name: str | Path,
    region: str,
    motifs: list[str],
    thresh: int | float | None = None,
    single_strand: bool = False,
    sort_by: str | list[str] = "shuffle",
    hover: bool = True,
    subset_parameters: dict | None = None,
    span_full_window: bool = False,
    width: int = 1,
    size: int = 4,
    colorscales: dict = utils.DEFAULT_COLORSCALES,
    **kwargs,
) -> plotly.graph_objs.Figure:
    """
    Plot base modifications on single reads in a high-quality, interactive-enabled fashion.

    This method returns a plotly Figure object, which can be used in a number of ways to view and save
    the figure in different formats. To view the figure interactively (in a notebook or python script),
    simply call the show() method of the returned Figure object. See the helper methods below for saving
    figures.

    Additional keyword arguments will be passed down to collapse_rows() if sort_by == "collapse". See that
    method for details.

    Args:
        mod_file_name: path to file containing modification data for single reads
        region: region string specifying the region to plot
        motifs: list of modifications to extract; expected to match mods available in the relevant mod_files
        thresh: A modification calling threshold. While the browser always displays float probabilities, setting
            this to a value will limit display to only modification events over the given threshold. Else, display
            all modifications regardless of probability.
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        sort_by: ordered list for hierarchical sort; see load_processed.read_vectors_from_hdf5() for details.
            Can also pass the argument "collapse" to allow multiple reads on single rows of the browser, for a
            more condensed visualization. Note that "collapse" is mutually exclusive with all other sorting options,
            and is only allowed to be passed as a single string option.
        hover: if False, disables display of information on mouse hover
        subset_parameters: Parameters to pass to the utils.random_sample() method, to subset the
            reads to be returned. If not None, at least one of n or frac must be provided.
        span_full_window: if True, only plot reads that fully span the window defined by region_start-region_end
        width: width of the read lines in the browser
        size: size of the modification event markers in the browser
        colorscales: dictionary mapping motif names to plotly colorscale specifications

    Returns:
        plotly Figure object containing the plot
    """
```
`plot_depth_profile` module for motif modification-tagged read depth profiles
```
def plot_depth_profile(
    mod_file_names: list[str | Path],
    regions_list: list[str | Path | list[str | Path]],
    motifs: list[str],
    sample_names: list[str],
    window_size: int | None = None,
    single_strand: bool = False,
    regions_5to3prime: bool = False,
    smooth_window: int | None = None,
    quiet: bool = False,
    cores: int | None = None,
    **kwargs,
) -> Axes:
    """
    Plot depth profiles, overlaying the resulting traces on top of each other.

    Each input list is expected to be parallel and the same length. Each index represents one analysis condition across the lists.
    Using the same file for multiple conditions requires adding the same file multiple times, in the appropriate indices.

    This is the most flexible method for depth profile plotting. For most use cases, consider
    using one of the plot_depth_profile.by_* methods.

    Args:
        mod_file_names: list of paths to modified base data files
        bed_file_names: list of paths to bed files specifying centered equal-length regions
        mod_names: list of modifications to extract; expected to match mods available in the relevant mod_files
        sample_names: list of names to use for labeling traces in the output; legend entries
        window_size: half-size of the desired window to plot; how far the window stretches on either side of the center point
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        regions_5to3prime: True means negative strand regions get flipped, False means no flipping
        smooth_window: size of the moving window to use for smoothing. If set to None, no smoothing is performed
        quiet: disables progress bars
        cores: CPU cores across which to parallelize processing. Default to None, which means all available.
        kwargs: other keyword parameters passed through to utils.line_plot

    Returns:
        Axes object containing the plot
    """
def by_modification(
    mod_file_name: str | Path,
    regions: str | Path,
    motifs: list[str],
    **kwargs,
) -> Axes:
    """
    Plot depth profile, holding modification file and regions constant, varying modification types

    See plot_depth_profile for details.
    """
def by_regions(
    mod_file_name: str | Path,
    regions_list: list[str | Path | list[str | Path]],
    motif: str,
    sample_names: list[str] | None = None,
    **kwargs,
) -> Axes:
    """
    Plot depth profile, holding modification file and modification types constant, varying regions

    Note: Sample names default to the names of the bed files.

    See plot_depth_profile for details.
    """
def by_dataset(
    mod_file_names: list[str | Path],
    regions: str | Path | list[str | Path],
    motif: str,
    sample_names: list[str] | None = None,
    **kwargs,
) -> Axes:
    """
    Plot depth profile, holding modification types and regions constant, varying modification files

    Note: Sample names default to the names of the modification files.

    See plot_depth_profile for details.
    """
```
`plot_depth_histogram` module for distribution of motif modification-tagged reads
```
def plot_depth_histogram(
    mod_file_names: list[str | Path],
    regions_list: list[str | Path | list[str | Path]],
    motifs: list[str],
    sample_names: list[str],
    window_size: int | None = None,
    single_strand: bool = False,
    one_depth_per_region: bool = False,
    quiet: bool = False,
    cores: int | None = None,
    split_large_regions: bool = False,
    **kwargs,
) -> Axes:
    """
    Plot depth histograms, overlaying the results on top of each other.

    Each input list is expected to be parallel and the same length. Each index represents one analysis condition across the lists.
    Using the same file for multiple conditions requires adding the same file multiple times, in the appropriate indices.

    This is the most flexible method for depth histogram plotting. For most use cases, consider
    using one of the plot_depth_histogram.by_* methods.

    Args:
        mod_file_names: list of paths to modified base data files
        bed_file_names: list of paths to bed files specifying centered equal-length regions
        mod_names: list of modifications to extract; expected to match mods available in the relevant mod_files
        sample_names: list of names to use for labeling traces in the output; legend entries
        window_size: half-size of the desired window to plot; how far the window stretches on either side of the center point
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        one_depth_per_region: if True, each region will only report a single depth value, averaging across all non-zero depths. If False
            depths will be reported separately for all nonzero count positions in each region for a more granular view of depth distribution.
        quiet: disables progress bars
        cores: CPU cores across which to parallelize processing. Default to None, which means all available.
        split_large_regions: if True, regions will be run sequentially in parallelized chunks. If False,
            each individual region's chunks will be run sequentially but there will be parallelization across
            regions, i.e. each core will be assigned one region at a time by the executor. Set to True if you
            are running a small number of very large regions (e.g. one or two chromosomes), otherwise to to False (default).
        kwargs: other keyword parameters passed through to utils.line_plot

    Returns:
        Axes object containing the plot
    """
def by_modification(
    mod_file_name: str | Path,
    regions: str | Path,
    motifs: list[str],
    **kwargs,
) -> Axes:
    """
    Plot depth histogram, holding modification file and regions constant, varying modification types

    See plot_depth_histogram for details.
    """
def by_regions(
    mod_file_name: str | Path,
    regions_list: list[str | Path | list[str | Path]],
    motif: str,
    sample_names: list[str] | None = None,
    **kwargs,
) -> Axes:
    """
    Plot depth histogram, holding modification file and modification types constant, varying regions

    Note: Sample names default to the names of the bed files.

    See plot_depth_histogram for details.
    """
def by_dataset(
    mod_file_names: list[str | Path],
    regions: str | Path | list[str | Path],
    motif: str,
    sample_names: list[str] | None = None,
    **kwargs,
) -> Axes:
    """
    Plot depth histogram, holding modification types and regions constant, varying modification files

    Note: Sample names default to the names of the modification files.

    See plot_depth_histogram for details.
    """
```
## Load values from processed files

`load_processed.pileup_counts_from_bedmethyl` for valid/modified counts from a specified region or set of regions
```
def pileup_counts_from_bedmethyl(
    bedmethyl_file: str | Path,
    motif: str,
    regions: str | Path | list[str | Path],
    window_size: int | None = None,
    single_strand: bool = False,
    quiet: bool = False,
    cores: int | None = None,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> tuple[int, int]:
    """
    User-facing function.

    Extract number of modified bases and total number of bases from the given bedmethyl file.
    Called by plotters or by the user.

    This function loops through all the provided regions and pulls those regions up in the input
    sorted and indexed bedmethyl file. For rows within those regions, checks that the motif
    is correct (i.e. sequence context, modified base, mod code, and optionally strand). All
    correct locations are included in the sum counts that get returned.

    If no regions are specified, returns the sum total for the motif of interest across the
    entire bedmethyl file.

    Args:
        bedmethyl_file: Path to bedmethyl file
        regions: Path to bed file specifying regions
        motif: type of modification to extract data for
        window_size: (currently disabled) window around center of region, +-window_size
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        quiet: disables progress bars
        cores: CPU cores across which to parallelize processing. Default to None, which means all available.
        chunk_size: size of genomic subregions to assign out to each process

    Returns:
        tuple containing counts of (modified_bases, total_bases)
    """
```
`load_processed.pileup_vectors_from_bedmethyl` for valid and modified counts-by-position from a specified region or set of length-matched regions
```
def pileup_vectors_from_bedmethyl(
    bedmethyl_file: str | Path,
    motif: str,
    regions: str | Path | list[str | Path],
    window_size: int | None = None,
    single_strand: bool = False,
    regions_5to3prime: bool = False,
    quiet: bool = False,
    cores: int | None = None,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> tuple[np.ndarray, np.ndarray]:
    """
    User-facing function.

    Extract per-position pileup counts at valid motifs across one or more superimposed regions.
    Called by profile plotters, can also be used by a user directly.

    Returns two vectors:
    * Total number of times a modified base in the motif was found at each position
    * Total number of times the motif was found at each position

    This function loops through all the provided regions and fetches those regions from the
    bedmethyl file. For rows within those regions, it checks that the motif
    is correct (i.e. sequence context, modified base, mod code, and optionally strand). It then adds
    to two vectors (mod and valid). By default all regions are assumed to
    be the same size (the size of the first region).

    If regions_5to3prime is set to True, then negative strand regions are flipped to that all regions
    are superimposed along the 5 prime to 3 prime direction, which can be helpful if there is
    directionality to the signal (e.g. upstream v downstream relative to TSSs, TF binding sites, and so on).
    A region must be provided because otherwise there is no way to know what vector to return.
    However, a region can be a whole chromosome if desired.

    Args:
        bedmethyl_file: Path to bedmethyl file
        regions: Path to bed file specifying centered equal-length regions
        motif: type of modification to extract data for
        window_size: the extent in either direction for windows around the center of regions.
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        regions_5to3prime: True means negative strand regions get flipped, False means no flipping
        quiet: disables progress bars
        cores: CPU cores across which to parallelize processing. Default to None, which means all available.
        chunk_size: size of genomic subregions to assign out to each process

    Returns:
        tuple containing (modified_base_counts, valid_base_counts)
    """
```
`load_processed.read_vectors_from_hdf5` for read-by-basemod lists of valid and modified positions
```
def read_vectors_from_hdf5(
    file: str | Path,
    motifs: list[str],
    regions: str | Path | list[str | Path] | None = None,
    window_size: int | None = None,
    single_strand: bool = False,
    sort_by: str | list[str] = ["chromosome", "region_start", "read_start"],
    calculate_mod_fractions: bool = True,
    quiet: bool = True,  # currently unused; change to default False when pbars are implemented
    cores: int | None = None,  # currently unused
    subset_parameters: dict | None = None,
    span_full_window: bool = False,
) -> tuple[list[tuple], list[str], dict | None]:
    """
    User-facing function.

    Pulls a list of read data out of an .h5 file containing processed read vectors, formatted
    for read-by-read vector processing downstream use cases.

    The flow of operation here is we load up the h5 file then loop through our regions and pick
    out reads corresponding to our criteria. Criteria include chromosome, read starts and ends
    (compared to region starts and ends), motif, and strand (if single_strand is True). The indices
    for the desired reads are identified region-by-region, then all the reads for the region (or
    the whole h5, if no region is passed) are loaded using the process_data function and put into
    a list. The bytes are then decoded for the array entries, which are manually compressed because
    h5py wasn't behaving.

    There's some adjustment for the raw probability (no thresh) to match modkit extract outputs.
    Specifically, the 0-255 8bit int has 0.5 added before dividing by 256, giving mod qualities
    between 0.001953 and 0.99805 for bases in valid motifs. (Invalid positions have zeros.)

    After this processing, we calculate modification fractions, sort, and return.

    Args:
        file: Path to an hdf5 (.h5) file containing modification data for single reads,
            stored in datasets read_name, chromosome, read_start,
            read_end, base modification motif, mod_vector, and val_vector.
        regions: Single or list of Path objects or strings. Path objects must point to .bed
            files, strings can be .bed paths or region string in the format chrX:XXX-XXX.
            All should all be regions for which your original .bam file had reads extracted,
            although by design this method will not raise an error if any region contains
            zero reads, as this may simply be a matter of low read depth.
            If no regions are specified, the entire .h5 file will be returned. This may cause
            memory issues.
        motifs: types of modification to extract data for. Motifs are specified as
            {DNA_sequence},{position_of_modification}. For example, a methylated adenine is specified
            as 'A,0' and CpG methylation is specified as 'CG,0'.
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        window_size: An optional parameter for creating centered windows for the provided regions.
            If provided, all regions will be adjusted to be the same size and centered. If not provided,
            all regions should already be the same size, or there should be only one.
        sort_by: Read properties by which to sort, either one string or a list of strings. Options
            include chromosome, region_start, region_end, read_start, read_end, and motif. More to
            be added in future.
        quiet: silences progress bars (currently unused)
        cores: cores across which to parallelize processes (currently unused)
        subset_parameters: Parameters to pass to the utils.random_sample() method, to subset the
            reads to be returned. If not None, at least one of n or frac must be provided. The array
            parameter should not be provided here.
        span_full_window: If True, only load reads that fully span the window defined by region_start-region_end

    Returns:
        a list of tuples, each tuple containing all datasets corresponding to an individual read that
        was within the specified regions.
        a list of strings, naming the datasets returned.
        a regions_dict, containing lists of (region_start,region_end) coordinates by chromosome/contig.
    """
```
`load_processed.regions_to_list` for parallel loading of many loci into a list for easy downstream analysis.
```
def regions_to_list(
    function_handle,
    regions,
    window_size: int | None = None,
    quiet: bool = True,
    cores: int | None = None,
    split_large_regions: bool = False,
    **kwargs,
):
    """
    User-facing function.

    Run any standard load_processed pileup or extract loader loading each region from the region
    specifier into a new element of a list.

    Args:
        function_handle: the loader function you want to run.
        regions: the region specifier. Typically we expect to get many regions for this function, in the form of a list
            of strings or bed file paths. regions_to_list will run across all of these one-by-one returning a separate
            function return for each independent region.
        window_size: window around centers of regions, defaults to None
        quiet: disables progress bars
        cores: CPU cores across which to parallelize processing. Default to None, which means all available.
        split_large_regions: if True, regions will be run sequentially in parallelized chunks. If False,
            each individual region's chunks will be run sequentially but there will be parallelization across
            regions, i.e. each core will be assigned one region at a time by the executor. Set to True if you
            are running a small number of very large regions (e.g. one or two chromosomes), otherwise to to False (default).
        **kwargs: all necessary keyword arguments to pass down to the loader

    Returns:
        List(function_handle return objects per region)
    """
```

## Clustering scaffolding

`dimelo.cluster` wraps the loading utilities above so you can quickly build matrices for downstream clustering or dimensionality reduction.

For the higher-level shared-boundary workflows added on top of these loaders, see [docs/shared-clustering.md](docs/shared-clustering.md). The short version is:

- run `parse_bam.extract()` for `mode="read_global"`
- run `parse_bam.pileup()` for `mode="region_anchored"`
- run both when you want region-level summaries plus read-level follow-up

```python
from dimelo import workflows
from dimelo.models import SampleSpec

result = workflows.shared_cluster_distribution(
    samples=[
        SampleSpec(sample_id="s1", condition="NS", extract_h5="output/s1.extract.h5"),
        SampleSpec(sample_id="s2", condition="15min", extract_h5="output/s2.extract.h5"),
    ],
    mode="read_global",
    motifs=["A,0"],
    n_clusters=8,
)
```

```
from dimelo import cluster

# Region-level features: rows are regions, columns are aligned positions.
pileup_matrix, region_info = cluster.region_feature_matrix_from_pileup(
    bedmethyl_file="output/pileup.sorted.bed.gz",
    motif="A,0",
    regions="regions.bed",
    window_size=500,
    regions_5to3prime=True,
)
labels, model = cluster.cluster_features(pileup_matrix, n_clusters=4)

# Read-level feature pipeline using extract outputs.
read_windows = cluster.extract_read_windows(
    hdf5_file="output/extract.h5",
    motifs=["A,0", "CG,0"],
    regions="regions.bed",
    config=cluster.ReadWindowExtractionConfig(window_size=2000, orientation_aware=True),
)
feature_matrix, feature_names = cluster.read_window_feature_matrix(read_windows)
read_clusters = cluster.cluster_read_windows(feature_matrix, method="kmeans", n_clusters=8)

# Optional visualization/export
cluster.plot_region_cluster_profiles(pileup_matrix, labels, window_bp=1000)
cluster.export_region_clusters_to_bed(region_info, labels, "region_clusters.bed")
cluster.plot_cluster_karyotype("region_clusters.bed", "ref.fasta.fai")
# If pileup_matrix concatenates multiple motifs, pass motif_index to select which slice to visualize.

# Multi-site read raster (paired windows) for exploratory QC
# fig, stats = cluster.plot_multisite_read_raster(
#     read_windows,
#     n_windows=2,
#     min_separation_bp=5000,
#     motif_index=0,
#     smoothing="gaussian",
# )

# Binary classification of read features across two samples
clf_result = cluster.classify_read_features_binary(
    feature_matrix,
    sample_labels=my_sample_labels,  # e.g., ["sampleA", "sampleA", "sampleB", ...]
    classifier="xgboost",           # or logreg/sgd/random_forest/svc/knn
    random_state=42,
)
print(clf_result["metrics"])
# QC plots
cluster.plot_confusion_matrices(clf_result["predictions"])
# If you have per-read windows aligned to predictions for this split:
# cluster.plot_classification_profiles(read_windows.data_matrix, clf_result["predictions"], split="test")

# Optional: sample rows for faster exploratory runs
sampled_features, sampled_labels, idx = cluster.sample_rows(feature_matrix, labels, frac=0.2, stratify=True)
cluster.cluster_read_windows(sampled_features, method="kmeans", n_clusters=4)
```

`region_feature_matrix_from_pileup` delegates to `load_processed.regions_to_list`, so you inherit the same controls for multi-core execution and windowing. For read-level work you can either keep using `read_mod_fraction_table` (motif-wide fractions) or run the new window-based pipeline above: `extract_read_windows` builds a matrix of per-base modification calls (optionally enforcing strand orientation and dropping reads spanning multiple regions), `read_window_feature_matrix` augments that matrix with PCA, autocorrelation, density, and summary statistics, and `cluster_read_windows` exposes the same clustering options demonstrated in the exploratory notebook (k-means, GMM, spectral, OPTICS/HDBSCAN, etc.).

`cluster_features` currently offers a lightweight k-means wrapper powered by `scikit-learn`, while `cluster_read_windows` can switch between `kmeans`, `minibatch_kmeans`, `gmm`, `agglomerative`, `spectral`, `birch`, `dbscan`, `optics`, `hdbscan`, or `umap_kmeans`. These dependencies (`scikit-learn`, `scipy`, `hdbscan`, `umap-learn`) live in the `clustering` extra (`pip install dimelo[clustering]`), or feed the matrices above into your own estimator.

**Motif semantics:** `extract.h5` stores one row per motif per read. When passing multiple motifs, you will see multiple rows per read; use `build_multimotif_read_windows` to group and concatenate motifs per read (set `require_all_motifs=True` to drop reads missing motifs), or filter by `metadata["motif"]` for per-motif analysis. Visualization helpers take a `motif_index` argument to choose which motif slice to display when windows are concatenated.

## Parameters and what they mean

Many of the parsing, loading, and plotting functions share parameters. The common ones and their meanings/defaults are listed below.

### Functional Parameters

`input_file` for parsing functions is a mandatory parameter providing the path to an aligned .bam file with modification calls, a .bam.bai index, and tags following the latest .bam specifications. parse_bam will check whether this .bam file meets the specifications and tell you what to do if it doesn't.

`output_name` for parsing functions is a mandatory string parameter that will be given to the new folder containing the processed outputs.

`ref_genome` for parsing functions is a mandatory parameter providing the path to the reference genome .fasta file to which the `input_file` .bam is aligned.

`output_directory` for parsing functions is an optional parameter specifying a parent directory for your outputs. By default, this will simply be the directory in which you `input_file` resides.

`mod_file_name` and `mod_file_names` for plotting functions are mandatory parameters providing a path to processed/parsed files that are ready for plotting. These paths are returns by the parsing functions but can also be provided manually by the user as a string or Path object. If providing manually, the path should point to a .bed.gz file with an accompanying .bed.gz.tbi index for profile and enrichment plots and to an .h5 file for read plots. `mod_file_name` points to a single file whereas `mod_files_names` is a list of files.

`regions` and `regions_list` are used for specifying subsets of the genome to parse, load, or plot. A `region` is defined as a range of genomic coordinates, and `regions` can refer to any number of `region` specifications. Thus for the `regions` parameter one can pass a single region specified as a string, `chrX:XXX-XXX`, many regions defined in a .bed tab-separated-value file with each line containing at miniimum chromosome, start, and end coordinates (plus optionally a strand, + or -), or a list of strings specifiers or bed files. The entire list will be rolled into a single `regions` set to be passed down for subsequent processing. In the case of regions-wise comparisons in plotting functions, `regions_list` is a *list of regions objects*, where each element of the list is a string, Path, or list of strings or Paths.

`regions_5to3prime` is used by loader/plotter functions that align data from many loci. It ensures that all loaded data is oriented in the right way with respect to functional elements by flipping region data to all orient 5' to 3' based on the strand specified in the region str or bed file line.

`single_strand` limits loaded/plotted data to only the strand specified in the region str or bed file line. Any modification on a given nucleotide has a different set of possible locations on the forward vs reverse strand which can be relevant for many biological phenomena.

`motif` and `motifs` are used to specify what base modifications you are interested in and what their sequence context is for parse, load, and plot functions. A single `motif` is a string containing several canonical bases (using the [IUPAC nucleic acid notation](https://en.wikipedia.org/wiki/Nucleic_acid_notation), e.g. **H** refers to "not a G"), followed by a comma, and then an integer specifying which coordinate in the string is your modified base. For example, 6mA is denoted "A,0" and CpG is denoted "CG,0" whereas GpC *excluding CpGs* is denoted "GCH,1". `motifs` is a list of such strings for functions that can work on multiple base modifications at once.

`thresh` for parsing and some loading/plotting functions refers to a base modification probability threshold, used to transform the the output of most basecalling pipelines into a binary call for any given read position. For parsing pileup calls, this defaults to `None` which allows `modkit` to pick its own threshold based on the data. For other calls, the parameter is mandatory. The normal use is specifying between 0 and 1, but 1-255 is also supported to make the inputs more backwards compatible with old `dimelo` package versions and with examination of the raw .bam file contents. A value between and 255 will simply be converted into a 0-1 probability before being handed down to subsequent processing.

`window_size` for parsing and most loading and plotting functions is a *modification to your regions* that will redefine them to be all the same size, i.e. 2 x window_size, centered around the centers of your original regions. This is important for the parsing and plotting applications that show many genomic regions at once, but should be left blank if you don't want your regions modified. The default is `None` for functions where the parameter is optional. 

`relative` is a boolean input that specifies whether loading and plotting operations adjust coordinates to be relative to some center point or simple plot in absolute genomic coordinates. This is not currently implemented for some methods, but exists as an interface to allow future implementation.

`sort_by` for plot_reads will sort reads by any of `region_start`, `region_end`, `read_name`, `read_start`, `read_end`, `chromosome`, `strand`, and `MOTIF_mod_fraction` for any extracted motif. New sorting options are planned in the future. The default is `shuffle`, which will put the reads in a random order. `sort_by` can be passed as one string or as a list of strings. If a list is passed, the reads will be sorted hierarchically i.e. first by the first list element, then the second, and so on. The exception is that if any of the list elements are `shuffle`, the reads will *first* be shuffled and then sorted by the rest of the elements in order of priority.

`**kwargs` for all plotting functions except for `plot_read_browser`, which uses `plotly` and doesn't provide as convenient of an interface, `**kwargs` get passed down to the underlying `matplotlib` / `seaborn` plotting functions. See `matplotlib` and `seaborn` documentation for more details.

### Execution Parameters  

`cores` allows you to specify that `modkit` uses only a fraction of all the compute resources of your machine, rather than all.

`log` will save the modkit logs for troubleshooting.

`cleanup` will keep the (often large) human-readable outputs that are inefficient for plotting and vector extraction but may be helpful for other use cases.

`quiet` suppressed progress bars and other outputs.

`override_checks` lets you run modkit even if the bam format checking and reference alignment checking are anomalous.

`CHUNK_SIZE` sets the genomic chunks processed in a single step for parallelized operations.

# Known Issues
## No progress bars
The most common culprit for progress bar issues in notebooks (Jupyter or Colab) is an incompatibility between your notebooks interfaces and your `ipywidgets` version. The latest jupyter notebooks or jupyter lab install and the latest ipywidgets should work together, but on Google Colab, VS Code, Open On Demand, and other jupyter interfaces this may not be the case. [setup.py](setup.py) contains details on which versions you can try downgrading to for different platforms. The following code run in your activated conda environment will downgrade `ipywidgets` to your specified version. **Our Colab instructions in the [Colab Installation](#google-colab-installation) section and the [tutorial](tutorial.ipynb) already handle this for you.**

```
pip install ipywidgets==X.XX.X
```
