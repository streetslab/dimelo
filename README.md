# Install Instructions
## Local Install via Conda

### Prerequisites

Platforms: Mac and Linux operating systems, ARM (e.g. M1/M2 mac) and x86 (e.g. Intel mac) architectures. 

**For Windows,  we recommend using Google Colab.** 

*The only option to install on Windows locally is currently to build modkit yourself from the rust [source code](https://github.com/nanoporetech/modkit), or find a pre-built executable somewhere. Without modkit, you can run the plotting functionality of the package but not the bam_parsing functionality.*

Conda installation: https://www.anaconda.com/download

### Load source code from the modkit_parsing_main branch

Open your terminal or command line and navigate to wherever you want to keep the `dimelo` source code (e.g. your Documents folder, `cd Documents`) and clone the repo

```
git clone -b modkit_parsing_beta https://github.com/streetslab/dimelo
```

### Set up virtual environment

Navigate into the dimelo directory

```
cd dimelo
```

Create a conda environment using environment.yml. This will make a new conda environment with the name `dimelo_modkit_parsing` and python 3.11 and nanportech::modkit.

```
conda env create -f environment.yml
```

Alternatively, you can install modkit into any conda environment you like, including the base environment. If you want to, you can install modkit some other way, and then add it to the path of your notebook or script. *NOTE: if you are creating the environment yourself, be sure to use python 3.10 or greater. Some dimelo package features require relatively new python releases.*

```
conda install nanoporetech::modkit==0.2.4
OR
# install modkit some other way
# add to path in python
import sys
sys.path.append('path_to_modkit_executable_directory')
```

### Install pip dependencies and core dimelo package

Activate your conda environment, which should now contain python 3.11 and a modkit executable on the path and executable on your system.

```
conda activate dimelo_modkit_parsing
```

Ensure that you are still in the top-level dimelo directory, and navigate back to it if you are not. You can check by running `ls`; you should see the following files and directories:

```
dimelo
README.md
environment.yml
setup.py
```

Install the dimelo package and its dependencies from source.

```
pip install .
```

## Google Colab

Run the following code in the first cell of your notebook to grab `modkit v0.2.4` from conda and install the `dimelo modkit_parsing_main` branch. This will have to be run whenever you make a new Colab instance, unless you have a better way of managing this, in which case please reach out.

```
from google.colab import drive
drive.mount('/content/drive')
!pip install -q condacolab
import condacolab
condacolab.install()
!conda install nanoporetech::modkit==0.2.4
!git clone -b modkit_parsing_main https://github.com/streetslab/dimelo
!cd dimelo && pip install .
import dimelo
```

# Basic Use

See the [tutorial](tutorial.ipynb) as a starting point. Interface design is not finalized at this time (Jan 23, 2024).

If you want to run the tutorial on Google colab, you can download [tutorial.ipynb](tutorial.ipynb), upload it to your drive, and follow the instructions in the cells.

For local operation on Mac or Linux, you will already have cloned the repo to disk in the installation step. Activate your conda environment, make sure you have jupyter installed, and then launch a jupyter notebook server and navigate to `tutorial.ipynb`. You can also use other tools to open the jupyter notebook or you can simply reference it as an example.

```
conda activate dimelo_modkit_parsing
pip install jupyter
jupyter notebook
```

## Parsing and processing

Both pileup and extract are typically run with a .bed file of regions, which can then be also passed to the plotting functions. One can also pass regions as a string, in the format `chrX:XXX-XXX`, and pass more than one .bed file or string in a list. All regions are processed into a file called `regions.processed.bed` which follows the format required by `modkit`:
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
    cleanup: bool = True,) -> Path:
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
    cleanup: bool = True,) -> Path:
```

For human-readable pileups (bedmethyl files, .bed) and extracted reads (.txt tab-separated values), run with `cleanup=False`. `cleanup=True` will clear these outputs because they can take up a lot of space.

## Plotting

`plot_enrichment_profile` module for pileup line plot profiles across one or more region
```
def plot_enrichment_profile(mod_file_names: list[str | Path],
                            regions_list: list[str | Path | list[str | Path]],
                            motifs: list[str],
                            sample_names: list[str],
                            window_size: int,
                            smooth_window: int | None = None,
                            **kwargs) -> Axes:
def by_modification(mod_file_name: str | Path,
                    regions: str | Path,
                    motifs: list[str],
                    *args,
                    **kwargs) -> Axes:
def by_regions(mod_file_name: str | Path,
               regions_list: list[str | Path | list[str | Path]],
               motif: str,
               sample_names: list[str] = None,
               *args,
               **kwargs) -> Axes:
def by_dataset(mod_file_names: list[str | Path],
               regions: str | Path | list[str | Path],
               motif: str,
               sample_names: list[str] = None,
               *args,
               **kwargs) -> Axes:
```
`plot_enrichment` module for enrichment (e.g. mA/A) bar plot comparisons
```

def plot_enrichment(mod_file_names: list[str | Path],
                    regions_list: list[str | Path | list[str | Path]],
                    motifs: list[str],
                    sample_names: list[str],
                    **kwargs) -> Axes:
def by_modification(mod_file_name: str | Path,
                    regions: str | Path | list[str | Path],
                    motifs: list[str],
                    *args,
                    **kwargs) -> Axes:
def by_regions(mod_file_name: str | Path,
               regions_list: list[str | Path | list[str | Path]],
               motif: str,
               sample_names: list[str] = None,
               *args,
               **kwargs) -> Axes:
def by_dataset(mod_file_names: list[str | Path],
               regions: str | Path | list[str | Path],
               motif: str,
               sample_names: list[str] = None,
               *args,
               **kwargs) -> Axes:
```
`plot_reads` module for single read plots
```
def plot_reads(mod_file_name: str | Path,
               regions: str | Path | list[str | Path],
               motifs: list[str],
               window_size: int = None,
               sort_by: str | list[str] = 'shuffle',
               thresh: float = None,
               relative: bool = True,
               **kwargs
              ) -> Axes:
```
## Load values from processed files

`load_processed.pileup_counts_from_bedmethyl` for valid/modified counts from a specified region or set of regions
```
def pileup_counts_from_bedmethyl(bedmethyl_file: Path,
                          motif: str,
                          regions: str | Path | list[str | Path] = None,
                          ) -> tuple[int, int]:
```
`load_processed.pileup_vectors_from_bedmethyl` for valid over modified fraction from a specified region or set of regions
```
def pileup_vectors_from_bedmethyl(bedmethyl_file: str | Path,
                          motif: str,
                          regions: str | Path | list[str | Path],
                          window_size: int = None) -> (np.ndarray,np.ndarray):
```
`load_processed.read_vectors_from_hdf5` for read-by-basemod lists of valid and modified positions
```
def read_vectors_from_hdf5(
        file: str | Path,
        motifs: list[str],
        regions: str | Path | list[str | Path] = None,
        window_size: int = None,
        sort_by: str | list[str] = ['chromosome','region_start','read_start'],
        calculate_mod_fractions: bool = True,
) -> (list[tuple],list[str],dict):
```
