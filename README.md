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

Both pileup and extract are typically run with a .bed file of regions, which can then be also passed to the plotting functions. However, for parsing, a list of many .bed files can be passed and all will be processed.

.bed format needs to have chromosome, start coordinate, end coordinate, strand, and two more values. However, the last three can simply be `.`

```
chr14	44123158	44123308	.  .  .
```

`parse_bam.pileup` for profiles and enrichment between regions/modifications

```
def pileup(
    input_file: str | Path,
    output_name: str,
    ref_genome: str | Path,
    output_directory: str | Path = None,
    region_str: str = None,
    bed_files: list[str | Path] = None,
    basemods: list = ['A,0','CG,0','GCH,1'],
    thresh: float = None,
    window_size: int = None,
    cores: int = None,
    log: bool = False,
    cleanup: bool = True,) -> Path:

    """
    Takes a file containing long read sequencing data aligned 
    to a reference genome with modification calls for one or more base/context 
    and creates a pileup. A pileup is a genome-position-wise sum of both reads with
    bases that could have the modification in question and of reads that are in
    fact modified.

    The current implementation of this method uses modkit, a tool built by 
    Nanopore Technologies, along with htslib tools compress and index the output
    bedmethyl file.
```

`parse_bam.extract` for single read plots

```
def extract(
    input_file: str | Path,
    output_name: str,
    ref_genome: str | Path,
    output_directory: str | Path = None,
    region_str: str = None,
    bed_files: list[str | Path] = None,
    basemods: list = ['A,0','CG,0','GCH,1'],
    thresh: float = 0,
    window_size: int = 0,
    cores: int = None,
    log: bool = False,
    cleanup: bool = True,) -> Path:

    """
    Takes a file containing long read sequencing data aligned 
    to a reference genome with modification calls for one or more base/context 
    and pulls out data from each individual read. The intermediate outputs contain
    a plain-text list of all base modifications, split out by type. The compressed
    and indexed output contains vectors of valid and modified positions within each
    read.

    The current implementation of this method uses modkit, a tool built by 
    Nanopore Technologies, along with h5py to build the final output file.

    https://github.com/nanoporetech/modkit/
```

For human-readable pileups (bedmethyl files, .bed) and extracted reads (.txt tab-separated values), run with `cleanup=False`. `cleanup=True` will clear these outputs because they can take up a lot of space.

## Plotting

`plot_enrichment_profile` module for pileup line plot profiles across one or more region

`plot_enrichment` module for enrichment (e.g. mA/A) bar plot comparisons

`plot_reads` module for single read plots

## Load values from processed files

`load_processed.counts_from_bedmethyl` for valid/modified counts from a specified region or set of regions

`load_processed.vector_from_bedmethyl` for valid over modified fraction from a specified region or set of regions

`load_processed.reads_from_hdf5` for read-by-basemod lists of valid and modified positions
