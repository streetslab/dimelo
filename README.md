# Install Instructions
## Local Install via Conda

### Prerequisites

Platforms: Mac and Linux operating systems, ARM (e.g. M1/M2 mac) and x86 (e.g. Intel mac) architectures. 

**For Windows,  we recommend using Google Colab.** 

*The only option to install on Windows locally is currently to build modkit yourself from the rust [source code](https://github.com/nanoporetech/modkit), or find a pre-built executable somewhere. Without modkit, you can run the plotting functionality of the package but not the bam_parsing functionality.*

Conda installation: https://www.anaconda.com/download

### Load source code from the modkit_parsing_main branch

Clone the repo

```
git clone -b modkit_parsing_main https://github.com/streetslab/dimelo
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

Run the following code in the first cell of your notebook to grab `modkit v0.2.4` from conda and install the `dimelo modkit_parsing_main` branch.

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

## Parsing and processing

`parse_bam.pileup` for profiles and enrichment between regions/modifications

`parse_bam.extract` for single read plots

## Plotting

`plot_enrichment_profile` module for pileup line plot profiles across one or more region

`plot_enrichment` module for enrichment (e.g. mA/A) bar plot comparisons

`plot_reads` module for single read plots

## Load values from processed files

`load_processed.counts_from_bedmethyl` for valid/modified counts from a specified region or set of regions

`load_processed.vector_from_bedmethyl` for valid over modified fraction from a specified region or set of regions

`load_processed.reads_from_hdf5` for read-by-basemod lists of valid and modified positions
