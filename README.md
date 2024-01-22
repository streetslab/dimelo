# Install Instructions
## Local Install via Conda

### Prerequisites

Platforms: Mac and Linux operating systems, ARM (e.g. M1/M2 mac) and x86 (e.g. Intel mac) architectures. *For Windows, the only option is currently to build modkit yourself from the rust source code, or find a pre-built executable somewhere. Without modkit, you can run the plotting functionality of the package but not the bam_parsing functionality.*

Conda installation: https://www.anaconda.com/download

### Load source code from the modkit_parsing_main branch

Clone the repo

```
git clone -b modkit_parsing_main https://github.com/streetslab/dimelo
```

### Download modkit and install python 3.11 in a new virtual environment

Navigate into the dimelo directory

```
cd dimelo
```

Create a conda environment using environment.yml. This will make a new conda environment with the name `dimelo_modkit_parsing`.

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

Run the following code to grab modkit from conda and install the dimelo modkit_parsing branch

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

## Parsing and processing

## Plotting

## Load values from processed files
