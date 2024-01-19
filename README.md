# Install Instructions
## Conda (Mac/Linux, ARM/x86)

NOTE: For Windows, the only option is currently to build modkit yourself from the rust source code, or find a pre-built executable somewhere. Without modkit, you can run the plotting functionality of the package but not the bam_parsing functionality.

Clone the repo

```
git clone -b modkit_parsing_main https://github.com/streetslab/dimelo
```

Navigate into the dimelo directory

```
cd dimelo
```

Create a conda environment using environment.yml

```
conda env create -f environment.yml
```

Alternatively, you can install modkit into any conda environment you like, including the base environment. If you want to, you can install modkit some other way, and then add it to the path of your notebook or script.

```
conda install nanoporetech::modkit==0.2.4
OR
# install modkit some other way
# add to path in python
import sys
sys.path.append('path_to_modkit_executable_directory')
```

Install the dimelo package and its dependencies from source (you should once again be inside the dimelo cloned repo folder for this to work

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
