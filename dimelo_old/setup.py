# The setup.py file is used as the build script for setuptools. Setuptools is a
# package that allows you to easily build and distribute Python distributions.

import os

import setuptools

from version import version as this_version

# write version to dimelo directory so it can be accessed by the command line
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(
    os.path.join(this_directory, "dimelo", "_version.py"), "wt"
) as fversion:
    fversion.write(f'__version__ = "{this_version}"\n')

# Define required packages. Alternatively, these could be defined in a separate
# file and read in here.
REQUIRED_PACKAGES = [
    "matplotlib==3.4.2",
    "pytest==6.2.4",
    "pandas==1.3.2",
    "numpy==1.20.3",
    "joblib==1.0.1",
    "plotly==4.14.3",
    "seaborn==0.11.2",
    "pyranges==0.0.104",
    "psutil==5.8.0",
    "pybedtools==0.8.1",
    "kaleido==0.2.1",
    "sorted-nearest==0.0.33",  # pin to 0.0.33 to avoid import error (https://github.com/pyranges/sorted_nearest/issues/5)
    "sphinx==4.4.0",
    "tqdm==4.64.0",
]

# Read in the project description. We define this in the README file.
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="dimelo",  # name of project
    install_requires=REQUIRED_PACKAGES,  # all requirements used by this package
    version=this_version,  # project version, read from version.py
    author="Annie Maslan & Reet Mishra & Jeremy Marcus",  # Author, shown on PyPI
    scripts=[
        "bin/dimelo-parse-bam",
        "bin/dimelo-plot-browser",
        "bin/dimelo-plot-enrichment",
        "bin/dimelo-plot-enrichment-profile",
        "bin/dimelo-qc-report",
    ],  # command line scripts installed
    author_email="amaslan@berkeley.edu",  # Author email
    description="Tool for analyzing modified bases from BAM files",  # Short description of project
    long_description=long_description,  # Long description, shown on PyPI
    long_description_content_type="text/markdown",  # Content type. Here, we used a markdown file.
    url="https://github.com/streetslab/dimelo",  # github path
    packages=setuptools.find_packages(),  # automatically finds packages in the current directory. You can also explictly list them.
    classifiers=[  # Classifiers give pip metadata about your project. See https://pypi.org/classifiers/ for a list of available classifiers.
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",  # python version requirement
)
