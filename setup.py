from setuptools import find_packages, setup

setup(
    name="dimelo",
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "seaborn",
        "pysam",
        "h5py",
        "pyBigWig",
        "notebook",
        "ipykernel",
        "ipywidgets",
        # ipywidgets is finicky on some platforms
        # Local jupyter notebooks: this works best with ipywidgets and jupyter running the latest versions, which is what this file will do.
        # Google Colab: ipywidgets 7.7.1 seems to be necessary for Google Colab needed for tqdm.auto progress bars as of Feb 11 2024.
        #    This says that was fixed in 2022 but it sure wasn't when I tried https://github.com/googlecolab/colabtools/issues/3023.
        #    You can simply downgrade after installing using pip install ., it won't break anything to our knowledge.
        # VS Code: using the VS code plugin to run jupyter, 7.7.1 also works while the latest version does not. Likely that
        #    some intermediate versions work too.
        # On the Berkeley High Performance Computing Cluster, Savio, I needed to install ipywidgets==7.6.5 for jupyter Open On-Demand.
        "tqdm",
        "plotly",
        "kaleido",
    ],
)
