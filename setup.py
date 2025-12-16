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
        "tqdm",
        "plotly",
        "kaleido",
    ],
    extras_require={
        "clustering": [
            "scikit-learn",
            "scipy",
            "hdbscan",
            "umap-learn",
            "pyranges",
            "xgboost",
        ]
    },
)
