from setuptools import find_packages, setup

setup(
    name="dimelo",
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        # scipy/pandas/matplotlib are imported at module scope by modules that
        # dimelo/__init__ pulls in (e.g. region_contrasts imports scipy.stats), so they
        # are required for `import dimelo`, not optional. Previously they only arrived
        # transitively via seaborn (scipy is not even a guaranteed seaborn dependency),
        # which could break a clean install; declare them explicitly.
        "scipy",
        "pandas",
        "matplotlib",
        # scikit-learn backs the core clustering/classification path (KMeans, etc.).
        "scikit-learn",
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
        # Advanced/optional clustering methods that are imported lazily only when used.
        "clustering": [
            "hdbscan",
            "umap-learn",
            "pyranges",
            "xgboost",
        ]
    },
)
