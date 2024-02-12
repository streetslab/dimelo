from setuptools import setup, find_packages

setup(
    name='dimelo',
    version='1.0.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'seaborn',
        'pysam',
        'h5py',
        'pyBigWig',
        'ipywidgets',
        'tqdm',
    ],
)