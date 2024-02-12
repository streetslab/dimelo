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
        'ipywidgets==7.6.5', 
        # ipywidgets 7.7.1 or earlier is needed because tqdm.auto / tqdm.notebook progress bars don't work with later versions on Google Colab, as of Feb 11 2024. 
        # this says that was fixed in 2022 but it sure wasn't when I tried today https://github.com/googlecolab/colabtools/issues/3023
        # on Savio I needed to install ipywidgets==7.6.5 because OOD is really old, so I'm just pushing back to that version because it works on Colab too.
        # No reason we need this version for functionality on other platforms though, to my knowledge.
        'tqdm',
    ],
)