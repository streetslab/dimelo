Installation
====================

Requirements
------------

* `conda <https://docs.conda.io/en/latest/miniconda.html>`__
* python 3.7

Installation from Source
------------------------

1. Get the code.

.. code:: bash

	git clone https://github.com/amaslan/dimelo
	cd dimelo

2. Create and activate a python 3.7 conda venv.

Recommended to avoid package dependency install issues:

Linux:

.. code:: bash

	conda env create -f environment_linux.yml
	conda activate dimelo

MacOS:

.. code:: bash

	conda env create -f environment_macOS.yml
	conda activate dimelo

Other option to create venv without using provided yml:

.. code:: bash

	conda create --name dimelo python=3.7
	conda activate dimelo

3. Install dimelo and its requirements.

.. code:: bash

	pip install .

4. Now the dimelo package is available for import. For example, you may import it like this: 

.. code:: bash

	python

>>> import dimelo as dm

And then call functions like this:

>>> dm.function_name(...)


Other notes
------------------------

* The dependency pybedtools requires that you have bedtools. Because this is not able to be installed with pip, if you did not create your environment from environment.yml, you must run:

.. code:: bash

	conda install bedtools

* If you would like to make changes to the package, you can install it in developer mode:

.. code:: bash

	pip install -e .


