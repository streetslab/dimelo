Installation
====================

Requirements
------------

* `conda <https://docs.conda.io/en/latest/miniconda.html>`__
* python 3.7

Installation from Pip
---------------------

1. Create and activate a pytion 3.7 conda venv:

.. code:: bash

	conda create --name dimelo python=3.7 pip
	source activate dimelo

2. Install dimelo from Pypi:

.. code:: bash

	pip install -i https://test.pypi.org/simple/dimelo

Installation from Source
------------------------

1. Create and activate a pytion 3.7 conda venv:

.. code:: bash

	conda create --name dimelo python=3.7 pip
	source activate dimelo


2. Get the code:

.. code:: bash

	git clone https://github.com/amaslan/dimelo
	cd dimelo


3. Install dimelo and its requirements

.. code:: bash

	pip install -e .
