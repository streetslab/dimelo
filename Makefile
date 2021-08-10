define help

Supported targets: prepare, develop, sdist, clean, test, and pypi.

The 'prepare' target installs this project's build requirements into the current virtualenv.

The 'develop' target creates an editable install of this project and its runtime requirements in the
current virtualenv. The install is called 'editable' because changes to the source code
immediately affect the virtualenv.

The 'clean' target undoes the effect of 'develop'.

The 'test' target runs unit tests. Set the 'tests' variable to run a particular test, e.g.

        make test tests=PlotMAPQtest/countMAPQ_test.py

The 'pypi' target publishes the current commit of this project to PyPI after enforcing that the working
copy and the index are clean, and tagging it as an unstable .dev build.

endef
export help
help:
	@printf "$$help"

SHELL=bash
python=python
pip=pip
tests=.
version:=$(shell $(python) version.py)
sdist_name:=dimelo-$(version).tar.gz

develop:
	$(pip) install -e .

clean_develop:
	- $(pip) uninstall -y dimelo
	- rm -rf *.egg-info

clean_sdist:
	- rm -rf dist

clean: clean_develop clean_pypi

check_build_reqs:
	@$(python) -c 'import pytest' \
                || ( printf "$(redpip)Build requirements are missing. Run 'make prepare' to install them.$(normal)" ; false )

test: check_build_reqs
	$(python) -m pytest -vv $(tests)

pypi: clean clean_sdist
	set -x \
	&& $(python) setup.py sdist bdist_wheel \
	&& twine check dist/* \
	&& twine upload --repository-url https://test.pypi.org/legacy/ dist/*
clean_pypi:
	- rm -rf build/