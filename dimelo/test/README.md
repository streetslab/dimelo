# test

This folder contains data and code to run basic tests on the `dimelo` package. 

`__init__.py` sets up a framework for running parsing, including downloading a reference genome and processing input files appropriately.

`dimelo_test.py` implements unit tests and integration tests using `pytest`

`generate_test_targets.ipynb` contains code to create the target outputs for the unit tests. These should not be updated unless you confirm that any change in behavior is actually correct. *Special care should be taken for the `load_processed` outputs, which should not change with interface changes. If those outputs don't match and need to be regenerated, that is* ***reason for concern.*** However, updates to e.g. the .h5 single read format and corresponding changes to the `load_processed` methods may require making a new target .h5 file while `load_processed.binarized_read_from_hdf5` still returns the right array values and so on.

The `data` folder contains .bam files and .bed files to use for testing. These files are also used by the tutorial.

The `output` folder stores the reference genome and processed outputs. This folder is included in `.gitignore` so its contents should never be included in commits or merges.

```
# Ignore tutorial output files
test/output
```