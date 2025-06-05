# test

This folder contains data and code to run basic tests on the `dimelo` package. `dimelo_test.py` currently contains all test code, and `data/test_targets/test_matrix.pickle` contains a mapping of test case kwargs to target values, including paths to files living in the same directory.

New versions of `test_matrix.pickle` and the accompanying target files to which it refers can be generated using `generate_targets.py`, which itself can
reference both existing results in `test_matrix.pickle` as well as a fresh set of test cases (i.e. kwarg combinations) defined in `cases.py`.

## files

`__init__.py` sets up a framework for running parsing, including downloading a reference genome and processing input files appropriately.

`cases.py` contains a set of test cases which will be run through all tests in `dimelo_test.py`, and used to generate corresponding targets in `generate_targets.py`. The format schema contains kwargs by name in a dictionary for each test case: these are passed *directly* to the parse, load, plot, and export functions, simply filtering out un-needed ones. Formats must match the requirements of the functions that take these arguments.

`dimelo_test.py` implements unit tests and integration tests using `pytest`. Tests are split into classes to handle temporary directories cleanly and separate different types of tests. 

`generate_targets.py` contains code to create the target outputs for the unit tests, ultimately creating a test_matrix.pkl file containing a pickled directionary of test kwargs and results. These should not be updated unless you confirm that any change in behavior is actually correct, i.e. if any test is failing make sure you know why before considering replacing the target values. *Special care should be taken for the `load_processed` outputs, which should not change with interface changes. If those outputs don't match and need to be regenerated, that is* ***reason for concern.*** However, updates to e.g. the .h5 single read format and corresponding changes to the `load_processed` methods may require making a new target .h5 file while `load_processed.binarized_read_from_hdf5` still returns the right array values and so on. Run `python generate_targets.py --help` for assistance with arguments to update only a subset of target value or test cases.

The `data` folder contains .bam files and .bed files to use for testing. These files are also used by the tutorial.

The `output` folder stores the reference genome and processed outputs. This folder is included in `.gitignore` so its contents should never be included in commits or merges.

```
# Ignore tutorial output files
test/output
```