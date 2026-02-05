# Notes on running unit tests.

- Due to the large size of spectrum files, these files are not added to version control rather they should be downloaded by the developer by following this: [ProteomeXchange - PXD005590](https://www.ebi.ac.uk/pride/archive/projects/PXD005590).
- In order to reduce build time, PXD005590 folder will not be added to the test-data folder. Rather set the varibale in the code (download location of test data): /tests/unit-tests/crux_lfq_test.cpp
  > string spectrum_files = "";
