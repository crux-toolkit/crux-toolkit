# This makefile controls the smoke, unit and performance tests.

# NOTE: When debugging unit tests using the check library you usually
# want to set the environment variable CK_FORK to no. This keeps the
# check library from running the tests on multiple threads which makes
# them hard to debug.

# Run all tests
check: check-unit-tests check-cpp-unit-tests check-smoke-tests performance-tests

# Build the unit test executable
unit-tests:
	cd unit-tests && $(MAKE)

# Run just the unit tests
check-unit-tests:
	cd unit-tests && $(MAKE) check

# Build the cpp-unit test executable
cpp-unit-tests:
	cd cppunit-tests && ./configure && $(MAKE)

# Run just the cpp-unit tests
check-cpp-unit-tests:
	cd cppunit-tests && $(MAKE) check

# Run just the smoke tests
check-smoke-tests:
	cd smoke-tests && ./runall

# Run just the performance tests
performance-tests:
	cd performance-tests && ./run-performance-test.sh

# Remove generated files in unit-test
clean:
	cd unit-tests && $(MAKE) clean

.PHONY: check check-smoke-tests check-unit-tests clean unit-tests
