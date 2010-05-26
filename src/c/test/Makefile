# This makefile controls the smoke-tests and unit-tests

# NOTE: When debugging unit tests using the check library you usually
# want to set the environment variable CK_FORK to no. This keeps the
# check library from running the tests on multiple threads which makes
# them hard to debug.

# Build the unit test executable
unit-tests:
	cd unit-tests && $(MAKE)

# Run the unit tests and the smoke tests
check: check-unit-tests check-smoke-tests

# Run just the unit tests
check-unit-tests:
	cd unit-tests && $(MAKE) check

# Run just the smoke tests
check-smoke-tests:
	cd smoke-tests && ./runall

# Remove generated files in unit-test
clean:
	cd unit-tests && $(MAKE) clean

.PHONY: check check-smoke-tests check-unit-tests clean unit-tests
