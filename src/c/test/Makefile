# This makefile controls the smoke-tests and unit-tests

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
