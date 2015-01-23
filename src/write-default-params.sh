#!/bin/bash

# A script to create the default.params file to distribute in the docs
# and to store in src.  Writes a param file with default values.
# Adds a header.  Copies to documents directory.
# It is up to the developer to run this script after any changes are made to 
# the parameters.

SCRIPT_DIR="$(dirname "$BASH_SOURCE")"
ROOT_DIR="$SCRIPT_DIR/.."

OUTPUT_DIR="crux-output"
TARGET_PARAM_FILE="$ROOT_DIR/src/default.params"
TARGET_DOC_PARAM_FILE="$ROOT_DIR/doc/user/default.params"

if [ -d "$OUTPUT_DIR" ]; then
  echo "Directory $OUTPUT_DIR already exists, remove it and re-run script"
  exit
fi

# Create param file
$ROOT_DIR/src/crux comet "$ROOT_DIR/test/smoke-tests/test.ms2" "$ROOT_DIR/test/smoke-tests/test.fasta" 1>/dev/null 2>&1

if [ ! -e "$OUTPUT_DIR/comet.params.txt" ]; then
  echo "Failed to write the param file."
  exit
fi

# Write a header
echo "####################################################################
# Sample parameter file
#
# Lines starting with '#' will be ignored.  Don't leave any space
# before or after a parameter setting.  The format is
#
# <parameter-name>=<value>
#
#####################################################################
" > "$TARGET_PARAM_FILE"
cat "$OUTPUT_DIR/comet.params.txt" >> "$TARGET_PARAM_FILE"

# Rename the value for parameter-file, and copy to doc directory.
# N.B. Don't use the "-i" option to 'sed' because it has different
# syntax under Linux and Darwin.  Instead, combine the sed and the copy,
# then copy the sed output back to the local directory.
sed 's/parameter-file=T/parameter-file=F/' "$TARGET_PARAM_FILE" > "$TARGET_DOC_PARAM_FILE"
cp "$TARGET_DOC_PARAM_FILE" "$TARGET_PARAM_FILE"

