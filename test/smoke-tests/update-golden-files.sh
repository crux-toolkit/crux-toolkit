#!/usr/bin/env bash
# update-golden-files.sh
# Run from the smoke-tests BUILD directory (e.g. Release/test/smoke-tests/).
# Sets up the build directory for running cucumber smoke tests and updates
# golden files in both the source tree and the build directory.

set -euo pipefail

CRUX=../../src/crux
SMOKE_SRC="$(cd "$(dirname "$0")" && pwd)"
GOOD_RESULTS_SRC="$SMOKE_SRC/good_results"
BUILD_DIR="$(pwd)"
GOOD_RESULTS_BUILD="$BUILD_DIR/good_results"

echo "=== Setting up build directory ==="

# Symlink features/ so fileComparator.py and step definitions are accessible
if [ ! -e "$BUILD_DIR/features" ]; then
  ln -sf "$SMOKE_SRC/features" "$BUILD_DIR/features"
  echo "Created symlink: features -> $SMOKE_SRC/features"
fi

# Copy input files not handled by cmake
cp "$SMOKE_SRC"/sample2.* "$BUILD_DIR/"
cp -r "$SMOKE_SRC/test_results" "$BUILD_DIR/"

echo "=== Updating LFQ golden files ==="

$CRUX lfq --no-analytics T --overwrite T --output-dir lfq-no-norm \
  --psm-file-format assign-confidence --is-psm-filtered T \
  lfq-psms-no-norm.tsv sliced-mzml.mzML

$CRUX lfq --no-analytics T --overwrite T --output-dir lfq-with-norm \
  --normalize T --psm-file-format assign-confidence --is-psm-filtered T \
  --specfile-replicates lfq-specfile-replicates.tsv \
  lfq-psms-with-norm.tsv \
  20100614_Velos1_TaGe_SA_K562_3.mzML \
  20100614_Velos1_TaGe_SA_K562_4.mzML

cp lfq-no-norm/crux-lfq-peaks.txt    "$GOOD_RESULTS_SRC/lfq-no-norm-peaks.txt"
cp lfq-no-norm/crux-lfq-mod-pep.txt  "$GOOD_RESULTS_SRC/lfq-no-norm-mod-pep.txt"
cp lfq-with-norm/crux-lfq-peaks.txt  "$GOOD_RESULTS_SRC/lfq-with-norm-peaks.txt"
cp lfq-with-norm/crux-lfq-mod-pep.txt "$GOOD_RESULTS_SRC/lfq-with-norm-mod-pep.txt"

echo "=== Updating make-pin golden files ==="

$CRUX make-pin --no-analytics T --overwrite T \
  --output-file make-pin_txt.pin sample2.search.target.txt
$CRUX make-pin --no-analytics T --overwrite T \
  --output-file make-pin_pep.pin sample2.search.target.pep.xml

cp crux-output/make-pin_txt.pin "$GOOD_RESULTS_SRC/make-pin_txt.pin"
cp crux-output/make-pin_pep.pin "$GOOD_RESULTS_SRC/make-pin_pep.pin"

echo "=== Updating psm-convert golden files ==="

$CRUX psm-convert --no-analytics T --overwrite T \
  test_results/results1.tide-search.txt pin
cp crux-output/psm-convert.pin "$GOOD_RESULTS_SRC/psmconv-from-txt1.pin"

$CRUX psm-convert --no-analytics T --overwrite T \
  test_results/results2.tide-search.txt pin
cp crux-output/psm-convert.pin "$GOOD_RESULTS_SRC/psmconv-from-txt2.pin"

$CRUX psm-convert --no-analytics T --overwrite T \
  test_results/results1.tide-search.txt pepxml
cp crux-output/psm-convert.pep.xml "$GOOD_RESULTS_SRC/psmconv-from-txt1.pep.xml"

$CRUX psm-convert --no-analytics T --overwrite T \
  test_results/results2.tide-search.txt pepxml
cp crux-output/psm-convert.pep.xml "$GOOD_RESULTS_SRC/psmconv-from-txt2.pep.xml"

echo "=== Syncing build good_results directory ==="

cp "$GOOD_RESULTS_SRC/lfq-no-norm-peaks.txt"     "$GOOD_RESULTS_BUILD/"
cp "$GOOD_RESULTS_SRC/lfq-no-norm-mod-pep.txt"   "$GOOD_RESULTS_BUILD/"
cp "$GOOD_RESULTS_SRC/lfq-with-norm-peaks.txt"   "$GOOD_RESULTS_BUILD/"
cp "$GOOD_RESULTS_SRC/lfq-with-norm-mod-pep.txt" "$GOOD_RESULTS_BUILD/"
cp "$GOOD_RESULTS_SRC/make-pin_txt.pin"           "$GOOD_RESULTS_BUILD/"
cp "$GOOD_RESULTS_SRC/make-pin_pep.pin"           "$GOOD_RESULTS_BUILD/"
cp "$GOOD_RESULTS_SRC/psmconv-from-txt1.pin"      "$GOOD_RESULTS_BUILD/"
cp "$GOOD_RESULTS_SRC/psmconv-from-txt2.pin"      "$GOOD_RESULTS_BUILD/"
cp "$GOOD_RESULTS_SRC/psmconv-from-txt1.pep.xml"  "$GOOD_RESULTS_BUILD/"
cp "$GOOD_RESULTS_SRC/psmconv-from-txt2.pep.xml"  "$GOOD_RESULTS_BUILD/"

echo "=== Done. Run cucumber to verify: ==="
echo "cucumber \\"
echo "  --require $SMOKE_SRC/features/support \\"
echo "  --require $SMOKE_SRC/features/step_definitions \\"
echo "  $SMOKE_SRC/features/lfq.feature \\"
echo "  $SMOKE_SRC/features/make-pin.feature \\"
echo "  $SMOKE_SRC/features/psm-convert.feature"
