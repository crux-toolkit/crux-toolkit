#!/bin/bash

rm -f good_results/*.observed
cucumber "$@"
echo

for f in $(ls good_results/*.observed 2>/dev/null); do
  echo "$f --> ${f%.*}"
  mv "$f" "${f%.*}"
done
echo

