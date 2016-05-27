#!/bin/bash
function usage {
    echo "Usage:"
    echo " $0 [ -l ] <file_to_parse>"
    echo "  [ -l ] Run the lint with lenient checking"
    exit
}
lflag='+'
while getopts 'l' flag; do
    case "${flag}" in
	l) lflag='-' ;;
	*) usage ;;
    esac
done
shift "$((OPTIND-1))"
[[ -n "$1" ]] || usage

include="include_rules.txt"
exclude="exclude_rules.txt"
lax="include_lax_rules.txt"
root="../../src/app"
filter=""
file_to_parse=$1

while IFS='' read -r line || [[ -n "$line" ]]; do
    filter+="+$line,"
done < "$lax"
while IFS='' read -r line || [[ -n "$line" ]]; do
    filter+="$lflag$line,"
done < "$include"
while IFS='' read -r line || [[ -n "$line" ]]; do
    filter+="-$line,"
done < "$exclude"
$(python cpplint/cpplint.py "--root=$root" "--filter=$filter" $file_to_parse)
