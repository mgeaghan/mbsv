#!/bin/bash
awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]=$2; next} NR>FNR {if ($0 in a) print a[$0]}' convert.vars.txt $1 > $2

