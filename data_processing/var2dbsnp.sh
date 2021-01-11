#!/bin/bash
awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$2]=$1; b[$3]=$1; c[$4]=$1; d[$5]=$1; next} NR>FNR {if ($0 in a) print a[$0]; else if ($0 in b) print b[$0]; else if ($0 in c) print c[$0]; else if ($0 in d) print d[$0]}' ../db/g1000_eur.magma.convert $1 > $2

