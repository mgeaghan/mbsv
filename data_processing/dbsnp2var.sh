#!/bin/bash
awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$1]=$2; b[$1]=$3; c[$1]=$4; d[$1]=$5; next} NR>FNR {if ($0 in a) print a[$0]; print b[$0]; print c[$0]; print d[$0]}' ../db/g1000_eur.magma.convert $1 | sort | uniq | grep -vP "^$" > $2

