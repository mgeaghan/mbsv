#!/bin/bash
awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {a[$2]=$1; b[$3]=$1; c[$4]=$1; d[$5]=$1; next} NR>FNR && FNR==1 {print $0, "variant", "dbsnp"; next} NR>FNR {v=($1":"$2":"$4":"$5); if (v in a) print $0, v, a[v]; else if (v in b) print $0, v, b[v]; else if (v in c) print $0, v, c[v]; else if (v in d) print $0, v, d[v]; else print $0, v, "."}' g1000_eur.magma.convert $1 > $2

