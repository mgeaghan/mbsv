#!/bin/bash
DB=$1
POS=$2
ANNDIR=$3
OUT=$4
if [[ "$DB" == "" || "$POS" == "" || "$ANNDIR" == "" || "$OUT" == "" || ! -f "$DB" || ! -f "${DB}.tbi" || ! -f "$POS" || ! -d "$ANNDIR" ]]; then
	echo "Usage: ./extract_annot.sh <CADD_DATABASE_FILE> <POSITION_FILE> <ANNOVAR_ANNOTATION_DIR> <OUTPUT_PREFIX>"
	echo "Position file must contain three tab-delimited columns: CHR, START, END. END should match START to select a single variant; if END is greater than START, all variants in the range will be selected."
else
	tabix $DB -R $POS > ${OUT}.cadd
	awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, ($2 + (length($3) - 1)), $3, $4, $5, $6}' ${OUT}.cadd > ${OUT}.avinput
	table_annovar.pl ${OUT}.avinput ${ANNDIR} -buildver hg19 -out ${OUT}.annovar -remove -otherinfo -protocol gnomad211_exome,gnomad211_genome -operation f,f -nastring .
fi

