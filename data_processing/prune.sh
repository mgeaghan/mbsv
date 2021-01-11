#!/bin/bash

# Step 1 - prune all variants
plink --bfile 1kg.auto.eur.maf_0.01.nolrld --keep-allele-order --extract gwas.nonmhc.1kgvar.txt --make-bed --out 1kg.auto.eur.maf_0.01.nolrld.gwas
plink --bfile 1kg.auto.eur.maf_0.01.nolrld.gwas --keep-allele-order --indep-pairwise 50 5 0.1 --out 1kg.auto.eur.maf_0.01.nolrld.gwas.r2_0.1.50_5.ld-pruned

# Step 2 - prune mbsvs
plink --bfile 1kg.auto.eur.maf_0.01.nolrld --keep-allele-order --extract mbsv.0.2.nonmhc.1kgvar.txt --make-bed --out 1kg.auto.eur.maf_0.01.nolrld.mbsv
plink --bfile 1kg.auto.eur.maf_0.01.nolrld.mbsv --keep-allele-order --indep-pairwise 50 5 0.1 --out 1kg.auto.eur.maf_0.01.nolrld.mbsv.r2_0.1.50_5.ld-pruned

# Step 3 - get list of pruned variants in LD with pruned mbsvs
plink --bfile 1kg.auto.eur.maf_0.01.nolrld.gwas --keep-allele-order --ld-snp-list 1kg.auto.eur.maf_0.01.nolrld.mbsv.r2_0.1.50_5.ld-pruned.prune.in --ld-window 50 --ld-window-r2 0.1 --out 1kg.auto.eur.maf_0.01.nolrld.gwas.r2_0.1.50_5.ld-pruned.ld.mbsv.r2_0.1_50 --r2
sed 's/^\ \+//g' 1kg.auto.eur.maf_0.01.nolrld.gwas.r2_0.1.50_5.ld-pruned.ld.mbsv.r2_0.1_50.ld | sed 's/\ \+/\t/g' | cut -f 6 | sort | uniq | awk 'NR>1' > 1kg.auto.eur.maf_0.01.nolrld.gwas.r2_0.1.50_5.ld-pruned.ld.mbsv.r2_0.1_50.ld_remove

# Step 4 - remove these variants from the list of all pruned variants and replace with the pruned mbsvs
grep -vwFf 1kg.auto.eur.maf_0.01.nolrld.gwas.r2_0.1.50_5.ld-pruned.ld.mbsv.r2_0.1_50.ld_remove 1kg.auto.eur.maf_0.01.nolrld.gwas.r2_0.1.50_5.ld-pruned.prune.in > gwas.nonmhc.pruned.1kgvar.txt
cp 1kg.auto.eur.maf_0.01.nolrld.mbsv.r2_0.1.50_5.ld-pruned.prune.in mbsv.0.2.nonmhc.pruned.1kgvar.txt
cat gwas.nonmhc.pruned.1kgvar.txt mbsv.0.2.nonmhc.pruned.1kgvar.txt > gwas.mbsv.0.2.nonmhc.pruned.1kgvar.txt

