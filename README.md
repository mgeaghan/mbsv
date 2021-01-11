# mbsv
Scripts for microRNA binding site variant (MBSV) data analysis

## dbmts2long.py

This is a python script used to convert the output from dbMTS to a long-table format ready for analysis with the downstream R scripts.

## analysis/

The R scripts in this directory form the main data analysis steps from reading and filtering the raw dbMTS database output to running the statistical analyses and processing the MAGMA results.

## data_processing

These bash scripts are responsible for converting variant IDs between dbSNP IDs and CHR:POS:REF:ALT IDs as well as variant pruning and annotation.
