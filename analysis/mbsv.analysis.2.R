# Merge dbMTS results with GWAS data
# Author: Michael Geaghan (2020)

library(reshape2)
library(ggplot2)
load("dbmts.RData")

# Read in GWAS data
gwas_files <- list(adhd.2018 = '../gwas/adhd-2018.gwas',
                   an.2019 = '../gwas/an-2019.gwas',
                   asd.2019 = '../gwas/asd-2019.gwas',
                   bip.2018 = '../gwas/bip-2018.gwas',
                   mdd.2019 = '../gwas/mdd-2019.gwas',
                   ocd.2018 = '../gwas/ocd-2018.gwas',
                   ptsd.2019 = '../gwas/ptsd-2019.gwas',
                   scz.clozuk2018 = '../gwas/scz-clozuk2018.gwas',
                   ts.2019 = '../gwas/ts-2019.gwas')

gwas_es <- list(adhd.2018 = '',
                an.2019 = '',
                asd.2019 = '',
                bip.2018 = '',
                mdd.2019 = '',
                ocd.2018 = '',
                ptsd.2019 = '',
                scz.clozuk2 = '',
                ts.2019 = '')

# Add MHC information
mhc_file <- '../db/g1000_eur.mhc.snps'
mhc <- scan(mhc_file, character())

for(g in names(gwas_files)) {
  gwas <- read.delim(gwas_files[[g]], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Filter out unwanted GWAS columns and standardise names
  gwas_filter_cols <- c("SNP", "OR", "BETA", "LOG_ODDS", "LogOR", "P", "PVAL")
  gwas <- gwas[colnames(gwas) %in% gwas_filter_cols]
  if(colnames(gwas)[2] == "OR") {
    gwas$OR <- log2(gwas$OR)
    gwas_es[[g]] <- "log OR"
  } else if(colnames(gwas)[2] == "LOG_ODDS") {
    gwas_es[[g]] <- "log OR"
  } else if (colnames(gwas)[2] == "LogOR") {
    gwas_es[[g]] <- "log OR"
  } else if(colnames(gwas)[2] == "BETA") {
    gwas_es[[g]] <- "Beta"
  }
  colnames(gwas) <- c("SNP", "ES", "P")
  # Add dbMTS information to GWAS data frames
  gwas <- merge(gwas, dbmts.agree.best.best.short.var.score, by.x = 'SNP', by.y = 'dbsnp', all.x = TRUE)
  colnames(gwas) <- c("SNP", "ES", "P", "variant", "mbsv.gtex.eqtl", "mbsv.gtex.brain.eqtl", "mbsv.psychencode.eqtl", "mbsv.any.eqtl", "diffScore", "non.mhc")
  gwas$mbsv.gtex.eqtl[is.na(gwas$mbsv.gtex.eqtl)] <- FALSE
  gwas$mbsv.gtex.brain.eqtl[is.na(gwas$mbsv.gtex.brain.eqtl)] <- FALSE
  gwas$mbsv.psychencode.eqtl[is.na(gwas$mbsv.psychencode.eqtl)] <- FALSE
  gwas$mbsv.any.eqtl[is.na(gwas$mbsv.any.eqtl)] <- FALSE
  gwas$non.mhc <- !(gwas$SNP %in% mhc)
  gwas$GWAS <- TRUE
  gwas$mbsv <- (!is.na(gwas$diffScore))
  # Generate reshaped gWAS data frame
  gwas.df <- melt(gwas, id.vars = c("SNP", "ES", "P", "variant", "diffScore", "non.mhc"))
  gwas.df <- gwas.df[gwas.df$value,]
  # Save data to disk
  gwas_files_char = gwas_files[[g]]
  gwas_es_char = gwas_es[[g]]
  save(gwas, gwas.df, gwas_files_char, gwas_es_char, file = paste(g, ".gwas.RData", sep = ""))
}
