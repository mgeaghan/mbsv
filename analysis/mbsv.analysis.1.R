# Process and score MBSVs from dbMTS
# Author: Michael Geaghan (2020)

dbmts_file <- '../dbmts/psych.gwas.out.dbMTS.long'
dbmts <- read.delim(dbmts_file, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
dbmts_orig <- dbmts
cols <- c('chr_38', 'pos_38', 'ref_38', 'alt_38', 'chr', 'pos', 'ref', 'alt',
          'vep_cons', 'vep_trans_id', 'vep_gene_name', 'vep_gene_id', 'vep_cann',
          'dbsnp_150', 'gtex_gene_id.version', 'gtex_tissue',
          'algorithm', 'allele', 'score', 'mirna', 'transcript',
          'corr_tum', 'tum_tiss', 'corr_norm', 'norm_tiss')
colnames(dbmts) <- cols

# Filter out X chromosome
keep <- dbmts$chr != 'X'
dbmts <- dbmts[keep, ]

# Filter out rows without a miRNA/score
keep <- dbmts$mirna != '.'
dbmts <- dbmts[keep, ]

# Filter out non-brain-expressed miRNAs
mirnaList <- scan('../db/filteredmirs.conf.txt', character())
keep <- dbmts$mirna %in% mirnaList
dbmts <- dbmts[keep, ]

# Filter out non-brain-expressed genes
geneList <- scan('../db/filteredgenes.txt', character())
# Download gene and transcript IDs and transcript types from biomart (GRCh37)
library(biomaRt)
ensembl <- useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl', host = 'http://grch37.ensembl.org')
bm <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_gene_id_version', 'ensembl_transcript_id_version', 'transcript_biotype'), mart = ensembl)
bm.all <- bm
bm <- bm[bm$transcript_biotype=='protein_coding',]
# Keep only transcripts whose genes are expressed in the brain
transcriptList <- bm$ensembl_transcript_id[bm$ensembl_gene_id %in% geneList]
keep <- dbmts$transcript %in% transcriptList
dbmts <- dbmts[keep, ]

# Add a CHR:POS:REF:ALT variant name column
dbmts$variant <- paste(dbmts$chr, dbmts$pos, dbmts$ref, dbmts$alt, sep = ":")

# Add a gene ID column
dbmts$gene <- bm$ensembl_gene_id[match(dbmts$transcript, bm$ensembl_transcript_id)]

# Add eQTL information
dbmts.gtex <- unique(data.frame(variant = dbmts$variant, gtex_gene = dbmts$gtex_gene_id.version, gtex_tiss = dbmts$gtex_tissue))
gtex.df <- do.call('rbind', apply(dbmts.gtex, 1, function(x) {
  var <- as.character(x[1])
  gene.str <- as.character(x[2])
  tiss.str <- as.character(x[3])
  if(gene.str == ".") {
    gene.list <- NA
  } else {
    gene.list <- gsub("\\.\\d+$", "", unlist(strsplit(gene.str, split = "|", fixed = TRUE)))
  }
  if(tiss.str == ".") {
    tiss.list <- NA
  } else {
    tiss.list <- unlist(strsplit(tiss.str, split = "|", fixed = TRUE))
  }
  if(is.na(gene.list) && is.na(tiss.list)) {
    return(data.frame(variant = character(), gene = character(), tissue = character()))
  } else {
    return(data.frame(variant = var, gene = gene.list, tissue = tiss.list))
  }
}))
gtex.any.eqtl <- paste(gtex.df$variant, gtex.df$gene, sep = "_")
gtex.brain.eqtl <- gtex.any.eqtl[grepl("[Bb]rain", gtex.df$tissue, perl = TRUE)]
gtex.any.eqtl <- unique(gtex.any.eqtl)
gtex.brain.eqtl <- unique(gtex.brain.eqtl)
# Read in PSYCHENCODE data
psychencode_file <- '../db/psychencode.gene.snp.txt'
psychencode <- read.delim(psychencode_file, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
colnames(psychencode) <- c('gene', 'variant')
psychencode.eqtl <- unique(paste(psychencode$variant, psychencode$gene, sep = "_"))
# Add a TRUE/FALSE column for each of: 'gtex.eqtl', 'gtex.brain.eqtl', and 'psychencode.eqtl'
dbmts.var.gen <- paste(dbmts$variant, dbmts$gene, sep = "_")
dbmts$gtex.eqtl <- dbmts.var.gen %in% gtex.any.eqtl
dbmts$gtex.brain.eqtl <- dbmts.var.gen %in% gtex.brain.eqtl
dbmts$psychencode.eqtl <- dbmts.var.gen %in% psychencode.eqtl
dbmts$any.eqtl <- dbmts$gtex.eqtl | dbmts$gtex.brain.eqtl | dbmts$psychencode.eqtl

# Keep targetscan entries that agree with at least one other algorithm
dbmts.agree <- dbmts
dbmts.agree$var.all.mir.gen <-  paste(dbmts.agree$variant, dbmts.agree$allele, dbmts.agree$mirna, dbmts.agree$gene, sep = "_")
dbmts.agree.t <- dbmts.agree$var.all.mir.gen[dbmts.agree$algorithm == 't']
dbmts.agree.m <- dbmts.agree$var.all.mir.gen[dbmts.agree$algorithm == 'm']
dbmts.agree.r <- dbmts.agree$var.all.mir.gen[dbmts.agree$algorithm == 'r']
keep <- (dbmts.agree.t %in% dbmts.agree.m) | (dbmts.agree.t %in% dbmts.agree.r)
dbmts.agree.t <- dbmts.agree.t[keep]
dbmts.agree <- dbmts.agree[dbmts.agree$var.all.mir.gen %in% dbmts.agree.t & dbmts.agree$algorithm == 't',]

# Calculate the diffScores for every variant-mirna-gene combination
dbmts.agree.diff <- dbmts.agree
dbmts.agree.diff$score <- as.numeric(dbmts.agree.diff$score)
dbmts.agree.diff$var.mir.gen <- paste(dbmts.agree.diff$variant, dbmts.agree.diff$mirna, dbmts.agree.diff$gene, sep = "_")
dbmts.agree.diff.var.mir.gen <- unique(dbmts.agree.diff$var.mir.gen)
dbmts.agree.diff.ref <- dbmts.agree.diff[dbmts.agree.diff$allele == 'ref',]
dbmts.agree.diff.alt <- dbmts.agree.diff[dbmts.agree.diff$allele == 'alt',]
dbmts.agree.diff.var.mir.gen.diffScore <- do.call('c', lapply(dbmts.agree.diff.var.mir.gen, function(x) {
  alt <- unique(dbmts.agree.diff.alt$score[dbmts.agree.diff.alt$var.mir.gen == x])
  ref <- unique(dbmts.agree.diff.ref$score[dbmts.agree.diff.ref$var.mir.gen == x])
  if (length(alt) == 0) {
    alt = 0
  } else if (length(alt) > 1) {
    alt = min(alt)  # Take the best scores across all transcripts
  }
  if (length(ref) == 0) {
    ref = 0
  } else if (length(ref) > 1) {
    ref = min(ref)  # Take the best scores across all transcripts
  }
  return(alt - ref)
}))
dbmts.agree.diff.var.mir.gen.diffScore.df <- data.frame(var.mir.gen = dbmts.agree.diff.var.mir.gen, diffScore = dbmts.agree.diff.var.mir.gen.diffScore)
dbmts.agree.diff$diffScore <- dbmts.agree.diff.var.mir.gen.diffScore.df$diffScore[match(dbmts.agree.diff$var.mir.gen, dbmts.agree.diff.var.mir.gen.diffScore.df$var.mir.gen)]

# Calculate the best score for every variant-allele-gene combination
dbmts.agree.best <- dbmts.agree
dbmts.agree.best$score <- as.numeric(dbmts.agree.best$score)
dbmts.agree.best$var.all.gen <- paste(dbmts.agree.best$variant, dbmts.agree.best$allele, dbmts.agree.best$gene, sep = "_")
dbmts.agree.best.var.all.gen <- unique(dbmts.agree.best$var.all.gen)
dbmts.agree.best.var.all.gen.score <- do.call('c', lapply(dbmts.agree.best.var.all.gen, function(x) min(dbmts.agree.best$score[dbmts.agree.best$var.all.gen == x])))
dbmts.agree.best.var.all.gen.score.df <- data.frame(var.all.gen = dbmts.agree.best.var.all.gen, score = dbmts.agree.best.var.all.gen.score)
dbmts.agree.best$best.score <- dbmts.agree.best.var.all.gen.score.df$score[match(dbmts.agree.best$var.all.gen, dbmts.agree.best.var.all.gen.score.df$var.all.gen)]
keep <- dbmts.agree.best$score == dbmts.agree.best$best.score
dbmts.agree.best <- dbmts.agree.best[keep,]

# Calculate the best-vs-best score for every variant-gene combination
dbmts.agree.best$var.gen <- paste(dbmts.agree.best$variant, dbmts.agree.best$gene, sep = "_")
dbmts.agree.best.var.gen <- unique(dbmts.agree.best$var.gen)
dbmts.agree.best.ref <- dbmts.agree.best[dbmts.agree.best$allele == 'ref',]
dbmts.agree.best.alt <- dbmts.agree.best[dbmts.agree.best$allele == 'alt',]
dbmts.agree.best.var.gen.diffScore <- do.call('c', lapply(dbmts.agree.best.var.gen, function(x) {
  alt <- unique(dbmts.agree.best.alt$score[dbmts.agree.best.alt$var.gen == x])
  ref <- unique(dbmts.agree.best.ref$score[dbmts.agree.best.ref$var.gen == x])
  if (length(alt) == 0) {
    alt = 0
  }
  if (length(ref) == 0) {
    ref = 0
  }
  return(alt - ref)
}))
dbmts.agree.best.var.gen.diffScore.df <- data.frame(var.gen = dbmts.agree.best.var.gen, diffScore = dbmts.agree.best.var.gen.diffScore)
dbmts.agree.best$diffScore <- dbmts.agree.best.var.gen.diffScore.df$diffScore[match(dbmts.agree.best$var.gen, dbmts.agree.best.var.gen.diffScore.df$var.gen)]

# Keep the best best-vs-best score for every variant (in the event of one variant affecting multiple genes)
dbmts.agree.best.best <- dbmts.agree.best
dbmts.agree.best.best.var <- unique(dbmts.agree.best.best$variant)
dbmts.agree.best.best.var.bestDiffScore <- do.call('c', lapply(dbmts.agree.best.best.var, function(x) max(abs(dbmts.agree.best.best$diffScore[dbmts.agree.best.best$variant == x]))))
dbmts.agree.best.best.var.bestDiffScore.df <- data.frame(var = dbmts.agree.best.best.var, bestDiffScore = dbmts.agree.best.best.var.bestDiffScore)
dbmts.agree.best.best$best.diff.score <- dbmts.agree.best.best.var.bestDiffScore.df$bestDiffScore[match(dbmts.agree.best.best$variant, dbmts.agree.best.best.var.bestDiffScore.df$var)]
keep <- abs(dbmts.agree.best.best$diffScore) == dbmts.agree.best.best$best.diff.score
dbmts.agree.best.best <- dbmts.agree.best.best[keep,]

# Filter gene and transcript diffScores, best-vs-best diffScores and best best-vs-best diffScores by an absolute value of 0.2
keep <- abs(dbmts.agree.diff$diffScore) >= 0.2
dbmts.agree.diff.0.2 <- dbmts.agree.diff[keep,]
keep <- abs(dbmts.agree.best$diffScore) >= 0.2
dbmts.agree.best.0.2 <- dbmts.agree.best[keep,]
keep <- abs(dbmts.agree.best.best$diffScore) >= 0.2
dbmts.agree.best.best.0.2 <- dbmts.agree.best.best[keep,]

# Create short versions of each data frame
short.colnames <- c("variant", "allele", "mirna", "gene", "diffScore", "gtex.eqtl", "gtex.brain.eqtl", "psychencode.eqtl", "any.eqtl")
dbmts.agree.diff.short <- unique(dbmts.agree.diff[colnames(dbmts.agree.diff) %in% short.colnames])
dbmts.agree.best.short <- unique(dbmts.agree.best[colnames(dbmts.agree.best) %in% short.colnames])
dbmts.agree.best.best.short <- unique(dbmts.agree.best.best[colnames(dbmts.agree.best.best) %in% short.colnames])
dbmts.agree.diff.0.2.short <- unique(dbmts.agree.diff.0.2[colnames(dbmts.agree.diff.0.2) %in% short.colnames])
dbmts.agree.best.0.2.short <- unique(dbmts.agree.best.0.2[colnames(dbmts.agree.best.0.2) %in% short.colnames])
dbmts.agree.best.best.0.2.short <- unique(dbmts.agree.best.best.0.2[colnames(dbmts.agree.best.best.0.2) %in% short.colnames])

# Create unique variant-score data frame based on the best best-vs-best scores
short.var.score.colnames <- c("variant", "diffScore", "gtex.eqtl", "gtex.brain.eqtl", "psychencode.eqtl", "any.eqtl")
dbmts.agree.best.best.short.var.score <- unique(dbmts.agree.best.best.short[colnames(dbmts.agree.best.best.short) %in% short.var.score.colnames])
dbmts.agree.best.best.0.2.short.var.score <- unique(dbmts.agree.best.best.0.2.short[colnames(dbmts.agree.best.best.0.2.short) %in% short.var.score.colnames])
eqtl_or <- function(x, df, eqtl_column) {
  col_select = (colnames(df) == eqtl_column)
  eqtl = df[df$variant == x, col_select]
  return(TRUE %in% eqtl)
}
dbmts.agree.best.best.short.var.score$gtex.eqtl <- do.call('c', lapply(dbmts.agree.best.best.short.var.score$variant, eqtl_or, df = dbmts.agree.best.best.short.var.score, eqtl_column = "gtex.eqtl"))
dbmts.agree.best.best.short.var.score$gtex.brain.eqtl <- do.call('c', lapply(dbmts.agree.best.best.short.var.score$variant, eqtl_or, df = dbmts.agree.best.best.short.var.score, eqtl_column = "gtex.brain.eqtl"))
dbmts.agree.best.best.short.var.score$psychencode.eqtl <- do.call('c', lapply(dbmts.agree.best.best.short.var.score$variant, eqtl_or, df = dbmts.agree.best.best.short.var.score, eqtl_column = "psychencode.eqtl"))
dbmts.agree.best.best.short.var.score$any.eqtl <- do.call('c', lapply(dbmts.agree.best.best.short.var.score$variant, eqtl_or, df = dbmts.agree.best.best.short.var.score, eqtl_column = "any.eqtl"))
dbmts.agree.best.best.short.var.score <- unique(dbmts.agree.best.best.short.var.score)
dbmts.agree.best.best.0.2.short.var.score$gtex.eqtl <- do.call('c', lapply(dbmts.agree.best.best.0.2.short.var.score$variant, eqtl_or, df = dbmts.agree.best.best.0.2.short.var.score, eqtl_column = "gtex.eqtl"))
dbmts.agree.best.best.0.2.short.var.score$gtex.brain.eqtl <- do.call('c', lapply(dbmts.agree.best.best.0.2.short.var.score$variant, eqtl_or, df = dbmts.agree.best.best.0.2.short.var.score, eqtl_column = "gtex.brain.eqtl"))
dbmts.agree.best.best.0.2.short.var.score$psychencode.eqtl <- do.call('c', lapply(dbmts.agree.best.best.0.2.short.var.score$variant, eqtl_or, df = dbmts.agree.best.best.0.2.short.var.score, eqtl_column = "psychencode.eqtl"))
dbmts.agree.best.best.0.2.short.var.score$any.eqtl <- do.call('c', lapply(dbmts.agree.best.best.0.2.short.var.score$variant, eqtl_or, df = dbmts.agree.best.best.0.2.short.var.score, eqtl_column = "any.eqtl"))
dbmts.agree.best.best.0.2.short.var.score <- unique(dbmts.agree.best.best.0.2.short.var.score)

# Load miRNA families
mirfams_file <- '../db/mirfams.txt'
mirfams <- read.delim(mirfams_file, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
colnames(mirfams) <- c('mature', 'family')

# Add miRNA family information
dbmts.agree.diff.short$mirna.family <- mirfams$family[match(dbmts.agree.diff.short$mirna, mirfams$mature)]
dbmts.agree.best.short$mirna.family <- mirfams$family[match(dbmts.agree.best.short$mirna, mirfams$mature)]
dbmts.agree.best.best.short$mirna.family <- mirfams$family[match(dbmts.agree.best.best.short$mirna, mirfams$mature)]
dbmts.agree.diff.0.2.short$mirna.family <- mirfams$family[match(dbmts.agree.diff.0.2.short$mirna, mirfams$mature)]
dbmts.agree.best.0.2.short$mirna.family <- mirfams$family[match(dbmts.agree.best.0.2.short$mirna, mirfams$mature)]
dbmts.agree.best.best.0.2.short$mirna.family <- mirfams$family[match(dbmts.agree.best.best.0.2.short$mirna, mirfams$mature)]

# Add dbSNP rsID to data
# Read in conversion file
convert_file <- '../db/g1000_eur.magma.convert'
convert <- read.delim(convert_file, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
library(reshape2)
convert <- melt(convert, id.vars = "V1")
colnames(convert) <- c("rsID", "variable", "variant")
convert <- convert[colnames(convert) %in% c("rsID", "variant")]
convert <- convert[!convert$variant == "",]
convert <- unique(convert)
# Add dbSNP column
dbmts.agree.diff.short$dbsnp <- convert$rsID[match(dbmts.agree.diff.short$variant, convert$variant)]
dbmts.agree.best.short$dbsnp <- convert$rsID[match(dbmts.agree.best.short$variant, convert$variant)]
dbmts.agree.best.best.short$dbsnp <- convert$rsID[match(dbmts.agree.best.best.short$variant, convert$variant)]
dbmts.agree.diff.0.2.short$dbsnp <- convert$rsID[match(dbmts.agree.diff.0.2.short$variant, convert$variant)]
dbmts.agree.best.0.2.short$dbsnp <- convert$rsID[match(dbmts.agree.best.0.2.short$variant, convert$variant)]
dbmts.agree.best.best.0.2.short$dbsnp <- convert$rsID[match(dbmts.agree.best.best.0.2.short$variant, convert$variant)]
dbmts.agree.best.best.short.var.score$dbsnp <- convert$rsID[match(dbmts.agree.best.best.short.var.score$variant, convert$variant)]
dbmts.agree.best.best.0.2.short.var.score$dbsnp <- convert$rsID[match(dbmts.agree.best.best.0.2.short.var.score$variant, convert$variant)]

# Add MHC information
mhc_file <- '../db/g1000_eur.mhc.snps'
mhc <- scan(mhc_file, character())

# Add non.MHC column
dbmts.agree.diff.short$non.mhc <- !(dbmts.agree.diff.short$dbsnp %in% mhc)
dbmts.agree.diff.short$non.mhc[is.na(dbmts.agree.diff.short$dbsnp)] <- NA

dbmts.agree.best.short$non.mhc <- !(dbmts.agree.best.short$dbsnp %in% mhc)
dbmts.agree.best.short$non.mhc[is.na(dbmts.agree.best.short$dbsnp)] <- NA

dbmts.agree.best.best.short$non.mhc <- !(dbmts.agree.best.best.short$dbsnp %in% mhc)
dbmts.agree.best.best.short$non.mhc[is.na(dbmts.agree.best.best.short$dbsnp)] <- NA

dbmts.agree.diff.0.2.short$non.mhc <- !(dbmts.agree.diff.0.2.short$dbsnp %in% mhc)
dbmts.agree.diff.0.2.short$non.mhc[is.na(dbmts.agree.diff.0.2.short$dbsnp)] <- NA

dbmts.agree.best.0.2.short$non.mhc <- !(dbmts.agree.best.0.2.short$dbsnp %in% mhc)
dbmts.agree.best.0.2.short$non.mhc[is.na(dbmts.agree.best.0.2.short$dbsnp)] <- NA

dbmts.agree.best.best.0.2.short$non.mhc <- !(dbmts.agree.best.best.0.2.short$dbsnp %in% mhc)
dbmts.agree.best.best.0.2.short$non.mhc[is.na(dbmts.agree.best.best.0.2.short$dbsnp)] <- NA

dbmts.agree.best.best.short.var.score$non.mhc <- !(dbmts.agree.best.best.short.var.score$dbsnp %in% mhc)
dbmts.agree.best.best.short.var.score$non.mhc[is.na(dbmts.agree.best.best.short.var.score$dbsnp)] <- NA

dbmts.agree.best.best.0.2.short.var.score$non.mhc <- !(dbmts.agree.best.best.0.2.short.var.score$dbsnp %in% mhc)
dbmts.agree.best.best.0.2.short.var.score$non.mhc[is.na(dbmts.agree.best.best.0.2.short.var.score$dbsnp)] <- NA

# Save dbmts data to disk
save.image("mbsv.analysis.1.RData")
save(dbmts.agree.diff.short,
     dbmts.agree.diff,
     dbmts.agree.best.short,
     dbmts.agree.best,
     dbmts.agree.best.best.short,
     dbmts.agree.best.best,
     dbmts.agree.diff.0.2.short,
     dbmts.agree.diff.0.2,
     dbmts.agree.best.0.2.short,
     dbmts.agree.best.0.2,
     dbmts.agree.best.best.0.2.short,
     dbmts.agree.best.best.0.2,
     dbmts.agree.best.best.short.var.score,
     dbmts.agree.best.best.0.2.short.var.score,
     file = "dbmts.RData")
save(convert, file = "convert.RData")
