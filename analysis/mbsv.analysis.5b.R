# Perform ECDF comparisons for P-values, absolute effect sizes, CADD scores and European MAFs. Compare MBSVs against non-MBSV 3'UTR variants.
# Author: Michael Geaghan (2020)

library(ggplot2)

# Read in GWAS data
gwas_files <- list(adhd.2018 = './adhd.2018.gwas.nonmhc.pruned.RData',
                   an.2019 = './an.2019.gwas.nonmhc.pruned.RData',
                   asd.2019 = './asd.2019.gwas.nonmhc.pruned.RData',
                   bip.2018 = './bip.2018.gwas.nonmhc.pruned.RData',
                   mdd.2019 = './mdd.2019.gwas.nonmhc.pruned.RData',
                   ocd.2018 = './ocd.2018.gwas.nonmhc.pruned.RData',
                   ptsd.2019 = './ptsd.2019.gwas.nonmhc.pruned.RData',
                   scz.clozuk2018 = './scz.clozuk2018.gwas.nonmhc.pruned.RData',
                   ts.2019 = './ts.2019.gwas.nonmhc.pruned.RData')

gwas.union.df <- data.frame(SNP = character(), variant = character(), diffScore = numeric(), non.mhc = logical(), variable = character(), value = logical())

# Download gene and transcript IDs and transcript types from biomart (GRCh37)
library(biomaRt)
ensembl <- useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl', host = 'http://grch37.ensembl.org')
bm <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_gene_id_version', 'ensembl_transcript_id_version', 'transcript_biotype', 'chromosome_name', '3_utr_start', '3_utr_end'), mart = ensembl)
chrs <- c(seq(1, 22))
bm.all <- bm
bm <- bm[bm$transcript_biotype=='protein_coding' & bm$chromosome_name %in% chrs & !is.na(bm$`3_utr_start`) & !is.na(bm$`3_utr_end`),]
utr_3_regions <- bm[c("chromosome_name", "3_utr_start", "3_utr_end")]
colnames(utr_3_regions) <- c("chr", "start", "end")
utr_3_regions <- unique(utr_3_regions)
utr_3_regions$diff <- utr_3_regions$end - utr_3_regions$start
utr_3_regions.list <- list()
for (c in chrs) {
  utr_3_regions.list[[c]] <- utr_3_regions[as.numeric(utr_3_regions$chr) == c,]
}

# Functions to generate cumulative density plots
ecdf.p <- function(df) {
  ret <- ggplot(data = df, aes(x = P, color = variable))
  ret <- ret + stat_ecdf(geom = "step") + labs(x = "p-value", y = "cumulative fraction", color = "variant subset")
  return(ret)
}

ecdf.es <- function(df, es) {
  lab.x <- paste("absolute", es, sep = " ")
  ret <- ggplot(data = df, aes(x = ABS.ES, color = variable))
  ret <- ret + stat_ecdf(geom = "step") + labs(x = lab.x, y = "cumulative fraction", color = "variant subset")
  return(ret)
}

ecdf.cadd <- function(df) {
  ret <- ggplot(data = df, aes(x = cadd_rank, color = variable))
  ret <- ret + stat_ecdf(geom = "step") + labs(x = "CADD rank (v1.6)", y = "cumulative fraction", color = "variant subset")
  return(ret)
}

ecdf.cadd.v1 <- function(df) {
  ret <- ggplot(data = df, aes(x = cadd_rank_v1, color = variable))
  ret <- ret + stat_ecdf(geom = "step") + labs(x = "CADD rank (v1.0)", y = "cumulative fraction", color = "variant subset")
  return(ret)
}

ecdf.maf <- function(df) {
  ret <- ggplot(data = df, aes(x = maf.eur, color = variable))
  ret <- ret + stat_ecdf(geom = "step") + labs(x = "minor allele frequency (Europeans)", y = "cumulative fraction", color = "variant subset")
  return(ret)
}

mbsv.ks <- function(df, header, ref) {
  vars <- levels(df$variable)
  vars <- vars[vars != ref]
  ref.med <- median(df[[header]][df$variable == ref], na.rm = TRUE)
  ret <- data.frame(ref.med = ref.med)
  colnames(ret) <- c(paste(gsub(" ", ".", ref), ".median", sep = ""))
  for (v in vars) {
    v.ks <- ks.test(df[[header]][df$variable == v], df[[header]][df$variable == ref])
    v.med <- median(df[[header]][df$variable == v], na.rm = TRUE)
    ret[[paste(gsub(" ", ".", v), ".median", sep = "")]] <- v.med
    ret[[paste(gsub(" ", ".", v), ".ks.d", sep = "")]] <- v.ks$statistic
    ret[[paste(gsub(" ", ".", v), ".ks.p", sep = "")]] <- v.ks$p.value
    ret[[paste(gsub(" ", ".", v), ".ks.fdr", sep = "")]] <- NA
  }
  p <- list()
  for (v in vars) {
    p[[v]] <- ret[[paste(gsub(" ", ".", v), ".ks.p", sep = "")]]
  }
  fdr <- p.adjust(p, method = "BH")
  for (v in vars) {
    ret[[paste(gsub(" ", ".", v), ".ks.fdr", sep = "")]] <- fdr[[v]]
  }
  return(ret)
}

# P-value and ES ECDF comparisons
for (g in names(gwas_files)) {
  load(gwas_files[[g]])
  g2 <- gsub(".", "-", g, fixed = TRUE)
  gwas.file <- paste("../gwas/", g2, ".gwas", sep = "")
  gwas.scb <- read.delim(gwas.file)
  gwas.scb.cols <- colnames(gwas.scb)[colnames(gwas.scb) %in% c("SNP", "CHR", "CHROM", "BP", "POS")]
  gwas.scb <- gwas.scb[gwas.scb.cols]
  colnames(gwas.scb) <- c("SNP", "CHR", "BP")
  gwas.df <- gwas.df[!is.na(gwas.df$non.mhc) &
                       !is.na(gwas.df$variable) &
                       gwas.df$non.mhc &
                       gwas.df$variable != "mbsv.any.eqtl" &
                       (gwas.df$variable == "GWAS" | (!is.na(gwas.df$diffScore) & abs(gwas.df$diffScore) >= 0.2)),]
  gwas.scb <- gwas.scb[gwas.scb$SNP %in% unique(gwas.df$SNP),]
  gwas.scb$UTR3 <- apply(gwas.scb, 1, function(x) {
    c <- as.numeric(x[2])
    b <- as.numeric(x[3])
    d <- b - utr_3_regions.list[[c]]$start
    d2 <- utr_3_regions.list[[c]]$diff
    return(any((d >= 0) & (d <= d2)))
  })
  gwas.scb.utr.snp <- gwas.scb$SNP[gwas.scb$UTR3]
  mbsv.snp <- unique(gwas.df$SNP[gwas.df$variable == "mbsv"])
  gwas.df.nonmbsv <- gwas.df[gwas.df$variable == "GWAS" & !(gwas.df$SNP %in% mbsv.snp) & gwas.df$SNP %in% gwas.scb.utr.snp,]
  gwas.df.nonmbsv$variable <- "utr3.non.mbsv"
  gwas.df <- rbind(gwas.df, gwas.df.nonmbsv)
  gwas.df <- gwas.df[gwas.df$variable != "GWAS",]
  gwas.df$variable <- as.character(gwas.df$variable)
  # gwas.df$variable[gwas.df$variable == "GWAS"] <- "all GWAS variants"
  gwas.df$variable[gwas.df$variable == "utr3.non.mbsv"] <- "non-MBSV 3'UTR variants"
  gwas.df$variable[gwas.df$variable == "mbsv"] <- "all MBSVs"
  gwas.df$variable[gwas.df$variable == "mbsv.gtex.eqtl"] <- "GTEx eQTL MBSVs"
  gwas.df$variable[gwas.df$variable == "mbsv.gtex.brain.eqtl"] <- "GTEx Brain eQTL MBSVs"
  gwas.df$variable[gwas.df$variable == "mbsv.psychencode.eqtl"] <- "PSYCHENCODE eQTL MBSVs"
  gwas.df$ABS.ES = abs(gwas.df$ES)
  # gwas.df$variable = factor(gwas.df$variable, levels = c("all GWAS variants", "non-MBSV GWAS variants", "all MBSVs", "GTEx eQTL MBSVs", "GTEx Brain eQTL MBSVs", "PSYCHENCODE eQTL MBSVs"))
  gwas.df$variable = factor(gwas.df$variable, levels = c("non-MBSV 3'UTR variants", "all MBSVs", "GTEx eQTL MBSVs", "GTEx Brain eQTL MBSVs", "PSYCHENCODE eQTL MBSVs"))
  ggsave(paste(g, ".ecdf.p.mbsv.0.2.nonmhc.nonmbsv3utr.pruned.png", sep = ""), ecdf.p(gwas.df), device = "png", width = 16, height = 9, dpi = 600)
  ggsave(paste(g, ".ecdf.es.mbsv.0.2.nonmhc.nonmbsv3utr.pruned.png", sep = ""), ecdf.es(gwas.df, gwas_es_char), device = "png", width = 16, height = 9, dpi = 600)
  if ("ks.p.df" %in% ls()) {
    ks.p.df.rows <- rownames(ks.p.df)
    ks.p.df <- rbind(ks.p.df, mbsv.ks(gwas.df, "P", "non-MBSV 3'UTR variants"))
    rownames(ks.p.df) <- c(ks.p.df.rows, g)
    ks.es.df.rows <- rownames(ks.es.df)
    ks.es.df <- rbind(ks.es.df, mbsv.ks(gwas.df, "ABS.ES", "non-MBSV 3'UTR variants"))
    rownames(ks.es.df) <- c(ks.es.df.rows, g)
  } else {
    ks.p.df <- mbsv.ks(gwas.df, "P", "non-MBSV 3'UTR variants")
    rownames(ks.p.df) <- g
    ks.es.df <- mbsv.ks(gwas.df, "ABS.ES", "non-MBSV 3'UTR variants")
    rownames(ks.es.df) <- g
  }
  # Combine current gwas into the gwas.union.df data.frame
  gwas.union.df <- unique(rbind(gwas.union.df, gwas.df[c("SNP", "variant", "diffScore", "non.mhc", "variable", "value")]))
}

# CADD and MAF ECDF comparisons
# load annotation data (CADD v1.6)
annot <- read.delim("gwas.nonmhc.annot.snv.dbsnp.txt", header = TRUE, sep = "\t", na.strings = ".", stringsAsFactors = FALSE)
annot <- annot[!is.na(annot$dbsnp),]
annot_cadd <- data.frame(PHRED = annot$PHRED, row.names = annot$dbsnp)
annot_af <- data.frame(AF_nfe_exome = annot$AF_nfe_exome, AF_nfe_genome = annot$AF_nfe_genome, row.names = annot$dbsnp)
annot_af$af <- as.numeric(apply(annot_af, 1, function(x) {
  eAf = as.numeric(x[1])
  gAf = as.numeric(x[2])
  if(is.na(eAf)) {
    if(is.na(gAf)) {
      # return NA
      return(NA)
    } else {
      # return gAf
      return(gAf)
    }
  } else {
    # return eAf (more samples)
    return(eAf)
  }
}))
annot_af$maf <- sapply(annot_af$af, function(x) {
  if (is.na(x)) {
    return(x)
  } else if (x <= 0.5) {
    return(x)
  } else {
    return(1 - x)
  }
})
gwas_vars <- gwas.union.df$SNP
gwas.union.df$maf.eur <- annot_af[gwas_vars,]$maf
gwas.union.df$cadd_phred <- annot_cadd[gwas_vars,]
gwas.union.df$cadd_rank <- 1 - (10^(gwas.union.df$cadd_phred/(-10)))
# load annotation data (CADD v1.0)
annot <- read.delim("gwas.nonmhc.annot.v1.snv.dbsnp.txt", header = TRUE, sep = "\t", na.strings = ".", stringsAsFactors = FALSE)
annot <- annot[!is.na(annot$dbsnp),]
annot_cadd <- data.frame(PHRED = annot$PHRED, row.names = annot$dbsnp)
gwas.union.df$cadd_phred_v1 <- annot_cadd[gwas_vars,]
gwas.union.df$cadd_rank_v1 <- 1 - (10^(gwas.union.df$cadd_phred_v1/(-10)))
# Plots and K-S statistics
ggsave("gwas.ecdf.cadd.mbsv.0.2.nonmhc.nonmbsv3utr.pruned.png", ecdf.cadd(gwas.union.df), device = "png", width = 16, height = 9, dpi = 600)
ggsave("gwas.ecdf.cadd.v1.mbsv.0.2.nonmhc.nonmbsv3utr.pruned.png", ecdf.cadd.v1(gwas.union.df), device = "png", width = 16, height = 9, dpi = 600)
ggsave("gwas.ecdf.maf.eur.mbsv.0.2.nonmhc.nonmbsv3utr.pruned.png", ecdf.maf(gwas.union.df), device = "png", width = 16, height = 9, dpi = 600)
ks.cadd.df <- mbsv.ks(gwas.union.df, "cadd_rank", "non-MBSV 3'UTR variants")
rownames(ks.cadd.df) <- "CADD Rank (v1.6)"
ks.cadd.v1.df <- mbsv.ks(gwas.union.df, "cadd_rank_v1", "non-MBSV 3'UTR variants")
rownames(ks.cadd.v1.df) <- "CADD Rank (v1.0)"
ks.maf.eur.df <- mbsv.ks(gwas.union.df, "maf.eur", "non-MBSV 3'UTR variants")
rownames(ks.maf.eur.df) <- "European MAF"

# Write tables
write.table(ks.p.df, "ks.p.nonmbsv3utr.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(ks.es.df, "ks.es.nonmbsv3utr.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(ks.cadd.df, "ks.cadd.nonmbsv3utr.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(ks.cadd.v1.df, "ks.cadd.v1.nonmbsv3utr.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(ks.maf.eur.df, "ks.maf.eur.nonmbsv3utr.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
