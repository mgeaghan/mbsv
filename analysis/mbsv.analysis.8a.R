# Get number of genome-wide and suggestive-significant SNPs per GWAS and compare MBSV subsets with Fisher's exact test. Compare MBSVs against non-MBSV 3'UTR variants.
# Author: Michael Geaghan (2020)

# Get list of brain-expressed genes
gene_list <- scan("../db/filteredgenes.txt", character())

# Download gene and transcript IDs and transcript types from biomart (GRCh37)
library(biomaRt)
ensembl <- useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl', host = 'http://grch37.ensembl.org')
bm <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_gene_id_version', 'ensembl_transcript_id_version', 'transcript_biotype', 'chromosome_name', '3_utr_start', '3_utr_end'), mart = ensembl)
chrs <- c(seq(1, 22))
bm.all <- bm
bm <- bm[bm$transcript_biotype=='protein_coding' & bm$chromosome_name %in% chrs & !is.na(bm$`3_utr_start`) & !is.na(bm$`3_utr_end`),]
bm <- bm[bm$ensembl_gene_id %in% gene_list,]
utr_3_regions <- bm[c("chromosome_name", "3_utr_start", "3_utr_end")]
colnames(utr_3_regions) <- c("chr", "start", "end")
utr_3_regions <- unique(utr_3_regions)
utr_3_regions$diff <- utr_3_regions$end - utr_3_regions$start
utr_3_regions.list <- list()
for (c in chrs) {
  utr_3_regions.list[[c]] <- utr_3_regions[as.numeric(utr_3_regions$chr) == c,]
}

gwas_list <- c('adhd.2018', 'an.2019', 'asd.2019', 'bip.2018',
               'mdd.2019', 'ocd.2018', 'ptsd.2019',
               'scz.clozuk2018', 'ts.2019')

getNumVar <- function(x, d) {
  if (x == "UTR3") {
    return(dim(d[d$variable == x & d$non.mhc,])[1])
  } else {
    return(dim(d[d$variable == x & d$non.mhc & abs(d$diffScore) >= 0.2,])[1])
  }
}
getNumSigVar <- function(x, d, alpha) {
  if (x == "UTR3") {
    return(dim(d[d$variable == x & d$non.mhc & d$P < alpha,])[1])
  } else {
    return(dim(d[d$variable == x & d$non.mhc & abs(d$diffScore) >= 0.2 & d$P < alpha,])[1])
  }
}

for (threshold in c(5e-8, 1e-5)) {
  df <- data.frame(gwas = character(), num.var = numeric(), num.var.sig = numeric(), prop.var.sig = numeric(),
                   num.mbsv = numeric(), num.mbsv.sig = numeric(), prop.mbsv.sig = numeric(), or.mbsv.sig = numeric(), p.mbsv.sig = numeric(), fdr.mbsv.sig = numeric(),
                   num.gtex.eqtl.mbsv = numeric(), num.gtex.eqtl.mbsv.sig = numeric(), prop.gtex.eqtl.mbsv.sig = numeric(), or.gtex.eqtl.mbsv.sig = numeric(), p.gtex.eqtl.mbsv.sig = numeric(), fdr.gtex.eqtl.mbsv.sig = numeric(),
                   num.gtex.brain.eqtl.mbsv = numeric(), num.gtex.brain.eqtl.mbsv.sig = numeric(), prop.gtex.brain.eqtl.mbsv.sig = numeric(), or.gtex.brain.eqtl.mbsv.sig = numeric(), p.gtex.brain.eqtl.mbsv.sig = numeric(), fdr.gtex.brain.eqtl.mbsv.sig = numeric(),
                   num.psych.eqtl.mbsv = numeric(), num.psych.eqtl.mbsv.sig = numeric(), prop.psych.eqtl.mbsv.sig = numeric(), or.psych.eqtl.mbsv.sig = numeric(), p.psych.eqtl.mbsv.sig = numeric(), fdr.psych.eqtl.mbsv.sig = numeric())
  for (g in gwas_list) {
    gwas.pruned <- paste(g, ".gwas.nonmhc.pruned.RData", sep = "")
    load(gwas.pruned)
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
    gwas.df.utr3 <- gwas.df[gwas.df$variable == "GWAS" & gwas.df$SNP %in% gwas.scb.utr.snp,]
    gwas.df.utr3$variable <- "UTR3"
    gwas.df <- rbind(gwas.df, gwas.df.utr3)
    gwas.df <- gwas.df[gwas.df$variable != "GWAS",]
    gwas.df$variable <- as.character(gwas.df$variable)
    gwas.df <- gwas.df[gwas.df$SNP %in% gwas.df$SNP[gwas.df$variable=="UTR3"],]
    # GWAS vars
    num.var <- getNumVar("UTR3", gwas.df)
    num.var.sig <- getNumSigVar("UTR3", gwas.df, threshold)
    prop.var.sig <- num.var.sig/num.var
    # MBSVs
    num.mbsv <- getNumVar("mbsv", gwas.df)
    num.mbsv.sig <- getNumSigVar("mbsv", gwas.df, threshold)
    num.nonmbsv <- num.var - num.mbsv
    num.nonmbsv.sig <- num.var.sig - num.mbsv.sig
    prop.mbsv.sig <- num.mbsv.sig/num.mbsv
    ft <- fisher.test(data.frame(non.mbsv = c(num.nonmbsv, num.nonmbsv.sig), mbsv = c(num.mbsv, num.mbsv.sig), row.names = c("all", "mbsv")))
    or.mbsv.sig <- ft$estimate
    p.mbsv.sig <- ft$p.value
    # eQTLs
    num.gtex.eqtl.mbsv <- getNumVar("mbsv.gtex.eqtl", gwas.df)
    num.gtex.eqtl.mbsv.sig <- getNumSigVar("mbsv.gtex.eqtl", gwas.df, threshold)
    num.gtex.eqtl.nonmbsv <- num.var - num.mbsv
    num.gtex.eqtl.nonmbsv.sig <- num.var.sig - num.gtex.eqtl.mbsv.sig
    prop.gtex.eqtl.mbsv.sig <- num.gtex.eqtl.mbsv.sig/num.gtex.eqtl.mbsv
    ft <- fisher.test(data.frame(non.mbsv = c(num.gtex.eqtl.nonmbsv, num.gtex.eqtl.nonmbsv.sig), mbsv = c(num.gtex.eqtl.mbsv, num.gtex.eqtl.mbsv.sig), row.names = c("all", "mbsv")))
    or.gtex.eqtl.mbsv.sig <- ft$estimate
    p.gtex.eqtl.mbsv.sig <- ft$p.value
    # Brain eQTLs
    num.gtex.brain.eqtl.mbsv <- getNumVar("mbsv.gtex.brain.eqtl", gwas.df)
    num.gtex.brain.eqtl.mbsv.sig <- getNumSigVar("mbsv.gtex.brain.eqtl", gwas.df, threshold)
    num.gtex.brain.eqtl.nonmbsv <- num.var - num.gtex.brain.eqtl.mbsv
    num.gtex.brain.eqtl.nonmbsv.sig <- num.var.sig - num.gtex.brain.eqtl.mbsv.sig
    prop.gtex.brain.eqtl.mbsv.sig <- num.gtex.brain.eqtl.mbsv.sig/num.gtex.brain.eqtl.mbsv
    ft <- fisher.test(data.frame(non.mbsv = c(num.gtex.brain.eqtl.nonmbsv, num.gtex.brain.eqtl.nonmbsv.sig), mbsv = c(num.gtex.brain.eqtl.mbsv, num.gtex.brain.eqtl.mbsv.sig), row.names = c("all", "mbsv")))
    or.gtex.brain.eqtl.mbsv.sig <- ft$estimate
    p.gtex.brain.eqtl.mbsv.sig <- ft$p.value
    # Psych eQTLs
    num.psych.eqtl.mbsv <- getNumVar("mbsv.psychencode.eqtl", gwas.df)
    num.psych.eqtl.mbsv.sig <- getNumSigVar("mbsv.psychencode.eqtl", gwas.df, threshold)
    num.psych.eqtl.nonmbsv <- num.var - num.psych.eqtl.mbsv
    num.psych.eqtl.nonmbsv.sig <- num.var.sig - num.psych.eqtl.mbsv.sig
    prop.psych.eqtl.mbsv.sig <- num.psych.eqtl.mbsv.sig/num.psych.eqtl.mbsv
    ft <- fisher.test(data.frame(non.mbsv = c(num.psych.eqtl.nonmbsv, num.psych.eqtl.nonmbsv.sig), mbsv = c(num.psych.eqtl.mbsv, num.psych.eqtl.mbsv.sig), row.names = c("all", "mbsv")))
    or.psych.eqtl.mbsv.sig <- ft$estimate
    p.psych.eqtl.mbsv.sig <- ft$p.value
    p <- c(p.mbsv.sig, p.gtex.eqtl.mbsv.sig, p.gtex.brain.eqtl.mbsv.sig, p.psych.eqtl.mbsv.sig)
    q <- p.adjust(p, method = "BH")
    fdr.mbsv.sig <- q[1]
    fdr.gtex.eqtl.mbsv.sig <- q[2]
    fdr.gtex.brain.eqtl.mbsv.sig <- q[3]
    fdr.psych.eqtl.mbsv.sig <- q[4]
    # Add to DF
    rm(gwas, gwas.df)
    df <- rbind(df, data.frame(gwas = g, num.var = num.var, num.var.sig = num.var.sig, prop.var.sig = prop.var.sig,
                               num.mbsv = num.mbsv, num.mbsv.sig = num.mbsv.sig, prop.mbsv.sig = prop.mbsv.sig, or.mbsv.sig = or.mbsv.sig, p.mbsv.sig = p.mbsv.sig, fdr.mbsv.sig = fdr.mbsv.sig,
                               num.gtex.eqtl.mbsv = num.gtex.eqtl.mbsv, num.gtex.eqtl.mbsv.sig = num.gtex.eqtl.mbsv.sig, prop.gtex.eqtl.mbsv.sig = prop.gtex.eqtl.mbsv.sig, or.gtex.eqtl.mbsv.sig = or.gtex.eqtl.mbsv.sig, p.gtex.eqtl.mbsv.sig = p.gtex.eqtl.mbsv.sig, fdr.gtex.eqtl.mbsv.sig = fdr.gtex.eqtl.mbsv.sig,
                               num.gtex.brain.eqtl.mbsv = num.gtex.brain.eqtl.mbsv, num.gtex.brain.eqtl.mbsv.sig = num.gtex.brain.eqtl.mbsv.sig, prop.gtex.brain.eqtl.mbsv.sig = prop.gtex.brain.eqtl.mbsv.sig, or.gtex.brain.eqtl.mbsv.sig = or.gtex.brain.eqtl.mbsv.sig, p.gtex.brain.eqtl.mbsv.sig = p.gtex.brain.eqtl.mbsv.sig, fdr.gtex.brain.eqtl.mbsv.sig = fdr.gtex.brain.eqtl.mbsv.sig,
                               num.psych.eqtl.mbsv = num.psych.eqtl.mbsv, num.psych.eqtl.mbsv.sig = num.psych.eqtl.mbsv.sig, prop.psych.eqtl.mbsv.sig = prop.psych.eqtl.mbsv.sig, or.psych.eqtl.mbsv.sig = or.psych.eqtl.mbsv.sig, p.psych.eqtl.mbsv.sig = p.psych.eqtl.mbsv.sig, fdr.psych.eqtl.mbsv.sig = fdr.psych.eqtl.mbsv.sig))
  }
  rownames(df) <- df$gwas
  numVar.df <- df
  save(numVar.df, file = paste("gwas.num_sig_vars.", as.character(threshold), ".nonmbsvbrain3utr.RData", sep = ""))
  write.table(numVar.df, file = paste("gwas.num_sig_vars.", as.character(threshold), ".nonmbsvbrain3utr.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
