# Analyse MAGMA meta analysis results (all variants - no p-value threshold)
# Author: Michael Geaghan (2020)

# Gene set-level

library(reshape2)
setwd("magma/meta/")
metas <- c('scz.bip', 'scz.mdd', 'scz.bip.mdd', 'adhd.an.asd.bip.mdd.ocd.ptsd.scz.ts')
subsets <- c('mbsv.0.2', paste('mbsv.0.2', c('eqtl', 'brain.eqtl', 'psych.eqtl'), sep = "."))

gene_sets <- list()
gene_sets_merged <- list()
for (m in metas) {
  setwd(m)
  gene_sets_merged[[m]] <- data.frame(VARIABLE = character(), TYPE = character(), NGENES = numeric(), BETA = numeric(), BETA_STD = numeric(), SE = numeric(), P = numeric(), SUBSET = character())
  for (s in subsets) {
    n = paste(m, s, sep = "_")
    file <- paste(m, 'magma.all', s, 'sets.go.gsa.out', sep = '.')
    gene_sets[[n]] <- read.delim(file, header = TRUE, comment.char = "#", sep = "", stringsAsFactors = FALSE, row.names = "VARIABLE")
    gene_sets[[n]]$SUBSET <- n
    gene_sets_merged[[m]] <- rbind(gene_sets_merged[[m]], gene_sets[[n]])
  }
  gene_sets_merged[[m]]$BONF.P <- p.adjust(gene_sets_merged[[m]]$P, method = "bonferroni")
  gene_sets_merged[[m]]$FDR <- p.adjust(gene_sets_merged[[m]]$P, method = "BH")
  setwd("../")
}

# only calculate for sets with 5 genes or more
gene_sets_5 <- gene_sets
for (n in names(gene_sets_5)) {
  gene_sets_5[[n]] <- gene_sets_5[[n]][gene_sets_5[[n]]$NGENES >= 5,]
  gene_sets_5[[n]]$FDR <- p.adjust(gene_sets_5[[n]]$P, method = "BH")
  gene_sets_5[[n]]$BONFERRONI <- p.adjust(gene_sets_5[[n]]$P, method = "bonferroni")
}
gene_sets_merged_5 <- gene_sets_merged
for (n in names(gene_sets_merged_5)) {
  gene_sets_merged_5[[n]] <- gene_sets_merged_5[[n]][gene_sets_merged_5[[n]]$NGENES >= 5,]
  gene_sets_merged_5[[n]]$FDR <- p.adjust(gene_sets_merged_5[[n]]$P, method = "BH")
  gene_sets_merged_5[[n]]$BONFERRONI <- p.adjust(gene_sets_merged_5[[n]]$P, method = "bonferroni")
}

save.image("magma.meta.all.analyses.RData")

cols <- c("FULL_NAME", "SUBSET", "NGENES", "BETA", "BETA_STD", "SE", "P", "BONF.P", 'FDR')

# for (m in names(gene_sets_merged_5)) {
#   name = paste(m, "magma.meta.all.gene_set_analysis.txt", sep = ".")
#   write.table(gene_sets_merged_5[[m]][cols], name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
# }

for (m in names(gene_sets_merged)) {
  name = paste(m, "magma.meta.all.gene_set_analysis.txt", sep = ".")
  write.table(gene_sets_merged[[m]][cols], name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# Gene-level

meta_genes <- list()
meta_genes_merged <- list()
for (m in metas) {
  setwd(m)
  meta_genes_merged[[m]] <- data.frame(GENE = character(), CHR = character(), START = numeric(), STOP = numeric(), NSNPS = numeric(), NPARAM = numeric(), N = numeric(), ZSTAT = numeric(), P = numeric(), SUBSET = character())
  for (s in subsets) {
    n = paste(m, s, sep = "_")
    file <- paste(m, 'magma.all', s, 'out.genes.out', sep = '.')
    meta_genes[[n]] <- read.delim(file, header = TRUE, comment.char = "#", sep = "", stringsAsFactors = FALSE)
    rownames(meta_genes[[n]]) <- meta_genes[[n]]$GENE
    meta_genes[[n]]$SUBSET <- n
    meta_genes_merged[[m]] <- rbind(meta_genes_merged[[m]], meta_genes[[n]])
  }
  meta_genes_merged[[m]]$BONF.P <- p.adjust(meta_genes_merged[[m]]$P, method = "bonferroni")
  meta_genes_merged[[m]]$FDR <- p.adjust(meta_genes_merged[[m]]$P, method = "BH")
  setwd("../")
}

save.image("magma.meta.all.gene.analyses.RData")

cols <- c("GENE", "SUBSET", "CHR", "START", "STOP", "NSNPS", "NPARAM", "N", "DATASETS", "ZSTAT", "P", "BONF.P")

for (m in names(meta_genes_merged)) {
  name = paste(m, "magma.meta.all.gene_analysis.txt", sep = ".")
  write.table(meta_genes_merged[[m]][cols], name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
