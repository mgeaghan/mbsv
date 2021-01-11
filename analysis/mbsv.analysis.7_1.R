# Analyse MAGMA results (all variants - no p-value threshold)
# Author: Michael Geaghan (2020)

# Gene set-level

library(reshape2)
setwd("magma/")
gwas <- c('adhd.2018', 'an.2019', 'asd.2019', 'bip.2018', 'mdd.2019', 'ocd.2018', 'ptsd.2019', 'scz.clozuk2018', 'ts.2019')
subsets <- c('mbsv.0.2', paste('mbsv.0.2', c('eqtl', 'brain.eqtl', 'psych.eqtl'), sep = "."))

gene_sets <- list()
gene_sets_merged <- list()
for (g in gwas) {
  setwd(g)
  gene_sets_merged[[g]] <- data.frame(VARIABLE = character(), TYPE = character(), NGENES = numeric(), BETA = numeric(), BETA_STD = numeric(), SE = numeric(), P = numeric(), SUBSET = character())
  for (s in subsets) {
    n = paste(g, s, sep = "_")
    file <- paste(g, 'magma.all', s, 'sets.go.gsa.out', sep = '.')
    gene_sets[[n]] <- read.delim(file, header = TRUE, comment.char = "#", sep = "", stringsAsFactors = FALSE, row.names = "VARIABLE")
    gene_sets[[n]]$SUBSET <- n
    gene_sets_merged[[g]] <- rbind(gene_sets_merged[[g]], gene_sets[[n]])
  }
  gene_sets_merged[[g]]$BONF.P <- p.adjust(gene_sets_merged[[g]]$P, method = "bonferroni")
  gene_sets_merged[[g]]$FDR <- p.adjust(gene_sets_merged[[g]]$P, method = "BH")
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

save.image("magma.all.analyses.RData")

cols <- c("FULL_NAME", "SUBSET", "NGENES", "BETA", "BETA_STD", "SE", "P", "BONF.P", 'FDR')

# for (g in names(gene_sets_merged_5)) {
#   name = paste(g, "magma.all.gene_set_analysis.txt", sep = ".")
#   write.table(gene_sets_merged_5[[g]][cols], name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
# }

for (g in names(gene_sets_merged)) {
  name = paste(g, "magma.all.gene_set_analysis.txt", sep = ".")
  write.table(gene_sets_merged[[g]][cols], name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# Gene-level

gwas_genes <- list()
gwas_genes_merged <- list()
for (g in gwas) {
  setwd(g)
  gwas_genes_merged[[g]] <- data.frame(GENE = character(), CHR = character(), START = numeric(), STOP = numeric(), NSNPS = numeric(), NPARAM = numeric(), N = numeric(), ZSTAT = numeric(), P = numeric(), SUBSET = character())
  for (s in subsets) {
    n = paste(g, s, sep = "_")
    file <- paste(g, 'magma.all', s, 'genes.genes.out', sep = '.')
    gwas_genes[[n]] <- read.delim(file, header = TRUE, comment.char = "#", sep = "", stringsAsFactors = FALSE)
    rownames(gwas_genes[[n]]) <- gwas_genes[[n]]$GENE
    gwas_genes[[n]]$SUBSET <- n
    gwas_genes_merged[[g]] <- rbind(gwas_genes_merged[[g]], gwas_genes[[n]])
  }
  gwas_genes_merged[[g]]$BONF.P <- p.adjust(gwas_genes_merged[[g]]$P, method = "bonferroni")
  gwas_genes_merged[[g]]$FDR <- p.adjust(gwas_genes_merged[[g]]$P, method = "BH")
  setwd("../")
}

save.image("magma.all.gene.analyses.RData")

cols <- c("GENE", "SUBSET", "CHR", "START", "STOP", "NSNPS", "NPARAM", "N", "ZSTAT", "P", "BONF.P")

for (g in names(gwas_genes_merged)) {
  name = paste(g, "magma.all.gene_analysis.txt", sep = ".")
  write.table(gwas_genes_merged[[g]][cols], name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
