setwd('magma/meta')
load("magma.meta.all.analyses.RData")
gene_sets_merged_full <- do.call(rbind, gene_sets_merged)
gene_sets_merged_full$BONF.P <- p.adjust(gene_sets_merged_full$P, method = "bonferroni")
gene_sets_merged_full$FDR <- p.adjust(gene_sets_merged_full$P, method = "BH")
write.table(gene_sets_merged_full, "magma.meta.all.gene_set_analysis.full.txt", quote = FALSE, sep = "\t", row.names = FALSE)

gene_sets_merged_5_full <- do.call(rbind, gene_sets_merged_5)
gene_sets_merged_5_full$BONF.P <- p.adjust(gene_sets_merged_5_full$P, method = "bonferroni")
gene_sets_merged_5_full$FDR <- p.adjust(gene_sets_merged_5_full$P, method = "BH")
write.table(gene_sets_merged_5_full, "magma.meta.all.gene_set_5_analysis.full.txt", quote = FALSE, sep = "\t", row.names = FALSE)

rm(gene_sets, gene_sets_5, gene_sets_merged, gene_sets_merged_5, gene_sets_merged_full, gene_sets_merged_5_full)

load("magma.meta.all.gene.analyses.RData")
meta_genes_merged_full <- do.call(rbind, meta_genes_merged)
meta_genes_merged_full$BONF.P <- p.adjust(meta_genes_merged_full$P, method = "bonferroni")
meta_genes_merged_full$FDR <- p.adjust(meta_genes_merged_full$P, method = "BH")
write.table(meta_genes_merged_full, "magma.meta.all.gene_analysis.full.txt", quote = FALSE, sep = "\t", row.names = FALSE)
