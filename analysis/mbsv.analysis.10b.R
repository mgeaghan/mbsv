load("acat.meta.RData")
num_comparisons <- sum(sapply(acat_meta_df, dim)[1,])
acat_meta_df_corr <- acat_meta_df
for (x in names(acat_meta_df_corr)) {
  acat_meta_df_corr[[x]]$meta <- x
  acat_meta_df_corr[[x]]$mirna <- rownames(acat_meta_df_corr[[x]])
}
acat_meta_df_corr <- do.call(rbind, acat_meta_df_corr)
acat_meta_df_corr[1:4] <- p.adjust(as.matrix(acat_meta_df_corr[1:4]), method = "bonferroni")

acat_all_meta_df_corr <- acat_all_meta_df
acat_all_meta_df_corr[1:4] <- p.adjust(as.matrix(acat_all_meta_df_corr[1:4]), method = "bonferroni")

write.table(acat_meta_df_corr, "acat.meta.gwas.corr.txt", sep = "\t", col.names = TRUE, row.names = TRUE)
write.table(acat_all_meta_df_corr, "acat.all.meta.corr.txt", sep = "\t", col.names = TRUE, row.names = TRUE)
