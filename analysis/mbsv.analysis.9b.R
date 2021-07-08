load('acat.RData')
num_comparisons <- sum(sapply(acat_df, dim)[1,])
acat_df_corr <- acat_df
for (x in names(acat_df_corr)) {
  acat_df_corr[[x]]$gwas <- x
  acat_df_corr[[x]]$mirna <- rownames(acat_df_corr[[x]])
}
acat_df_corr <- do.call(rbind, acat_df_corr)
acat_df_corr[1:4] <- p.adjust(as.matrix(acat_df_corr[1:4]), method = "bonferroni")

acat_all_df_corr <- acat_all_df
acat_all_df_corr[1:4] <- p.adjust(as.matrix(acat_all_df_corr[1:4]), method = "bonferroni")

write.table(acat_df_corr, "acat.gwas.corr.txt", sep = "\t", col.names = TRUE, row.names = TRUE)
write.table(acat_all_df_corr, "acat.all.corr.txt", sep = "\t", col.names = TRUE, row.names = TRUE)
