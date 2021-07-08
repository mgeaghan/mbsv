ks_files <- dir(".", pattern = "ks\\..*\\.nonmbsvbrain3utr\\.txt")

for (ks in ks_files) {
  df <- read.delim(ks)
  df[c(5,9,13,17)] <- p.adjust(as.matrix(df[c(4,8,12,16)]), method = "BH")
  write.table(df, gsub(".txt", ".corr.txt", ks, fixed = TRUE), quote = FALSE, sep = "\t")
}

rq_files <- dir(".", pattern = "rq\\..*\\.nonmbsvbrain3utr\\.txt")

for (r in rq_files) {
  df <- read.delim(r)
  df$mbsv_fdr <- p.adjust(df$mbsv_p, method = "BH")
  write.table(df, gsub(".txt", ".corr.txt", r, fixed = TRUE), quote = FALSE, sep = "\t")
}
