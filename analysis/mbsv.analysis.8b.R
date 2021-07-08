load('gwas.num_sig_vars.5e-08.nonmbsvbrain3utr.RData')
p <- c(
  numVar.df$p.mbsv.sig,
  numVar.df$p.gtex.eqtl.mbsv.sig,
  numVar.df$p.gtex.brain.eqtl.mbsv.sig,
  numVar.df$p.psych.eqtl.mbsv.sig
)
q <- p.adjust(p, method = "BH")
numVar.df.full <- numVar.df
numVar.df.full$fdr.mbsv.sig <- q[1:9]
numVar.df.full$fdr.gtex.eqtl.mbsv.sig <- q[10:18]
numVar.df.full$fdr.gtex.brain.eqtl.mbsv.sig <- q[19:27]
numVar.df.full$fdr.psych.eqtl.mbsv.sig <- q[28:36]

save(numVar.df.full, file = "gwas.num_sig_vars.5e-08.nonmbsvbrain3utr.full.RData")
write.table(numVar.df.full, file = "gwas.num_sig_vars.5e-08.nonmbsvbrain3utr.full.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

load('gwas.num_sig_vars.1e-05.nonmbsvbrain3utr.RData')
p <- c(
  numVar.df$p.mbsv.sig,
  numVar.df$p.gtex.eqtl.mbsv.sig,
  numVar.df$p.gtex.brain.eqtl.mbsv.sig,
  numVar.df$p.psych.eqtl.mbsv.sig
)
q <- p.adjust(p, method = "BH")
numVar.df.full <- numVar.df
numVar.df.full$fdr.mbsv.sig <- q[1:9]
numVar.df.full$fdr.gtex.eqtl.mbsv.sig <- q[10:18]
numVar.df.full$fdr.gtex.brain.eqtl.mbsv.sig <- q[19:27]
numVar.df.full$fdr.psych.eqtl.mbsv.sig <- q[28:36]

save(numVar.df.full, file = "gwas.num_sig_vars.1e-05.nonmbsvbrain3utr.full.RData")
write.table(numVar.df.full, file = "gwas.num_sig_vars.1e-05.nonmbsvbrain3utr.full.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
