# Run ACAT meta-analyses using Stouffer's weighted Z method
# Author: Michael Geaghan (2020)

load("acat.RData")
load("dbmts.RData")

# Stouffer's method
stouffer <- function(p, w = rep(1, length(p))) {
  p_new <- p
  minP <- .Machine$double.xmin
  maxP <- 1 - .Machine$double.eps
  p_new[p_new < minP] <- minP
  p_new[p_new > maxP] <- maxP
  Zi <- qnorm(1 - p_new)
  Z <- sum(w * Zi)/sqrt(sum(w^2))
  p_new_2 <- 1 - pnorm(Z)
  return(p_new_2)
}

gwas_list <- c('adhd.2018', 'an.2019', 'asd.2019', 'bip.2018',
               'mdd.2019', 'ocd.2018', 'ptsd.2019',
               'scz.clozuk2018', 'ts.2019')

# meta-analysis lists
meta_list <- list(meta.scz.bip = c("scz.clozuk2018", "bip.2018"),
                  meta.scz.mdd = c("scz.clozuk2018", "mdd.2019"),
                  meta.scz.bip.mdd = c("scz.clozuk2018", "bip.2018", "mdd.2019"),
                  meta.all = gwas_list)

filters <- c("mbsv", "mbsv.gtex.eqtl", "mbsv.gtex.brain.eqtl", "mbsv.psychencode.eqtl")

sample.sizes <- list()
for(g in gwas_list) {
  tmp_file <- paste("../gwas/", gsub("\\.", "\\-", g), ".n.txt", sep = "")
  tmp_df <- read.delim(tmp_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  sample.sizes[[g]] <- tmp_df[1, 2]
}

mirna.families <- unique(dbmts.agree.best.best.0.2.short$mirna.family)

acat_meta_df <- list()
for(m in names(meta_list)) {
  acat_meta_df[[m]] <- data.frame(fam = mirna.families, mbsv = NA, mbsv.gtex.eqtl = NA, mbsv.gtex.brain.eqtl = NA, mbsv.psychencode.eqtl = NA, row.names = 1)
}
acat_all_meta_df <- data.frame(meta = names(meta_list), mbsv = NA, mbsv.gtex.eqtl = NA, mbsv.gtex.brain.eqtl = NA, mbsv.psychencode.eqtl = NA, row.names = 1)

# Run meta analyses
for(m in names(meta_list)) {
  for(f in filters) {
    n = sapply(meta_list[[m]], function(x) {
      return(sample.sizes[[x]])
    })
    # All MBSVs
    p <- acat_all_df[meta_list[[m]], f]
    acat_all_meta_df[m, f] <- stouffer(as.numeric(p), as.numeric(n))
    # Individual families
    tmp_df <- data.frame(fam = mirna.families, row.names = 1)
    for(g in meta_list[[m]]) {
      tmp_df[[g]] <- acat_df[[g]][mirna.families, f]
    }
    acat_meta_df[[m]][f] <- sapply(mirna.families, function(x) {
      p <- tmp_df[x,]
      return(stouffer(as.numeric(p), as.numeric(n)))
    })
  }
}

save(acat_all_meta_df, acat_meta_df, mirna.families, file = "acat.meta.RData")
write.table(acat_all_meta_df, "acat.all.meta.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
for (m in names(acat_meta_df)) {
  write.table(acat_meta_df[[m]], paste("acat.", m, ".txt", sep = ""), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
}
