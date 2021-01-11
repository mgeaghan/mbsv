# Run ACAT on miRNA families
# Author: Michael Geaghan(2020)

source("acat.R")

load("dbmts.RData")

# Read in GWAS data
gwas_data <- list(adhd.2018 = './adhd.2018.gwas.RData',
                  an.2019 = './an.2019.gwas.RData',
                  asd.2019 = './asd.2019.gwas.RData',
                  bip.2018 = './bip.2018.gwas.RData',
                  mdd.2019 = './mdd.2019.gwas.RData',
                  ocd.2018 = './ocd.2018.gwas.RData',
                  ptsd.2019 = './ptsd.2019.gwas.RData',
                  scz.clozuk2018 = './scz.clozuk2018.gwas.RData',
                  ts.2019 = './ts.2019.gwas.RData')

mirna.families <- unique(dbmts.agree.best.best.0.2.short$mirna.family)

acat_df <- list()
for(g in names(gwas_data)) {
  acat_df[[g]] <- data.frame(fam = mirna.families, mbsv = NA, mbsv.gtex.eqtl = NA, mbsv.gtex.brain.eqtl = NA, mbsv.psychencode.eqtl = NA, row.names = 1)
}
acat_all_df <- data.frame(gwas = names(gwas_data), mbsv = NA, mbsv.gtex.eqtl = NA, mbsv.gtex.brain.eqtl = NA, mbsv.psychencode.eqtl = NA, row.names = 1)

minP <- .Machine$double.xmin
maxP <- 1 - .Machine$double.eps

for(g in names(gwas_data)) {
  load(gwas_data[[g]])
  gwas.df <- gwas.df[!is.na(gwas.df$non.mhc) &
                       !is.na(gwas.df$variable) &
                       gwas.df$non.mhc &
                       gwas.df$variable != "mbsv.any.eqtl" &
                       (gwas.df$variable == "GWAS" | (!is.na(gwas.df$diffScore) & abs(gwas.df$diffScore) >= 0.2)),]
  filters <- c("mbsv", "mbsv.gtex.eqtl", "mbsv.gtex.brain.eqtl", "mbsv.psychencode.eqtl")
  for(f in filters) {
    tmp_g <- gwas.df[gwas.df$variable == f,]
    tmp_d <- dbmts.agree.best.best.0.2.short
    # All MBSVs (weighted by |diffScore|)
    p = tmp_g$P[tmp_g$SNP %in% tmp_d$dbsnp]
    p[p < minP] <- minP
    p[p > maxP] <- maxP
    w = abs(tmp_g$diffScore[tmp_g$SNP %in% tmp_d$dbsnp])
    acat_all_df[g,f] = ACAT(p, w)
    # Individual families (weighted by |diffScore|)
    acat_df[[g]][f] <- do.call(c, lapply(mirna.families, function(x) {
      tmp_d_m <- tmp_d[tmp_d$mirna.family == x,]
      p = tmp_g$P[tmp_g$SNP %in% tmp_d_m$dbsnp]
      p[p < minP] <- minP
      p[p > maxP] <- maxP
      w = abs(tmp_g$diffScore[tmp_g$SNP %in% tmp_d_m$dbsnp])
      return(ACAT(p, w))
    }))
  }
}

save(acat_all_df, acat_df, mirna.families, file = "acat.RData")
write.table(acat_all_df, "acat.all.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
for (g in names(acat_df)) {
  write.table(acat_df[[g]], paste("acat.", g, ".txt", sep = ""), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
}
