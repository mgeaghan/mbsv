# Export variants for MAGMA analyses
# Author: Michael Geaghan (2020)

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

for(g in names(gwas_data)) {
  load(gwas_data[[g]])
  gwas.df <- gwas.df[!is.na(gwas.df$non.mhc) &
                       !is.na(gwas.df$variable) &
                       gwas.df$non.mhc &
                       gwas.df$variable != "mbsv.any.eqtl" &
                       (gwas.df$variable == "GWAS" | (!is.na(gwas.df$diffScore) & abs(gwas.df$diffScore) >= 0.2)),]
  mbsv.snp <- unique(gwas.df$SNP[gwas.df$variable == "mbsv"])
  gwas.df.nonmbsv <- gwas.df[gwas.df$variable == "GWAS" & !(gwas.df$SNP %in% mbsv.snp),]
  gwas.df.nonmbsv$variable <- "GWAS.non.mbsv"
  gwas.df <- rbind(gwas.df, gwas.df.nonmbsv)
  # gwas.df <- gwas.df[gwas.df$variable != "GWAS",]
  # gwas.df$variable <- as.character(gwas.df$variable)
  for (p in c("all", "0.5", "0.05", "0.005")) {
    if (p == "all") {
      tmp.df <- gwas.df
    } else {
      tmp.df <- gwas.df[!is.na(gwas.df$P) & gwas.df$P < as.numeric(p),]
    }
    # all mbsvs
    vars.mbsv.0.2 <- unique(tmp.df$SNP[tmp.df$variable == "mbsv"])
    write.table(vars.mbsv.0.2, paste("./magma/", g, ".", p, ".mbsv.0.2.nonmhc.dbsnp.txt", sep = ""), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE, eol = "\n")
    # mbsv eqtls
    vars.mbsv.0.2.eqtl <- unique(tmp.df$SNP[tmp.df$variable == "mbsv.gtex.eqtl"])
    write.table(vars.mbsv.0.2.eqtl, paste("./magma/", g, ".", p, ".mbsv.0.2.eqtl.nonmhc.dbsnp.txt", sep = ""), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE, eol = "\n")
    # mbsv brain eqtls
    vars.mbsv.0.2.brain.eqtl <- unique(tmp.df$SNP[tmp.df$variable == "mbsv.gtex.brain.eqtl"])
    write.table(vars.mbsv.0.2.brain.eqtl, paste("./magma/", g, ".", p, ".mbsv.0.2.brain.eqtl.nonmhc.dbsnp.txt", sep = ""), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE, eol = "\n")
    # mbsv psychencode eqtls
    vars.mbsv.0.2.psych.eqtl <- unique(tmp.df$SNP[tmp.df$variable == "mbsv.psychencode.eqtl"])
    write.table(vars.mbsv.0.2.psych.eqtl, paste("./magma/", g, ".", p, ".mbsv.0.2.psych.eqtl.nonmhc.dbsnp.txt", sep = ""), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE, eol = "\n")
    # non-mbsvs
    vars.mbsv.0.2.nonmbsv <- unique(tmp.df$SNP[tmp.df$variable == "GWAS.non.mbsv"])
    write.table(vars.mbsv.0.2.nonmbsv, paste("./magma/", g, ".", p, ".nonmbsv.0.2.nonmhc.dbsnp.txt", sep = ""), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE, eol = "\n")
    # gwas
    vars.mbsv.0.2.gwas <- unique(tmp.df$SNP[tmp.df$variable == "GWAS"])
    write.table(vars.mbsv.0.2.gwas, paste("./magma/", g, ".", p, ".gwas.nonmhc.dbsnp.txt", sep = ""), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE, eol = "\n")
  }
}
