# Import all GWAS data and extract GWAS and MBSV dbSNP IDs
# Author: Michael Geaghan (2020)

# Read in GWAS data
gwas_files <- list(adhd.2018 = './adhd.2018.gwas.RData',
                   an.2019 = './an.2019.gwas.RData',
                   asd.2019 = './asd.2019.gwas.RData',
                   bip.2018 = './bip.2018.gwas.RData',
                   mdd.2019 = './mdd.2019.gwas.RData',
                   ocd.2018 = './ocd.2018.gwas.RData',
                   ptsd.2019 = './ptsd.2019.gwas.RData',
                   scz.clozuk2018 = './scz.clozuk2018.gwas.RData',
                   ts.2019 = './ts.2019.gwas.RData')

all_mbsvs <- c()
all_gwas.vars <- c()
for (g in names(gwas_files)) {
  load(gwas_files[[g]])
  mbsvs <- unique(gwas.df$SNP[!is.na(gwas.df$non.mhc) &
                                !is.na(gwas.df$variable) &
                                !is.na(gwas.df$diffScore) &
                                gwas.df$non.mhc &
                                gwas.df$variable == "mbsv" &
                                abs(gwas.df$diffScore) >= 0.2])
  gwas.vars <- unique(gwas.df$SNP[!is.na(gwas.df$non.mhc) &
                                    !is.na(gwas.df$variable) &
                                    gwas.df$non.mhc &
                                    gwas.df$variable == "GWAS"])
  all_mbsvs <- unique(c(all_mbsvs, mbsvs))
  all_gwas.vars <- unique(c(all_gwas.vars, gwas.vars))
}

write.table(all_mbsvs, "mbsv.0.2.nonmhc.dbsnp.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_gwas.vars, "gwas.nonmhc.dbsnp.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
