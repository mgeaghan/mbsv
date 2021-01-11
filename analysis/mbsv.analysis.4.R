# Import all GWAS data, prune, and export tables of variant counts
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

pruned_vars <- scan("gwas.mbsv.0.2.nonmhc.pruned.dbsnp.txt", character())

num.vars <- data.frame(num.gwas = numeric(), num.pruned.gwas = numeric(),
                       num.mbsv = numeric(), num.pruned.mbsv = numeric(), num.all.mbsv = numeric(), num.all.pruned.mbsv = numeric(),
                       num.gtex.eqtl.mbsv = numeric(), num.pruned.gtex.eqtl.mbsv = numeric(), num.all.gtex.eqtl.mbsv = numeric(), num.all.pruned.gtex.eqtl.mbsv = numeric(),
                       num.gtex.brain.eqtl.mbsv = numeric(), num.pruned.gtex.brain.eqtl.mbsv = numeric(), num.all.gtex.brain.eqtl.mbsv = numeric(), num.all.pruned.gtex.brain.eqtl.mbsv = numeric(),
                       num.psych.eqtl.mbsv = numeric(), num.pruned.psych.eqtl.mbsv = numeric(), num.all.psych.eqtl.mbsv = numeric(), num.all.pruned.psych.eqtl.mbsv = numeric())
unique.gwas <- c()
unique.mbsv <- c()
unique.gtex.eqtl.mbsv <- c()
unique.gtex.brain.eqtl.mbsv <- c()
unique.psych.eqtl.mbsv <- c()
unique.pruned.gwas <- c()
unique.pruned.mbsv <- c()
unique.pruned.gtex.eqtl.mbsv <- c()
unique.pruned.gtex.brain.eqtl.mbsv <- c()
unique.pruned.psych.eqtl.mbsv <- c()
unique.all.mbsv <- c()
unique.all.gtex.eqtl.mbsv <- c()
unique.all.gtex.brain.eqtl.mbsv <- c()
unique.all.psych.eqtl.mbsv <- c()
unique.all.pruned.mbsv <- c()
unique.all.pruned.gtex.eqtl.mbsv <- c()
unique.all.pruned.gtex.brain.eqtl.mbsv <- c()
unique.all.pruned.psych.eqtl.mbsv <- c()

for (g in names(gwas_files)) {
  load(gwas_files[[g]])
  gwas.df <- gwas.df[!is.na(gwas.df$non.mhc) &
                       !is.na(gwas.df$variable) &
                       gwas.df$non.mhc &
                       gwas.df$variable != "mbsv.any.eqtl" &
                       (gwas.df$variable == "GWAS" | !is.na(gwas.df$diffScore)),]
  g.unique.gwas <- unique(gwas.df$SNP[gwas.df$variable == "GWAS"])
  g.unique.mbsv <- unique(gwas.df$SNP[gwas.df$variable == "mbsv" & abs(gwas.df$diffScore) >= 0.2])
  g.unique.gtex.eqtl.mbsv <- unique(gwas.df$SNP[gwas.df$variable == "mbsv.gtex.eqtl" & abs(gwas.df$diffScore) >= 0.2])
  g.unique.gtex.brain.eqtl.mbsv <- unique(gwas.df$SNP[gwas.df$variable == "mbsv.gtex.brain.eqtl" & abs(gwas.df$diffScore) >= 0.2])
  g.unique.psych.eqtl.mbsv <- unique(gwas.df$SNP[gwas.df$variable == "mbsv.psychencode.eqtl" & abs(gwas.df$diffScore) >= 0.2])
  g.unique.all.mbsv <- unique(gwas.df$SNP[gwas.df$variable == "mbsv"])
  g.unique.all.gtex.eqtl.mbsv <- unique(gwas.df$SNP[gwas.df$variable == "mbsv.gtex.eqtl"])
  g.unique.all.gtex.brain.eqtl.mbsv <- unique(gwas.df$SNP[gwas.df$variable == "mbsv.gtex.brain.eqtl"])
  g.unique.all.psych.eqtl.mbsv <- unique(gwas.df$SNP[gwas.df$variable == "mbsv.psychencode.eqtl"])
  num.gwas <- length(g.unique.gwas)
  num.mbsv <- length(g.unique.mbsv)
  num.gtex.eqtl.mbsv <- length(g.unique.gtex.eqtl.mbsv)
  num.gtex.brain.eqtl.mbsv <- length(g.unique.gtex.brain.eqtl.mbsv)
  num.psych.eqtl.mbsv <- length(g.unique.psych.eqtl.mbsv)
  num.all.mbsv <- length(g.unique.all.mbsv)
  num.all.gtex.eqtl.mbsv <- length(g.unique.all.gtex.eqtl.mbsv)
  num.all.gtex.brain.eqtl.mbsv <- length(g.unique.all.gtex.brain.eqtl.mbsv)
  num.all.psych.eqtl.mbsv <- length(g.unique.all.psych.eqtl.mbsv)
  # Prune
  gwas <- gwas[gwas$SNP %in% pruned_vars,]
  gwas.df <- gwas.df[gwas.df$SNP %in% pruned_vars,]
  g.unique.pruned.gwas <- unique(gwas.df$SNP[gwas.df$variable == "GWAS"])
  g.unique.pruned.mbsv <- unique(gwas.df$SNP[gwas.df$variable == "mbsv" & abs(gwas.df$diffScore) >= 0.2])
  g.unique.pruned.gtex.eqtl.mbsv <- unique(gwas.df$SNP[gwas.df$variable == "mbsv.gtex.eqtl" & abs(gwas.df$diffScore) >= 0.2])
  g.unique.pruned.gtex.brain.eqtl.mbsv <- unique(gwas.df$SNP[gwas.df$variable == "mbsv.gtex.brain.eqtl" & abs(gwas.df$diffScore) >= 0.2])
  g.unique.pruned.psych.eqtl.mbsv <- unique(gwas.df$SNP[gwas.df$variable == "mbsv.psychencode.eqtl" & abs(gwas.df$diffScore) >= 0.2])
  g.unique.all.pruned.mbsv <- unique(gwas.df$SNP[gwas.df$variable == "mbsv"])
  g.unique.all.pruned.gtex.eqtl.mbsv <- unique(gwas.df$SNP[gwas.df$variable == "mbsv.gtex.eqtl"])
  g.unique.all.pruned.gtex.brain.eqtl.mbsv <- unique(gwas.df$SNP[gwas.df$variable == "mbsv.gtex.brain.eqtl"])
  g.unique.all.pruned.psych.eqtl.mbsv <- unique(gwas.df$SNP[gwas.df$variable == "mbsv.psychencode.eqtl"])
  num.pruned.gwas <- length(g.unique.pruned.gwas)
  num.pruned.mbsv <- length(g.unique.pruned.mbsv)
  num.pruned.gtex.eqtl.mbsv <- length(g.unique.pruned.gtex.eqtl.mbsv)
  num.pruned.gtex.brain.eqtl.mbsv <- length(g.unique.pruned.gtex.brain.eqtl.mbsv)
  num.pruned.psych.eqtl.mbsv <- length(g.unique.pruned.psych.eqtl.mbsv)
  num.all.pruned.mbsv <- length(g.unique.all.pruned.mbsv)
  num.all.pruned.gtex.eqtl.mbsv <- length(g.unique.all.pruned.gtex.eqtl.mbsv)
  num.all.pruned.gtex.brain.eqtl.mbsv <- length(g.unique.all.pruned.gtex.brain.eqtl.mbsv)
  num.all.pruned.psych.eqtl.mbsv <- length(g.unique.all.pruned.psych.eqtl.mbsv)
  # Tabulate
  num.vars.rows <- rownames(num.vars)
  num.vars <- rbind(num.vars, data.frame(num.gwas = num.gwas, num.pruned.gwas = num.pruned.gwas,
                                         num.mbsv = num.mbsv, num.pruned.mbsv = num.pruned.mbsv, num.all.mbsv = num.all.mbsv, num.all.pruned.mbsv = num.all.pruned.mbsv,
                                         num.gtex.eqtl.mbsv = num.gtex.eqtl.mbsv, num.pruned.gtex.eqtl.mbsv = num.pruned.gtex.eqtl.mbsv, num.all.gtex.eqtl.mbsv = num.all.gtex.eqtl.mbsv, num.all.pruned.gtex.eqtl.mbsv = num.all.pruned.gtex.eqtl.mbsv,
                                         num.gtex.brain.eqtl.mbsv = num.gtex.brain.eqtl.mbsv, num.pruned.gtex.brain.eqtl.mbsv = num.pruned.gtex.brain.eqtl.mbsv, num.all.gtex.brain.eqtl.mbsv = num.all.gtex.brain.eqtl.mbsv, num.all.pruned.gtex.brain.eqtl.mbsv = num.all.pruned.gtex.brain.eqtl.mbsv,
                                         num.psych.eqtl.mbsv = num.psych.eqtl.mbsv, num.pruned.psych.eqtl.mbsv = num.pruned.psych.eqtl.mbsv, num.all.psych.eqtl.mbsv = num.all.psych.eqtl.mbsv, num.all.pruned.psych.eqtl.mbsv = num.all.pruned.psych.eqtl.mbsv))
  rownames(num.vars) <- c(num.vars.rows, g)
  unique.gwas <- unique(c(unique.gwas, g.unique.gwas))
  unique.mbsv <- unique(c(unique.mbsv, g.unique.mbsv))
  unique.gtex.eqtl.mbsv <- unique(c(unique.gtex.eqtl.mbsv, g.unique.gtex.eqtl.mbsv))
  unique.gtex.brain.eqtl.mbsv <- unique(c(unique.gtex.brain.eqtl.mbsv, g.unique.gtex.brain.eqtl.mbsv))
  unique.psych.eqtl.mbsv <- unique(c(unique.psych.eqtl.mbsv, g.unique.psych.eqtl.mbsv))
  unique.pruned.gwas <- unique(c(unique.pruned.gwas, g.unique.pruned.gwas))
  unique.pruned.mbsv <- unique(c(unique.pruned.mbsv, g.unique.pruned.mbsv))
  unique.pruned.gtex.eqtl.mbsv <- unique(c(unique.pruned.gtex.eqtl.mbsv, g.unique.pruned.gtex.eqtl.mbsv))
  unique.pruned.gtex.brain.eqtl.mbsv <- unique(c(unique.pruned.gtex.brain.eqtl.mbsv, g.unique.pruned.gtex.brain.eqtl.mbsv))
  unique.pruned.psych.eqtl.mbsv <- unique(c(unique.pruned.psych.eqtl.mbsv, g.unique.pruned.psych.eqtl.mbsv))
  unique.all.mbsv <- unique(c(unique.all.mbsv, g.unique.all.mbsv))
  unique.all.gtex.eqtl.mbsv <- unique(c(unique.all.gtex.eqtl.mbsv, g.unique.all.gtex.eqtl.mbsv))
  unique.all.gtex.brain.eqtl.mbsv <- unique(c(unique.all.gtex.brain.eqtl.mbsv, g.unique.all.gtex.brain.eqtl.mbsv))
  unique.all.psych.eqtl.mbsv <- unique(c(unique.all.psych.eqtl.mbsv, g.unique.all.psych.eqtl.mbsv))
  unique.all.pruned.mbsv <- unique(c(unique.all.pruned.mbsv, g.unique.all.pruned.mbsv))
  unique.all.pruned.gtex.eqtl.mbsv <- unique(c(unique.all.pruned.gtex.eqtl.mbsv, g.unique.all.pruned.gtex.eqtl.mbsv))
  unique.all.pruned.gtex.brain.eqtl.mbsv <- unique(c(unique.all.pruned.gtex.brain.eqtl.mbsv, g.unique.all.pruned.gtex.brain.eqtl.mbsv))
  unique.all.pruned.psych.eqtl.mbsv <- unique(c(unique.all.pruned.psych.eqtl.mbsv, g.unique.all.pruned.psych.eqtl.mbsv))
  # save(gwas, gwas.df, gwas_files_char, gwas_es_char, file = paste(g, ".gwas.nonmhc.pruned.RData", sep = ""))
}
write.table(num.vars, "num.vars.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
num.all.vars <- data.frame(subset = c("gwas", "mbsv", "gtex.eqtl.mbsv", "gtex.brain.eqtl.mbsv", "psychencode.eqtl.mbsv"),
                           num.vars = c(length(unique.gwas), length(unique.mbsv), length(unique.gtex.eqtl.mbsv), length(unique.gtex.brain.eqtl.mbsv), length(unique.psych.eqtl.mbsv)),
                           num.pruned.vars = c(length(unique.pruned.gwas), length(unique.pruned.mbsv), length(unique.pruned.gtex.eqtl.mbsv), length(unique.pruned.gtex.brain.eqtl.mbsv), length(unique.pruned.psych.eqtl.mbsv)),
                           num.all.vars = c(NA, length(unique.all.mbsv), length(unique.all.gtex.eqtl.mbsv), length(unique.all.gtex.brain.eqtl.mbsv), length(unique.all.psych.eqtl.mbsv)),
                           num.all.pruned.vars = c(NA, length(unique.all.pruned.mbsv), length(unique.all.pruned.gtex.eqtl.mbsv), length(unique.all.pruned.gtex.brain.eqtl.mbsv), length(unique.all.pruned.psych.eqtl.mbsv)))
write.table(num.all.vars, "num.all.vars.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
