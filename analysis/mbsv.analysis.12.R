# Investigate relationship between effect sizes of MBSVs and diffScore in each disorder, subset by significantly associated miRNA families
# Author: Michael Geaghan (2020)

load("dbmts.RData")
load("acat.RData")
acat_results <- read.delim("acat.gwas.corr.txt")
library(ggplot2)
library(reshape2)

result.list <- list(mir.df.list = list(),
                    linmod.score.list = list(), linmod.score.summ.list = list(), linmod.score.df.list = list(),
                    linmod.score_sign.list = list(), linmod.score_sign.summ.list = list(), linmod.score_sign.df.list = list(),
                    plot.score.list = list(), plot.score_sign.list = list())

for (g in names(acat_df)) {
  print(paste("=====", g, "=====", sep = " "))
  load(paste(g, ".gwas.nonmhc.pruned.RData", sep = ""))
  gwas_acat_results <- acat_results[acat_results$gwas == g,]
  sig_mirs <- gwas_acat_results$mirna[
    gwas_acat_results$mbsv < 0.05 |
      gwas_acat_results$mbsv.gtex.eqtl < 0.05 |
      gwas_acat_results$mbsv.gtex.brain.eqtl < 0.05 |
      gwas_acat_results$mbsv.psychencode.eqtl < 0.05
    ]
  gwas_summ <- read.delim(paste("../gwas/", gsub(".", "-", g, fixed = TRUE), ".gwas", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  gwas_summ <- gwas_summ[!(gwas_summ$SNP %in% unique(gwas_summ$SNP[duplicated(gwas_summ$SNP)])),]
  rownames(gwas_summ) <- gwas_summ$SNP
  # sig_mirs <- rownames(acat_df[[g]][acat_df[[g]]$mbsv < 2.5e-6 | acat_df[[g]]$mbsv.gtex.eqtl < 2.5e-6 | acat_df[[g]]$mbsv.gtex.brain.eqtl < 2.5e-6 | acat_df[[g]]$mbsv.psychencode.eqtl < 2.5e-6,])
  mir.df <- list()
  linmod.list <- list()
  linmod.summ.list <- list()
  linmod.df <- data.frame(mir.family = character(), beta = numeric(), se = numeric(), t = numeric(), p = numeric())
  linmod.2.list <- list()
  linmod.2.summ.list <- list()
  linmod.2.df <- data.frame(mir.family = character(), beta = numeric(), se = numeric(), t = numeric(), p = numeric())
  plot.list <- list()
  plot.2.list <- list()
  for (m in sig_mirs) {
    mir.dbsnp <- dbmts.agree.best.best.0.2.short[dbmts.agree.best.best.0.2.short$mirna.family == m,]$dbsnp
    mir.dbsnp <- mir.dbsnp[!is.na(mir.dbsnp)]
    mir.df[[m]] <- gwas.df[gwas.df$SNP %in% mir.dbsnp & gwas.df$variable == "GWAS",]
    if (g == "an.2019") {
      mir.df[[m]]$a1 <- gwas_summ[mir.df[[m]]$SNP,]$ALT
    } else {
      mir.df[[m]]$a1 <- gwas_summ[mir.df[[m]]$SNP,]$A1
    }
    mir.df[[m]]$ref <- as.data.frame(do.call(rbind, strsplit(mir.df[[m]]$variant, split = ":", fixed = TRUE)))$V3
    mir.df[[m]]$alt <- as.data.frame(do.call(rbind, strsplit(mir.df[[m]]$variant, split = ":", fixed = TRUE)))$V4
    mir.df[[m]]$alt_es <- apply(mir.df[[m]], 1, function(x) {
      ref = as.character(x[10])
      alt = as.character(x[11])
      a1 = as.character(x[9])
      es = as.numeric(x[2])
      if (is.na(a1) || is.na(alt) || is.na(ref)) {
        return(NA)
      } else if (a1 == alt) {
        return(es)
      } else if (a1 == ref) {
        return(-es)
      } else {
        return(NA)
      }
    })
    mir.df[[m]]$diffScore.sign <- sign(mir.df[[m]]$diffScore)
    mir.df[[m]]$negLog10p <- -log10(mir.df[[m]]$P)
    mir.df[[m]] <- mir.df[[m]][!grepl("((A\\:T)|(C\\:G)|(T\\:A)|(G\\:C))$", mir.df[[m]]$variant, perl = TRUE),]
    linmod.list[[m]] <- lm(alt_es ~ diffScore.sign, data = mir.df[[m]])
    linmod.summ.list[[m]] <- summary(linmod.list[[m]])
    linmod.df <- rbind(linmod.df, data.frame(mir.family = m,
                                             beta = linmod.summ.list[[m]]$coefficients["diffScore.sign", "Estimate"],
                                             se = linmod.summ.list[[m]]$coefficients["diffScore.sign", "Std. Error"],
                                             t = linmod.summ.list[[m]]$coefficients["diffScore.sign", "t value"],
                                             p = linmod.summ.list[[m]]$coefficients["diffScore.sign", "Pr(>|t|)"]))
    linmod.2.list[[m]] <- lm(alt_es ~ diffScore, data = mir.df[[m]])
    linmod.2.summ.list[[m]] <- summary(linmod.2.list[[m]])
    linmod.2.df <- rbind(linmod.2.df, data.frame(mir.family = m,
                                                  beta = linmod.2.summ.list[[m]]$coefficients["diffScore", "Estimate"],
                                                  se = linmod.2.summ.list[[m]]$coefficients["diffScore", "Std. Error"],
                                                  t = linmod.2.summ.list[[m]]$coefficients["diffScore", "t value"],
                                                  p = linmod.2.summ.list[[m]]$coefficients["diffScore", "Pr(>|t|)"]))
    plot.list[[m]] <- ggplot(data = mir.df[[m]], aes(x = diffScore.sign, y = alt_es, col = negLog10p)) +
      geom_point() + geom_smooth(method = "lm") +
      geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
    plot.2.list[[m]] <- ggplot(data = mir.df[[m]], aes(x = diffScore, y = alt_es, col = negLog10p)) +
      geom_point() + geom_smooth(method = "lm") +
      geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
  }
  linmod.df$fdr <- p.adjust(linmod.df$p, method = "BH")
  linmod.df$bonf.p <- p.adjust(linmod.df$p, method = "bonferroni")
  linmod.2.df$fdr <- p.adjust(linmod.2.df$p, method = "BH")
  linmod.2.df$bonf.p <- p.adjust(linmod.2.df$p, method = "bonferroni")
  result.list$mir.df.list[[g]] <- mir.df
  result.list$linmod.score.list[[g]] <- linmod.2.list
  result.list$linmod.score.summ.list[[g]] <- linmod.2.summ.list
  result.list$linmod.score.df.list[[g]] <- linmod.2.df
  result.list$linmod.score_sign.list[[g]] <- linmod.list
  result.list$linmod.score_sign.summ.list[[g]] <- linmod.summ.list
  result.list$linmod.score_sign.df.list[[g]] <- linmod.df
  result.list$plot.score.list[[g]] <- plot.2.list
  result.list$plot.score_sign.list[[g]] <- plot.list
}
save(result.list, file = "mirna.linmods.new.RData")

df <- do.call(rbind, result.list$linmod.score.df.list)
df$fdr <- p.adjust(df$p, method = "BH")
df$bonf.p <- p.adjust(df$p, method = "bonferroni")
df$gwas <- gsub("\\.\\d+$", "", rownames(df))
write.table(df, "mirna.linmod.new.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

mir.plot.list <- apply(df[df$p < 0.05,], 1, function(x) {
  y <- as.character(x[[8]])
  z <- as.character(x[[1]])
  return(
    result.list$plot.score.list[[y]][[z]] +
      theme(
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)
      ) +
      ggtitle(paste(z, " - ", y, sep = ""))
  )
})

library(ggpubr)
ggsave("mirna.linmods.new.scorevseffectsize.png", ggarrange(plotlist = mir.plot.list), device = "png", width = 8, height = 8, dpi = 600)


df <- do.call(rbind, result.list$linmod.score_sign.df.list)
df$fdr <- p.adjust(df$p, method = "BH")
df$bonf.p <- p.adjust(df$p, method = "bonferroni")
df$gwas <- gsub("\\.\\d+$", "", rownames(df))
write.table(df, "mirna.linmod_sign.new.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

mir.df.list <- apply(df[df$p < 0.05,], 1, function(x) {
  y <- as.character(x[[8]])
  z <- as.character(x[[1]])
  d <- result.list$mir.df.list[[y]][[z]]
  d$gwas <- y
  d$mirna <- z
  return(d)
})
mir.df <- do.call(rbind, mir.df.list)
mir.df$diffScore.direction <- "Decreased affinity"
mir.df$diffScore.direction[mir.df$diffScore < 0] <- "Increased affinity"
mir.df$sig <- apply(mir.df, 1, function(x) {
  e <- as.numeric(x[[12]])
  y <- as.character(x[[15]])
  z <- as.character(x[[16]])
  p <- df[df$gwas == y & df$mir.family == z, "p"]
  if(e == max(mir.df$alt_es[mir.df$gwas == y & mir.df$mirna == z])) {
    return(paste("p = ", as.character(round(p, 3)), sep = ""))
  } else {
    return("")
  }
})

g <- ggplot(data = mir.df, aes(x = mirna, y = alt_es)) +
  geom_boxplot(aes(color = diffScore.direction)) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  geom_text(aes(label = sig, group = mirna), position = position_dodge(0.9), vjust = -1, size = 6) +
  labs(x = "", y = "Alternate allele log odds ratio", color = "Predicted MBSV effect:") +
  ylim(c(NA, max(mir.df$alt_es) * 1.2)) +
  facet_wrap(~ gwas, scales = "free_x"); g

ggsave("mirna.linmods_sign.new.scorevseffectsize.png", g, device = "png", width = 10, height = 8, dpi = 600)
