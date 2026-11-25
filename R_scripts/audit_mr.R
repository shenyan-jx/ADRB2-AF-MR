# -*- coding: utf-8 -*-
# 审核校准：从本地 harmonised dat 重新计算异质性与多效性，产出可引用的权威数值
suppressMessages(library(TwoSampleMR))

audit <- function(path, label) {
  dat <- readRDS(path)
  dat <- dat[dat$mr_keep, ]
  cat("=====", label, " n =", nrow(dat), "=====\n")
  het <- mr_heterogeneity(dat)
  cat("-- heterogeneity --\n")
  print(as.data.frame(het[, c("method", "Q", "Q_df", "Q_pval")]), row.names = FALSE)
  ple <- mr_pleiotropy_test(dat)
  cat("-- pleiotropy (Egger intercept) --\n")
  print(as.data.frame(ple[, c("egger_intercept", "se", "pval")]), row.names = FALSE)
  cat("-- per-SNP Wald --\n")
  w <- dat$beta.outcome / dat$beta.exposure
  se <- abs(dat$se.outcome / dat$beta.exposure)
  print(data.frame(SNP = dat$SNP, b = round(w, 4), se = round(se, 4),
                   OR = round(exp(w), 3),
                   lci = round(exp(w - 1.96 * se), 3),
                   uci = round(exp(w + 1.96 * se), 3),
                   p = signif(2 * pnorm(-abs(w / se)), 3),
                   F = round(dat$beta.exposure^2 / dat$se.exposure^2, 1)), row.names = FALSE)
  write.csv(het, paste0("audit_", label, "_heterogeneity.csv"), row.names = FALSE)
  write.csv(ple, paste0("audit_", label, "_pleiotropy.csv"), row.names = FALSE)
  invisible(NULL)
}

audit("mr_af/dat_harmonised.rds", "ADRB2")
audit("mr_adrb1/dat_harmonised.rds", "ADRB1")
cat("DONE\n")
