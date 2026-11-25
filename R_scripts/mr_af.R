suppressMessages({library(TwoSampleMR)})
dir.create("mr_af", showWarnings = FALSE)

## ---------- 1. 暴露 ----------
exp <- readRDS("exp_ADRB2.rds")
cat("STEP1 instruments n =", nrow(exp), "\n")
cat("  SNPs:", paste(exp$SNP, collapse = ", "), "\n")
if (nrow(exp) < 3) { cat("STOP: insufficient instruments\n"); quit(save = "no", status = 2) }
exp$F <- exp$beta.exposure^2 / exp$se.exposure^2
cat("  meanF =", round(mean(exp$F), 1), " minF =", round(min(exp$F), 1), "\n")

## ---------- 2. 结局 ----------
out <- extract_outcome_data(snps = exp$SNP, outcomes = "finn-b-I9_AF", proxies = FALSE)
cat("STEP2 outcome rows =", ifelse(is.null(out), 0, nrow(out)), "\n")
if (is.null(out) || nrow(out) < 2) {
  cat("Overlap <2, falling back to ebi-a-GCST006414 (Nielsen 2018, N=1,030,836)\n")
  out <- extract_outcome_data(snps = exp$SNP, outcomes = "ebi-a-GCST006414", proxies = FALSE)
  cat("  fallback rows =", ifelse(is.null(out), 0, nrow(out)), "\n")
}
if (is.null(out) || nrow(out) < 2) { cat("STOP: outcome overlap <2\n"); quit(save = "no", status = 3) }
print(out[, c("SNP","effect_allele.outcome","other_allele.outcome","beta.outcome","se.outcome","pval.outcome","eaf.outcome","samplesize.outcome","outcome")])

## ---------- 3. harmonise ----------
dat <- harmonise_data(exposure_dat = exp, outcome_dat = out, action = 2)
cat("STEP3 harmonised rows =", nrow(dat), "\n")
print(dat[, c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure",
              "effect_allele.outcome","other_allele.outcome","eaf.outcome",
              "beta.exposure","beta.outcome","palindromic","ambiguous","mr_keep")])
dat_keep <- dat[dat$mr_keep, ]
cat("  mr_keep rows =", nrow(dat_keep), "\n")

## 单 SNP Wald ratio
cat("\n-- per-SNP Wald ratio --\n")
for (i in seq_len(nrow(dat_keep))) {
  b <- dat_keep$beta.outcome[i] / dat_keep$beta.exposure[i]
  se <- abs(dat_keep$se.outcome[i] / dat_keep$beta.exposure[i])
  cat(sprintf("  %-11s beta_exp=%7.4f  Wald_logOR=%7.4f (se %.4f)  OR=%.3f (%.3f-%.3f)  p=%.3g\n",
              dat_keep$SNP[i], dat_keep$beta.exposure[i], b, se, exp(b), exp(b-1.96*se), exp(b+1.96*se),
              2*pnorm(-abs(b/se))))
}

## ---------- 4. MR ----------
res <- mr(dat_keep, method_list = c("mr_ivw","mr_weighted_median","mr_egger_regression"))
res$OR <- exp(res$b); res$OR_lci <- exp(res$b - 1.96*res$se); res$OR_uci <- exp(res$b + 1.96*res$se)
cat("\nSTEP4 mr() results\n")
print(res[, c("method","nsnp","b","se","pval","OR","OR_lci","OR_uci")], row.names = FALSE)
write.csv(res, "mr_af/mr_results.csv", row.names = FALSE)

## IVW 随机效应（mr_ivw 默认报告固定效应；另算 RE）
ivw_re <- mr_ivw(dat_keep, model = "default", robust = FALSE)
cat("\nIVW detail:/n"); print(ivw_re[, c("method","nsnp","b","se","pval")], row.names = FALSE)

## ---------- 5. 敏感性 ----------
cat("\nSTEP5 heterogeneity\n"); het <- mr_heterogeneity(dat_keep)
print(het[, c("method","Q","Q_df","Q_pval")], row.names = FALSE)
cat("\nSTEP5 pleiotropy\n"); ple <- mr_pleiotropy_test(dat_keep)
print(ple[, c("egger_intercept","se","pval")], row.names = FALSE)
cat("\nSTEP5 leave-one-out\n"); loo <- mr_leaveoneout(dat_keep)
print(loo[, c("SNP","b","se","pval")], row.names = FALSE)
write.csv(het, "mr_af/heterogeneity.csv", row.names = FALSE)
write.csv(ple, "mr_af/pleiotropy.csv", row.names = FALSE)
write.csv(loo, "mr_af/leaveoneout.csv", row.names = FALSE)

## ---------- 6. 图 ----------
pdf("mr_af/fig_scatter.pdf", width = 7, height = 6)
print(mr_scatter_plot(res, dat_keep)[[1]]); dev.off()
pdf("mr_af/fig_forest.pdf", width = 7, height = 5)
print(mr_forest_plot(res)[[1]]); dev.off()
pdf("mr_af/fig_loo.pdf", width = 7, height = 5)
print(mr_leaveoneout_plot(loo)[[1]]); dev.off()
cat("\nSTEP6 plots: ", paste(list.files("mr_af", pattern="pdf$"), collapse=", "), "\n")
cat("DONE\n")
