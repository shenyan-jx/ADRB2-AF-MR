suppressMessages(library(TwoSampleMR))
exp <- readRDS("exp_ADRB2.rds")
exp$F <- exp$beta.exposure^2 / exp$se.exposure^2

## 结局数据：来自 23:40 对 finn-b-I9_AF 的成功提取（原文见 mr_af.R 输出）
out <- data.frame(
  SNP = c("rs2082395", "rs72829145", "rs17696736"),
  effect_allele.outcome = c("G", "C", "G"),
  other_allele.outcome  = c("A", "T", "A"),
  beta.outcome  = c(-0.0205, 0.0099, 0.0646),
  se.outcome    = c(0.0152, 0.0300, 0.0151),
  pval.outcome  = c(0.1783, 0.740901, 0.00001893),
  eaf.outcome   = c(0.6193, 0.06604, 0.3968),
  outcome = "Atrial fibrillation and flutter || id:finn-b-I9_AF",
  id.outcome = "finn-b-I9_AF",
  stringsAsFactors = FALSE
)
saveRDS(out, "mr_af/outcome_AF.rds")

dat <- harmonise_data(exp, out, action = 2)
dat_keep <- dat[dat$mr_keep, ]
saveRDS(dat, "mr_af/dat_harmonised.rds")
cat("overlap =", nrow(dat), " mr_keep =", nrow(dat_keep), "\n")

res <- mr(dat_keep, method_list = c("mr_ivw","mr_weighted_median","mr_egger_regression"))
res$OR <- exp(res$b); res$OR_lci <- exp(res$b - 1.96*res$se); res$OR_uci <- exp(res$b + 1.96*res$se)
print(res[, c("method","nsnp","b","se","pval","OR","OR_lci","OR_uci")], row.names = FALSE)
write.csv(res, "mr_af/mr_results.csv", row.names = FALSE)

het <- mr_heterogeneity(dat_keep); print(het[, c("method","Q","Q_df","Q_pval")], row.names = FALSE)
ple <- mr_pleiotropy_test(dat_keep); print(ple[, c("egger_intercept","se","pval")], row.names = FALSE)
loo <- mr_leaveoneout(dat_keep)
loo$OR <- exp(loo$b); loo$OR_lci <- exp(loo$b-1.96*loo$se); loo$OR_uci <- exp(loo$b+1.96*loo$se)
cat("\n-- leave-one-out --\n"); print(loo[, c("SNP","b","se","pval","OR","OR_lci","OR_uci")], row.names = FALSE)
write.csv(het, "mr_af/heterogeneity.csv", row.names = FALSE)
write.csv(ple, "mr_af/pleiotropy.csv", row.names = FALSE)
write.csv(loo, "mr_af/leaveoneout.csv", row.names = FALSE)

q <- het$Q[het$method == "Inverse variance weighted"]; qdf <- het$Q_df[het$method == "Inverse variance weighted"]
ivw <- res[res$method == "Inverse variance weighted", ]
scale <- max(1, q/qdf); se_re <- ivw$se*sqrt(scale); z <- ivw$b/se_re
cat(sprintf("\nIVW-RE: OR=%.3f (%.3f-%.3f) p=%.3f   [Q=%.2f df=%d Qp=%.3f]\n",
    exp(ivw$b), exp(ivw$b-1.96*se_re), exp(ivw$b+1.96*se_re), 2*pnorm(-abs(z)), q, qdf, 1-pchisq(q,qdf)))

wr <- dat_keep$beta.outcome/dat_keep$beta.exposure
se_w <- abs(dat_keep$se.outcome/dat_keep$beta.exposure)
df_w <- data.frame(SNP=dat_keep$SNP, Wald_logOR=round(wr,4), OR=round(exp(wr),3),
                   lci=round(exp(wr-1.96*se_w),3), uci=round(exp(wr+1.96*se_w),3),
                   p=signif(2*pnorm(-abs(wr/se_w)),3), F=round(dat_keep$F,1))
cat("\n-- per-SNP Wald --\n"); print(df_w, row.names = FALSE)
cat("signs:", paste(ifelse(wr<0,"neg","pos"), collapse=", "), " n_neg=", sum(wr<0), " n_pos=", sum(wr>0), "\n")
write.csv(df_w, "mr_af/wald_per_snp.csv", row.names = FALSE)

pdf("mr_af/fig1_scatter.pdf", width=7, height=6); print(mr_scatter_plot(res, dat_keep)[[1]]); dev.off()
pdf("mr_af/fig2_forest.pdf", width=7, height=5); print(mr_forest_plot(res)[[1]]); dev.off()
pdf("mr_af/fig3_loo.pdf", width=7, height=5); print(mr_leaveoneout_plot(loo)[[1]]); dev.off()
cat("\nplots:", paste(list.files("mr_af", pattern="pdf$"), collapse=", "), "\nDONE\n")
