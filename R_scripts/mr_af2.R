suppressMessages(library(TwoSampleMR))
dir.create("mr_af", showWarnings = FALSE)
exp <- readRDS("exp_ADRB2.rds")
exp$F <- exp$beta.exposure^2 / exp$se.exposure^2

out <- extract_outcome_data(snps = exp$SNP, outcomes = "finn-b-I9_AF", proxies = FALSE)
saveRDS(out, "mr_af/outcome_AF.rds")
dat <- harmonise_data(exp, out, action = 2)
saveRDS(dat, "mr_af/dat_harmonised.rds")
dat_keep <- dat[dat$mr_keep, ]
cat("overlap =", nrow(dat), " mr_keep =", nrow(dat_keep), "\n")

res <- mr(dat_keep, method_list = c("mr_ivw","mr_weighted_median","mr_egger_regression"))
res$OR <- exp(res$b); res$OR_lci <- exp(res$b - 1.96*res$se); res$OR_uci <- exp(res$b + 1.96*res$se)
print(res[, c("method","nsnp","b","se","pval","OR","OR_lci","OR_uci")], row.names = FALSE)
write.csv(res, "mr_af/mr_results.csv", row.names = FALSE)

het <- mr_heterogeneity(dat_keep)
cat("\n-- heterogeneity --\n"); print(het[, c("method","Q","Q_df","Q_pval")], row.names = FALSE)
ple <- mr_pleiotropy_test(dat_keep)
cat("\n-- pleiotropy (Egger intercept) --\n"); print(ple[, c("egger_intercept","se","pval")], row.names = FALSE)
loo <- mr_leaveoneout(dat_keep)
cat("\n-- leave-one-out --\n")
loo$OR <- exp(loo$b); loo$OR_lci <- exp(loo$b-1.96*loo$se); loo$OR_uci <- exp(loo$b+1.96*loo$se)
print(loo[, c("SNP","b","se","pval","OR","OR_lci","OR_uci")], row.names = FALSE)
write.csv(het, "mr_af/heterogeneity.csv", row.names = FALSE)
write.csv(ple, "mr_af/pleiotropy.csv", row.names = FALSE)
write.csv(loo, "mr_af/leaveoneout.csv", row.names = FALSE)

## IVW 随机效应（Q 校正）
q <- het$Q[het$method == "Inverse variance weighted"]; qdf <- het$Q_df[het$method == "Inverse variance weighted"]
ivw <- res[res$method == "Inverse variance weighted", ]
scale <- max(1, q / qdf)
se_re <- ivw$se * sqrt(scale)
z <- ivw$b / se_re
cat(sprintf("\nIVW random-effects: logOR=%.4f se=%.4f OR=%.3f (%.3f-%.3f) p=%.3f  [Q=%.2f, df=%d, p=%.2f, scale=%.2f]\n",
            ivw$b, se_re, exp(ivw$b), exp(ivw$b-1.96*se_re), exp(ivw$b+1.96*se_re), 2*pnorm(-abs(z)), q, qdf,
            1-pchisq(q, qdf), scale))

## 单 SNP Wald 与方向一致性
cat("\n-- per-SNP Wald ratios --\n")
wr <- dat_keep$beta.outcome / dat_keep$beta.exposure
se_w <- abs(dat_keep$se.outcome / dat_keep$beta.exposure)
df_w <- data.frame(SNP = dat_keep$SNP, Wald_logOR = round(wr,4), OR = round(exp(wr),3),
                   lci = round(exp(wr-1.96*se_w),3), uci = round(exp(wr+1.96*se_w),3),
                   p = signif(2*pnorm(-abs(wr/se_w)), 3), F = round(dat_keep$F, 1))
print(df_w, row.names = FALSE)
cat("direction signs:", paste(ifelse(wr < 0, "neg", "pos"), collapse = ", "), "\n")
cat("n_neg =", sum(wr < 0), " n_pos =", sum(wr > 0), "\n")
write.csv(df_w, "mr_af/wald_per_snp.csv", row.names = FALSE)

## 图
pdf("mr_af/fig1_scatter.pdf", width = 7, height = 6); print(mr_scatter_plot(res, dat_keep)[[1]]); dev.off()
pdf("mr_af/fig2_forest.pdf", width = 7, height = 5); print(mr_forest_plot(res)[[1]]); dev.off()
pdf("mr_af/fig3_loo.pdf", width = 7, height = 5); print(mr_leaveoneout_plot(loo)[[1]]); dev.off()
cat("\nplots:", paste(list.files("mr_af", pattern = "pdf$"), collapse = ", "), "\n")
cat("DONE\n")
