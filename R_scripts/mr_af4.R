suppressMessages(library(TwoSampleMR))
dat <- readRDS("mr_af/dat_harmonised.rds"); dat_keep <- dat[dat$mr_keep, ]
res <- mr(dat_keep, method_list = c("mr_ivw","mr_weighted_median","mr_egger_regression"))
res$OR <- exp(res$b); res$OR_lci <- exp(res$b-1.96*res$se); res$OR_uci <- exp(res$b+1.96*res$se)
het <- mr_heterogeneity(dat_keep)

loo <- mr_leaveoneout(dat_keep)
cat("loo names:", paste(names(loo), collapse=" | "), "\n")
loo$OR <- exp(loo$b); loo$OR_lci <- exp(loo$b-1.96*loo$se); loo$OR_uci <- exp(loo$b+1.96*loo$se)
keep_cols <- intersect(c("SNP","b","se","pval","OR","OR_lci","OR_uci"), names(loo))
print(loo[, keep_cols], row.names = FALSE)
write.csv(loo, "mr_af/leaveoneout.csv", row.names = FALSE)

q <- het$Q[het$method=="Inverse variance weighted"]; qdf <- het$Q_df[het$method=="Inverse variance weighted"]
ivw <- res[res$method=="Inverse variance weighted",]
scale <- max(1, q/qdf); se_re <- ivw$se*sqrt(scale)
cat(sprintf("\nIVW fixed : OR=%.3f (%.3f-%.3f) p=%.3f\n", exp(ivw$b), exp(ivw$b-1.96*ivw$se), exp(ivw$b+1.96*ivw$se), ivw$pval))
cat(sprintf("IVW random: OR=%.3f (%.3f-%.3f) p=%.3f  [scale=%.2f]\n",
    exp(ivw$b), exp(ivw$b-1.96*se_re), exp(ivw$b+1.96*se_re), 2*pnorm(-abs(ivw$b/se_re)), scale))

pdf("mr_af/fig1_scatter.pdf", width=7, height=6); print(mr_scatter_plot(res, dat_keep)[[1]]); dev.off()
pdf("mr_af/fig2_forest.pdf", width=7, height=5); print(mr_forest_plot(res)[[1]]); dev.off()
p <- try(mr_leaveoneout_plot(loo), silent=TRUE)
if (!inherits(p,"try-error")) { pdf("mr_af/fig3_loo.pdf", width=7, height=5); print(p[[1]]); dev.off() } else cat("loo plot err:", conditionMessage(attr(p,"condition")), "\n")
cat("\nfiles:", paste(list.files("mr_af", pattern="pdf$"), collapse=", "), "\n")
