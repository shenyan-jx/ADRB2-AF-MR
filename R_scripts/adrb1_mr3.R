suppressMessages(library(TwoSampleMR))
dat <- readRDS("mr_adrb1/dat_harmonised.rds"); dk <- dat[dat$mr_keep, ]
res <- mr(dk, method_list = c("mr_ivw","mr_weighted_median","mr_egger_regression"))
res$OR <- exp(res$b); res$OR_lci <- exp(res$b-1.96*res$se); res$OR_uci <- exp(res$b+1.96*res$se)
cat("== mr results ==\n"); print(res[, c("method","nsnp","b","se","pval","OR","OR_lci","OR_uci")], row.names=FALSE)
write.csv(res, "mr_adrb1/mr_results.csv", row.names=FALSE)

het <- mr_heterogeneity(dk); cat("\n== heterogeneity ==\n"); print(het[, c("method","Q","Q_df","Q_pval")], row.names=FALSE)
ple <- mr_pleiotropy_test(dk); cat("\n== pleiotropy ==\n"); print(ple[, c("egger_intercept","se","pval")], row.names=FALSE)
loo <- mr_leaveoneout(dk)
loo$OR <- exp(loo$b); loo$OR_lci <- exp(loo$b-1.96*loo$se); loo$OR_uci <- exp(loo$b+1.96*loo$se)
cat("\n== leave-one-out ==\n"); print(loo[, c("SNP","b","se","pval","OR","OR_lci","OR_uci")], row.names=FALSE)
write.csv(het,"mr_adrb1/heterogeneity.csv",row.names=FALSE)
write.csv(ple,"mr_adrb1/pleiotropy.csv",row.names=FALSE)
write.csv(loo,"mr_adrb1/leaveoneout.csv",row.names=FALSE)

wr <- dk$beta.outcome/dk$beta.exposure; sew <- abs(dk$se.outcome/dk$beta.exposure)
dfw <- data.frame(SNP=dk$SNP, Wald_logOR=round(wr,4), OR=round(exp(wr),3),
                  lci=round(exp(wr-1.96*sew),3), uci=round(exp(wr+1.96*sew),3),
                  p=signif(2*pnorm(-abs(wr/sew)),3), F=round(dk$F,1))
cat("\n== per-SNP Wald ==\n"); print(dfw, row.names=FALSE)
cat("signs:", paste(ifelse(wr<0,"neg","pos"),collapse=", "), " n_neg=",sum(wr<0)," n_pos=",sum(wr>0),"\n")
write.csv(dfw,"mr_adrb1/wald_per_snp.csv",row.names=FALSE)

q <- het$Q[het$method=="Inverse variance weighted"]; qdf <- het$Q_df[het$method=="Inverse variance weighted"]
ivw <- res[res$method=="Inverse variance weighted",]
sc <- max(1,q/qdf); se_re <- ivw$se*sqrt(sc)
cat(sprintf("\nIVW fixed : OR=%.3f (%.3f-%.3f) p=%.3f\n", exp(ivw$b), exp(ivw$b-1.96*ivw$se), exp(ivw$b+1.96*ivw$se), ivw$pval))
cat(sprintf("IVW random: OR=%.3f (%.3f-%.3f) p=%.3f [Q=%.2f df=%d Qp=%.4f scale=%.2f]\n",
    exp(ivw$b), exp(ivw$b-1.96*se_re), exp(ivw$b+1.96*se_re), 2*pnorm(-abs(ivw$b/se_re)), q, qdf, 1-pchisq(q,qdf), sc))

pdf("mr_adrb1/fig1_scatter.pdf", width=7, height=6); print(mr_scatter_plot(res, dk)[[1]]); dev.off()
pdf("mr_adrb1/fig2_forest.pdf", width=7, height=5); print(mr_forest_plot(res)[[1]]); dev.off()
pdf("mr_adrb1/fig3_loo.pdf", width=7, height=5); print(mr_leaveoneout_plot(loo)[[1]]); dev.off()
cat("\nplots:", paste(list.files("mr_adrb1", pattern="pdf$"), collapse=", "), "\nDONE\n")
