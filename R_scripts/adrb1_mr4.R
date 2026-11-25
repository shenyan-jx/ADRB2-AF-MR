suppressMessages(library(TwoSampleMR))
dat <- readRDS("mr_adrb1/dat_harmonised.rds"); dk <- dat[dat$mr_keep, ]
res <- mr(dk, method_list = c("mr_ivw","mr_weighted_median","mr_egger_regression"))
loo <- mr_leaveoneout(dk)
cat("loo cols:", paste(names(loo), collapse=" | "), "\n")
loo$OR <- exp(loo$b); loo$OR_lci <- exp(loo$b-1.96*loo$se); loo$OR_uci <- exp(loo$b+1.96*loo$se)
pc <- if ("pval" %in% names(loo)) "pval" else "p"
print(loo[, c("SNP","b","se",pc,"OR","OR_lci","OR_uci")], row.names=FALSE)
write.csv(loo, "mr_adrb1/leaveoneout.csv", row.names=FALSE)
pdf("mr_adrb1/fig1_scatter.pdf", width=7, height=6); print(mr_scatter_plot(res, dk)[[1]]); dev.off()
pdf("mr_adrb1/fig2_forest.pdf", width=7, height=5); print(mr_forest_plot(res)[[1]]); dev.off()
pdf("mr_adrb1/fig3_loo.pdf", width=7, height=5); print(mr_leaveoneout_plot(loo)[[1]]); dev.off()
cat("plots:", paste(list.files("mr_adrb1", pattern="pdf$"), collapse=", "), "\nDONE\n")
