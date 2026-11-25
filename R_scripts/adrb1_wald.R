dat <- readRDS("mr_adrb1/dat_harmonised.rds"); dk <- dat[dat$mr_keep, ]
wr <- dk$beta.outcome/dk$beta.exposure; sew <- abs(dk$se.outcome/dk$beta.exposure)
d <- data.frame(SNP=dk$SNP, b_exp=round(dk$beta.exposure,3), b_out=round(dk$beta.outcome,4),
                Wald=round(wr,4), OR=round(exp(wr),3), lci=round(exp(wr-1.96*sew),3), uci=round(exp(wr+1.96*sew),3),
                p=signif(2*pnorm(-abs(wr/sew)),3), F=round(dk$F,1))
print(d, row.names=FALSE)
cat("n_neg=", sum(wr<0), " n_pos=", sum(wr>0), "\n")
write.csv(d, "mr_adrb1/wald_per_snp.csv", row.names=FALSE)
print(data.frame(SNP=dk$SNP, eaf=round(dk$eaf.exposure,3), beta=round(dk$beta.exposure,3), se=round(dk$se.exposure,3),
                 p=signif(dk$pval.exposure,3), F=round(dk$F,1)))
