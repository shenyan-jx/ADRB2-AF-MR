suppressMessages(library(TwoSampleMR))
dir.create("mr_adrb1", showWarnings = FALSE)

exp <- extract_instruments(outcomes = "eqtl-a-ENSG00000043591",
                           p1 = 5e-6, clump = TRUE, r2 = 0.001, kb = 10000)
cat("STEP1 n_instruments =", nrow(exp), "\n")
cat("SNPs:", paste(exp$SNP, collapse=", "), "\n")
if (nrow(exp) < 3) { cat("STOP: insufficient instruments\n"); quit(save="no", status=2) }
exp$F <- exp$beta.exposure^2 / exp$se.exposure^2
cat("meanF =", round(mean(exp$F),1), " minF =", round(min(exp$F),1), "\n")
saveRDS(exp, "exp_ADRB1.rds")
print(as.data.frame(exp[, c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure",
                            "beta.exposure","se.exposure","pval.exposure")]), row.names=FALSE)
