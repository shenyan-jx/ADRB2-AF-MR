suppressMessages({library(TwoSampleMR); library(dplyr)})
options(width = 220)
instr   <- readRDS("_instruments.rds")
out_all <- readRDS("_outcomes.rds")

mk <- function(d, nm) data.frame(SNP = d$rsid, beta = d$beta, se = d$se,
                                 ea = d$ea, nea = d$nea, eaf = d$eaf, pval = d$p,
                                 phenotype = nm, stringsAsFactors = FALSE)

rows <- list(); het <- list(); plei <- list(); det <- list()
for (en in names(instr)) {
  th <- instr[[en]]
  ex <- format_data(mk(th, en), type = "exposure", snp_col = "SNP", beta_col = "beta",
                    se_col = "se", effect_allele_col = "ea", other_allele_col = "nea",
                    eaf_col = "eaf", pval_col = "pval", phenotype_col = "phenotype")
  for (on in unique(out_all$outcome)) {
    o <- out_all %>% filter(outcome == on, rsid %in% th$rsid)
    if (!nrow(o)) next
    ou <- format_data(mk(o, on), type = "outcome", snp_col = "SNP", beta_col = "beta",
                      se_col = "se", effect_allele_col = "ea", other_allele_col = "nea",
                      eaf_col = "eaf", pval_col = "pval", phenotype_col = "phenotype")
    d <- tryCatch(suppressMessages(harmonise_data(ex, ou, action = 2)), error = function(e) NULL)
    if (is.null(d)) next
    d <- d[d$mr_keep, , drop = FALSE]
    n <- nrow(d); if (n == 0) next
    meth <- if (n == 1) "mr_wald_ratio" else if (n == 2) c("mr_ivw", "mr_weighted_median") else
      c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
    r <- tryCatch(suppressMessages(mr(d, method_list = meth)), error = function(e) NULL)
    if (is.null(r)) next
    r$exposure_set <- en; r$outcome_set <- on; r$nsnp_used <- n
    rows[[paste(en, on)]] <- r
    det[[paste(en, on)]] <- d
    if (n >= 2) {
      h <- tryCatch(suppressMessages(mr_heterogeneity(d)), error = function(e) NULL)
      if (!is.null(h)) { h$exposure_set <- en; h$outcome_set <- on; het[[paste(en,on)]] <- h }
    }
    if (n >= 3) {
      p <- tryCatch(suppressMessages(mr_pleiotropy_test(d)), error = function(e) NULL)
      if (!is.null(p)) { p$exposure_set <- en; p$outcome_set <- on; plei[[paste(en,on)]] <- p }
    }
  }
}
res <- bind_rows(rows)
cat("完成的暴露×结局组合:", length(unique(paste(res$exposure_set, res$outcome_set))), "\n")
cat("MR 结果行数:", nrow(res), "\n\n")

res <- res %>% mutate(OR = exp(b), CI_low = exp(b - 1.96*se), CI_high = exp(b + 1.96*se),
                      p = as.numeric(pval))

main <- res %>% filter(method %in% c("Inverse variance weighted", "Wald ratio")) %>%
  mutate(P_FDR = p.adjust(p, "BH")) %>%
  select(exposure_set, outcome_set, method, nsnp_used, b, se, p, P_FDR, OR, CI_low, CI_high) %>%
  arrange(p)
cat("=========== 主分析（IVW / Wald ratio）===========\n")
print(as.data.frame(main %>% mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)

cat("\n=========== 敏感性分析 ===========\n")
sens <- res %>% filter(!method %in% c("Inverse variance weighted", "Wald ratio")) %>%
  select(exposure_set, outcome_set, method, nsnp_used, b, p) %>%
  mutate(across(c(b, p), ~round(.x, 4)))
print(as.data.frame(sens), row.names = FALSE)

if (length(plei)) {
  cat("\n=========== MR-Egger 截距（多效性检验）===========\n")
  print(as.data.frame(bind_rows(plei) %>% select(exposure_set, outcome_set, egger_intercept, se, pval) %>%
                       mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)
}
saveRDS(list(res = res, main = main, detail = det,
             het = bind_rows(het), plei = bind_rows(plei)), "_mr_results.rds")
write.csv(main, "MR_main_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(res,  "MR_all_methods.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("\n已保存\n")
