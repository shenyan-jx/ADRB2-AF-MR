suppressMessages({library(TwoSampleMR); library(dplyr)})
options(width = 200)
x <- readRDS("_mr_results.rds")

cat("=========== Cochran Q 异质性 ===========\n")
if (!is.null(x$het) && nrow(x$het)) {
  print(as.data.frame(x$het %>% filter(method == "Inverse variance weighted") %>%
    select(exposure_set, outcome_set, Q, Q_df, Q_pval) %>%
    mutate(across(where(is.numeric), ~round(.x, 3))) %>% arrange(Q_pval)), row.names = FALSE)
}

cat("\n=========== 留一法：VEGFA_eqtl ===========\n")
for (on in c("CHD", "MI", "heart_fail", "hypertension", "atrial_fib", "VTE")) {
  d <- x$detail[[paste("VEGFA_eqtl", on)]]
  if (is.null(d)) next
  loo <- tryCatch(suppressMessages(mr_leaveoneout(d)), error = function(e) NULL)
  if (is.null(loo)) next
  cat("\n--", on, "--\n")
  print(as.data.frame(loo %>% select(SNP, b, se, p) %>%
    mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)
}

cat("\n=========== Steiger 方向性检验 ===========\n")
st <- list()
for (k in names(x$detail)) {
  d <- x$detail[[k]]
  s <- tryCatch(suppressMessages(directionality_test(d)), error = function(e) NULL)
  if (is.null(s)) next
  s$set <- k; st[[k]] <- s
}
if (length(st)) {
  s <- bind_rows(st) %>% select(set, snp_r2.exposure, snp_r2.outcome, correct_causal_direction, steiger_pval)
  print(as.data.frame(s %>% mutate(across(where(is.numeric), ~round(.x, 5)))), row.names = FALSE)
}
