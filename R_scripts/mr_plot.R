suppressMessages({library(ggplot2); library(dplyr)})
options(width = 200)
m <- read.csv("MR_main_results.csv", fileEncoding = "UTF-8", stringsAsFactors = FALSE)

lab_exp <- c(VEGFA_eqtl = "VEGFA (eQTL)", VEGFA_pqtl_soma = "VEGFA (pQTL)",
             FLT1_eqtl = "FLT1/VEGFR1 (eQTL)", FLT4_eqtl = "FLT4/VEGFR3 (eQTL)",
             FLT4_pqtl = "FLT4/VEGFR3 (pQTL)", KDR_pqtl = "KDR/VEGFR2 (pQTL)",
             PDGFRB_eqtl = "PDGFRB (eQTL)")
lab_out <- c(hypertension = "Hypertension", CHD = "Coronary heart disease", MI = "Myocardial infarction",
             heart_fail = "Heart failure", atrial_fib = "Atrial fibrillation", VTE = "Venous thromboembolism")
m <- m %>% mutate(exp = lab_exp[exposure_set], out = lab_out[outcome_set],
                  sig = ifelse(P_FDR < 0.05, "FDR < 0.05",
                        ifelse(p < 0.05, "p < 0.05", "NS")))
m$sig <- factor(m$sig, levels = c("FDR < 0.05", "p < 0.05", "NS"))
m$exp <- factor(m$exp, levels = rev(c("VEGFA (eQTL)", "VEGFA (pQTL)", "KDR/VEGFR2 (pQTL)",
                                      "FLT1/VEGFR1 (eQTL)", "FLT4/VEGFR3 (eQTL)", "FLT4/VEGFR3 (pQTL)",
                                      "PDGFRB (eQTL)")))
m$out <- factor(m$out, levels = rev(c("Coronary heart disease", "Myocardial infarction", "Heart failure",
                                      "Hypertension", "Atrial fibrillation", "Venous thromboembolism")))

p <- ggplot(m, aes(x = OR, y = out, colour = sig, shape = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.22, linewidth = 0.55) +
  geom_point(size = 2.4) +
  facet_grid(exp ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_log10(breaks = c(0.6, 0.8, 1.0, 1.25, 1.6)) +
  scale_colour_manual(values = c("FDR < 0.05" = "#C0392B", "p < 0.05" = "#E67E22", "NS" = "#7F8C8D")) +
  scale_shape_manual(values = c("FDR < 0.05" = 18, "p < 0.05" = 16, "NS" = 16)) +
  labs(x = "Odds ratio per genetically predicted increase in protein/gene expression (95% CI)",
       y = NULL, colour = NULL, shape = NULL,
       title = "Drug-target Mendelian randomisation: angiogenic pathway proteins and cardiovascular outcomes",
       subtitle = "cis instruments (p<5e-8, r2<0.001); outcomes FinnGen; Wald ratio where n=1") +
  theme_bw(base_size = 10) +
  theme(strip.placement = "outside", strip.background = element_rect(fill = "grey94", colour = "grey60"),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8.5, face = "bold"),
        panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        legend.position = "bottom", legend.key.width = unit(1.2, "lines"),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        plot.title.position = "plot")

ggsave("fig_mr_forest.png", p, width = 8.2, height = 7.6, dpi = 300, bg = "white")
cat("已生成 fig_mr_forest.png\n")

# 附表：标记 FDR 显著
top <- m %>% filter(P_FDR < 0.10) %>% arrange(p) %>%
  transmute(Exposure = exp, Outcome = out, n = nsnp_used,
            OR = round(OR, 3), CI = sprintf("%.3f-%.3f", CI_low, CI_high),
            p = signif(p, 3), FDR = round(P_FDR, 3))
print(as.data.frame(top), row.names = FALSE)
