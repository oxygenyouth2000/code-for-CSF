library(tidyverse)
library(openxlsx)
library(grf)
library(scales)

set.seed(2026)

current_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
project_dir <- if (basename(current_dir) == "code_all") dirname(current_dir) else current_dir

df_overall_path <- file.path(project_dir, "df_overall.xlsx")
csf_fit_path <- file.path(project_dir, "results data", "csf_fit.rds")
toc_curve_output_path <- file.path(project_dir, "toc_curve_df.csv")
autoc_output_path <- file.path(project_dir, "autoc_summary.csv")
toc_plot_output_path <- file.path(project_dir, "TOC_Curve.png")

df_overall <- read.xlsx(df_overall_path, sheet = 1)

csf_fit <- readRDS(csf_fit_path)

priority_scores <- predict(csf_fit)$predictions

csf_fit$W.hat <- pmin(pmax(csf_fit$W.hat, 0.01), 0.99)

debiasing_weights <- with(
  csf_fit,
  (W.orig - W.hat) / (W.hat * (1 - W.hat))
)

q_grid <- c(0.01, seq(0.05, 0.95, by = 0.05), 1)

autoc_result <- rank_average_treatment_effect(
  forest = csf_fit,
  priorities = priority_scores,
  target = "AUTOC",
  q = q_grid,
  R = 1000,
  debiasing.weights = debiasing_weights
)

autoc_value <- unname(autoc_result$estimate[[1]])

autoc_se <- unname(autoc_result$std.err[[1]])

autoc_ci_lower <- autoc_value - 1.96 * autoc_se

autoc_ci_upper <- autoc_value + 1.96 * autoc_se

autoc_p_value <- if (is.na(autoc_se) || autoc_se == 0) NA_real_ else 2 * pnorm(-abs(autoc_value / autoc_se))

p_text <- case_when(
  is.na(autoc_p_value) ~ "P = NA",
  autoc_p_value < 0.001 ~ "P < 0.001",
  TRUE ~ paste0("P = ", sprintf("%.3f", autoc_p_value))
)

toc_curve_df <- autoc_result$TOC %>%
  transmute(
    q = q,
    toc = estimate,
    toc_se = std.err,
    toc_lower = toc - 1.96 * toc_se,
    toc_upper = toc + 1.96 * toc_se,
    n_selected = ceiling(length(priority_scores) * q)
  )

autoc_label <- paste0(
  "AUTOC = ", sprintf("%.3f", autoc_value), "\n",
  "95% CI: ", sprintf("%.3f", autoc_ci_lower), " to ", sprintf("%.3f", autoc_ci_upper), "\n",
  p_text
)

label_y <- quantile(toc_curve_df$toc_upper, probs = 0.25, na.rm = TRUE)

toc_plot <- ggplot(toc_curve_df, aes(x = q, y = toc)) +
  geom_ribbon(aes(ymin = toc_lower, ymax = toc_upper), fill = "#9EC9E2", alpha = 0.25) +
  geom_area(aes(y = pmax(toc, 0)), fill = "#9EC9E2", alpha = 0.70) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 1.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  scale_x_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.2),
    labels = label_percent(accuracy = 1),
    expand = c(0, 0)
  ) +
  labs(
    x = "Proportion of patients prioritized for BB treatment",
    y = "Difference in RMST (month)"
  ) +
  annotate(
    geom = "label",
    x = 0.22,
    y = label_y,
    label = autoc_label,
    hjust = 0,
    vjust = 1,
    size = 5.3,
    linewidth = 0,
    fill = "white",
    alpha = 0.90
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16)
  )

autoc_summary_df <- tibble(
  metric = "AUTOC",
  estimate = autoc_value,
  std_error = autoc_se,
  ci_lower = autoc_ci_lower,
  ci_upper = autoc_ci_upper,
  p_value = autoc_p_value,
  n_total = length(priority_scores)
)

write_csv(toc_curve_df, toc_curve_output_path)
write_csv(autoc_summary_df, autoc_output_path)
print(toc_plot)
ggsave(
  filename = toc_plot_output_path,
  plot = toc_plot,
  width = 8,
  height = 6,
  dpi = 600,
  bg = "white"
)
