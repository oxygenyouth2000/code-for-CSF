library(tidyverse)
library(openxlsx)
library(grf)
library(scales)

# 统一随机种子，保证 bootstrap 标准误结果可复现。
set.seed(2026)

# 兼容从项目根目录或 `code_all` 目录直接运行脚本。
current_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
project_dir <- if (basename(current_dir) == "code_all") dirname(current_dir) else current_dir

# 明确输入与输出路径，避免受当前工作目录影响。
df_overall_path <- file.path(project_dir, "df_overall.xlsx")
csf_fit_path <- file.path(project_dir, "中间分析结果数据", "csf_fit.rds")
toc_curve_output_path <- file.path(project_dir, "toc_curve_df.csv")
autoc_output_path <- file.path(project_dir, "autoc_summary.csv")
toc_plot_output_path <- file.path(project_dir, "TOC_Curve.png")

# 读取已保存的全样本分析数据。
df_overall <- read.xlsx(df_overall_path, sheet = 1)

# 读取已经建好的 CSF 模型对象。
csf_fit <- readRDS(csf_fit_path)

# 直接使用已建好的 CSF 模型预测值作为 treatment prioritization scores。
# 这里采用全样本模型性能评价口径，不再随机二分样本，也不重新训练 forest。
priority_scores <- predict(csf_fit)$predictions

# 校验 priorities 与 forest 样本量是否一致，避免错位。
if (length(priority_scores) != length(csf_fit$Y.orig)) {
  stop("priority_scores 与 csf_fit 样本量不一致，无法计算 TOC/AUTOC。")
}

# 若 `df_overall.xlsx` 中也保存了 ITE，则额外校验长度是否一致。
if ("ITE" %in% names(df_overall) && nrow(df_overall) != length(priority_scores)) {
  stop("df_overall.xlsx 与 csf_fit.rds 的样本量不一致，请先核对上游导出结果。")
}

# 为避免 `rank_average_treatment_effect()` 因极端 propensity 报错，
# 将 forest 内部的 W.hat 截断到 [0.01, 0.99]。
csf_fit$W.hat <- pmin(pmax(csf_fit$W.hat, 0.01), 0.99)

# 按 grf 二值处理公式手动构造 debiasing weights，
# 这样可以在保留现有模型的前提下顺利完成 RATE 估计。
debiasing_weights <- with(
  csf_fit,
  (W.orig - W.hat) / (W.hat * (1 - W.hat))
)

# 设定 TOC 曲线计算网格，并保证最后一个点为 1 以满足 grf 接口要求。
q_grid <- c(0.01, seq(0.05, 0.95, by = 0.05), 1)

# 基于已建好的 CSF 模型直接计算 AUTOC 与 TOC 曲线。
autoc_result <- rank_average_treatment_effect(
  forest = csf_fit,
  priorities = priority_scores,
  target = "AUTOC",
  q = q_grid,
  R = 1000,
  debiasing.weights = debiasing_weights
)

# 提取 AUTOC 点估计。
autoc_value <- unname(autoc_result$estimate[[1]])

# 提取 AUTOC 标准误。
autoc_se <- unname(autoc_result$std.err[[1]])

# 计算 95% CI 下限。
autoc_ci_lower <- autoc_value - 1.96 * autoc_se

# 计算 95% CI 上限。
autoc_ci_upper <- autoc_value + 1.96 * autoc_se

# 计算双侧 P 值。
autoc_p_value <- if (is.na(autoc_se) || autoc_se == 0) NA_real_ else 2 * pnorm(-abs(autoc_value / autoc_se))

# 将 P 值整理成图中标注文本。
p_text <- case_when(
  is.na(autoc_p_value) ~ "P = NA",
  autoc_p_value < 0.001 ~ "P < 0.001",
  TRUE ~ paste0("P = ", sprintf("%.3f", autoc_p_value))
)

# 提取 grf 返回的 TOC 数据框，并补充标准误与 95% CI。
toc_curve_df <- autoc_result$TOC %>%
  transmute(
    q = q,
    toc = estimate,
    toc_se = std.err,
    toc_lower = toc - 1.96 * toc_se,
    toc_upper = toc + 1.96 * toc_se,
    n_selected = ceiling(length(priority_scores) * q)
  )

# 构造 AUTOC 结果标签。
autoc_label <- paste0(
  "AUTOC = ", sprintf("%.3f", autoc_value), "\n",
  "95% CI: ", sprintf("%.3f", autoc_ci_lower), " to ", sprintf("%.3f", autoc_ci_upper), "\n",
  p_text
)

# 让标签尽量落在图内空白区域。
label_y <- quantile(toc_curve_df$toc_upper, probs = 0.25, na.rm = TRUE)

# 使用 grf 返回的标准 TOC 点估计与区间直接作图。
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

# 汇总 AUTOC 结果，便于论文撰写与复核。
autoc_summary_df <- tibble(
  metric = "AUTOC",
  estimate = autoc_value,
  std_error = autoc_se,
  ci_lower = autoc_ci_lower,
  ci_upper = autoc_ci_upper,
  p_value = autoc_p_value,
  n_total = length(priority_scores)
)

# 输出 TOC 曲线数据。
write_csv(toc_curve_df, toc_curve_output_path)

# 输出 AUTOC 汇总结果。
write_csv(autoc_summary_df, autoc_output_path)

# 在交互环境中直接显示图形。
print(toc_plot)

# 保存高分辨率 TOC 曲线图。
ggsave(
  filename = toc_plot_output_path,
  plot = toc_plot,
  width = 8,
  height = 6,
  dpi = 600,
  bg = "white"
)
