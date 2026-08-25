library(tidyverse)
library(openxlsx)
library(readr)

# 读取calibration汇总结果
calib_summary_df <- read.xlsx("分析结果数据/calib_plot_df.xlsx") %>%
  as_tibble() %>%
  rename(tertile = tritile) %>%
  mutate(
    tertile = factor(tertile, levels = c("Q1", "Q2", "Q3"))
  )

# 使用现成的预测值与观察值及其区间，直接绘制3个tertile的calibration图
calib_plot_summary <- ggplot(
  calib_summary_df,
  aes(x = mean_pred, y = obs_diff, color = tertile, group = 1)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "gray55",
    linewidth = 0.5
  ) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(size = 3.2) +
  scale_color_manual(
    values = c("Q1" = "steelblue", "Q2" = "gray50", "Q3" = "red"),
    labels = c("Q1" = "Lowest tertile", "Q2" = "Middle tertile", "Q3" = "Highest tertile")
  ) +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 2)) +
  labs(
    x = "Estimated individualized treatment effects",
    y = "Observed RMST difference by tertiles",
    color = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.45, color = "gray70"),
    panel.grid.minor.y = element_line(linewidth = 0.45, color = "gray70"),
    legend.position = c(0.94, 0.14),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.key = element_blank()
  )

print(calib_plot_summary)

write_csv(calib_summary_df, "calib_plot_summary_df.csv")
ggsave("calibration_plot_summary.png", width = 8, height = 6, dpi = 600)

rm(list = ls())
