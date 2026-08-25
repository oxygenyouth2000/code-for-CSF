library(openxlsx)
library(ggplot2)
library(dplyr)
library(forcats)
library(tidyverse)

df_overall <- read.xlsx("df_overall.xlsx")

# 将Killip分级整理为有序四分类，便于绘制ITE分布图
df_overall <- df_overall %>%
  mutate(
    killip_class = factor(
      as.character(killip_class),
      levels = c("1", "2", "3", "4"),
      labels = c("Class I", "Class II", "Class III", "Class IV")
    )
  )

plot_distribution <- function(data, fill_var, fill_label, output_name) {
  plot_obj <- ggplot(data, aes(x = ITE, fill = .data[[fill_var]])) +
    geom_histogram(bins = 50, color = "black", position = "stack", alpha = 1) +
    geom_vline(xintercept = quantile(data$ITE, probs = c(1 / 3, 2 / 3), na.rm = TRUE), color = "gray20", linetype = "dashed", linewidth = 1) +
    labs(
      x = "Estimated individualized treatment effect",
      y = "No. of Patients",
      fill = fill_label
    ) +
    scale_x_continuous(limits = c(-7, 7), breaks = seq(-7, 7, by = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.text = element_text(size = 16),
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.position = c(1, 0.05),
      legend.justification = c("right", "bottom"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.5, color = "gray80"),
      panel.grid.minor.y = element_line(linewidth = 0.5, color = "gray70"),
      legend.title = element_text(size = 16),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5)
    )

  print(plot_obj)
  ggsave(output_name, plot_obj, dpi = 600, width = 8, height = 6)
}

plot_distribution(df_overall, "sex", "sex", "Sex_Distribution.png")
plot_distribution(df_overall, "mi_type", "MI type", "MI_Distribution.png")
plot_distribution(df_overall, "hypertension", "hypertension", "HTN_Distribution.png")
plot_distribution(df_overall, "revascularization", "Revascularization", "Revascularization_Distribution.png")
plot_distribution(df_overall, "multi_vessel", "multivessel disease", "MVD_Distribution.png")
plot_distribution(df_overall, "smoking", "Current smoker", "SMOK_Distribution.png")
plot_distribution(df_overall, "diabetes", "Diabetes", "Diabetes_Distribution.png")
plot_distribution(df_overall, "killip_class", "Killip class", "Killip_Distribution.png")
