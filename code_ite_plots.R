library(ggplot2)
library(dplyr)
library(forcats)
library(tidyverse)

ite_values <- read.csv("中间分析结果数据/ite_values.csv")

ite_summary <- ite_values %>%
  summarise(
    ate = mean(ITE, na.rm = TRUE),
    sd_val = sd(ITE, na.rm = TRUE),
    valid_n = sum(!is.na(ITE)),
    se = sd_val / sqrt(valid_n),
    margin = qt(0.975, df = valid_n - 1) * se,
    ate_lci = ate - margin,
    ate_uci = ate + margin
  )

benefit_colors <- c(
  "Lower third" = "steelblue",
  "Middle third" = "gray50",
  "Upper third" = "red"
)

ite_dist_plot <- ggplot(ite_values, aes(x = ITE, fill = ITE_group)) +
  geom_histogram(bins = 50, color = "black", position = "stack", alpha = 1) +
  scale_fill_manual(values = benefit_colors) +
  geom_vline(xintercept = quantile(ite_values$ITE, probs = c(1 / 3), na.rm = TRUE), color = "gray20", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = quantile(ite_values$ITE, probs = c(1 / 2), na.rm = TRUE), color = "darkgreen", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = quantile(ite_values$ITE, probs = c(2 / 3), na.rm = TRUE), color = "gray20", linetype = "dashed", linewidth = 1) +
  annotate(
    "text",
    x = ite_summary$ate,
    y = Inf,
    label = sprintf(
      "Average Treatment Effect:\n%.3f (95%% CI: %.3f, %.3f)",
      ite_summary$ate,
      ite_summary$ate_lci,
      ite_summary$ate_uci
    ),
    vjust = 1.1,
    hjust = -0.2,
    color = "darkgreen",
    fontface = "italic",
    size = 6
  ) +
  labs(
    x = "Predicted individualized treatment effect (month)",
    y = "No. of Patients",
    fill = NULL
  ) +
  scale_x_continuous(limits = c(-3, 4), breaks = seq(-3, 4, by = 1)) +
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
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5)
  )

print(ite_dist_plot)
ggsave("ITE_Distribution.png", ite_dist_plot, dpi = 600, width = 8, height = 6)

scatter_summary <- ite_values %>%
  group_by(ITE_group) %>%
  summarise(
    mean_rank = mean(ITE_rank, na.rm = TRUE),
    mean = mean(ITE, na.rm = TRUE),
    sd_val = sd(ITE, na.rm = TRUE),
    valid_n = sum(!is.na(ITE)),
    se = sd_val / sqrt(valid_n),
    margin = qt(0.975, df = valid_n - 1) * se,
    LCI = mean - margin,
    UCI = mean + margin,
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("%.2f (%.2f, %.2f)", mean, LCI, UCI))

ite_scatter_plot <- ggplot(ite_values, aes(x = ITE_rank, y = ITE, color = ITE_group)) +
  geom_point(size = 3.2, stroke = 0) +
  scale_color_manual(values = benefit_colors) +
  annotate("text", x = Inf, y = Inf, label = "Favors benefit from BB", hjust = 1.3, vjust = 2.8, color = "red", fontface = "italic", size = 6) +
  annotate("text", x = -Inf, y = -Inf, label = "Favors no benefit from BB", hjust = -0.3, vjust = -1.1, color = "steelblue", fontface = "italic", size = 6) +
  geom_text(
    data = scatter_summary,
    aes(x = mean_rank, y = c(0.45, 1, 2.5), label = label, color = ITE_group),
    size = 6,
    fontface = "bold",
    show.legend = FALSE
  ) +
  labs(
    x = "Patients ranked by predicted individualized treatment effect",
    y = "Predicted individualized treatment effect (month)"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(-3, 4), breaks = seq(-3, 4, by = 2)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.5, color = "gray80"),
    panel.grid.minor.y = element_line(linewidth = 0.5, color = "gray70"),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5)
  )

print(ite_scatter_plot)
ggsave("ITE_Scatter.png", ite_scatter_plot, dpi = 600, width = 8, height = 6)
