library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggpubr)

df_overall <- read.xlsx("df_overall.xlsx", sheet = 1)

pdp_lvef <- ggplot(df_overall, aes(x = lvef, y = ITE)) +
  geom_point(color = "#9EC9E2", size = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "steelblue", fill = "#9EC9E2", linewidth = 1) +
  labs(x = "Left ventricular ejection fraction (%)", y = "Mean estimated individual treatment effect") +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 20), axis.text = element_text(size = 20), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_line(linewidth = 0.5, color = "gray80"), panel.grid.minor.y = element_line(linewidth = 0.5, color = "gray70"), axis.line.y.left = element_line(color = "black", linewidth = 0.5), axis.line.x.bottom = element_line(color = "black", linewidth = 0.5))

print(pdp_lvef)
ggsave(file.path(output_dir, "SHAP_pdp_lvef.png"), plot = pdp_lvef, width = 8, height = 6, dpi = 600)

pdp_age <- ggplot(df_overall, aes(x = age, y = ITE)) +
  geom_point(color = "#9EC9E2", size = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "steelblue", fill = "#9EC9E2", linewidth = 1) +
  labs(x = "Age (years)", y = "Mean estimated individual treatment effect") +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 20), axis.text = element_text(size = 20), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_line(linewidth = 0.5, color = "gray80"), panel.grid.minor.y = element_line(linewidth = 0.5, color = "gray70"), axis.line.y.left = element_line(color = "black", linewidth = 0.5), axis.line.x.bottom = element_line(color = "black", linewidth = 0.5))

ggsave(file.path(output_dir, "SHAP_pdp_age.png"), plot = pdp_age, width = 8, height = 6, dpi = 600)

pdp_heartrate <- ggplot(df_overall, aes(x = HR, y = ITE)) +
  geom_point(color = "#9EC9E2", size = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "steelblue", fill = "#9EC9E2", linewidth = 1) +
  labs(x = "Heart rate (bpm)", y = "Mean estimated individual treatment effect") +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 20), axis.text = element_text(size = 20), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_line(linewidth = 0.5, color = "gray80"), panel.grid.minor.y = element_line(linewidth = 0.5, color = "gray70"), axis.line.y.left = element_line(color = "black", linewidth = 0.5), axis.line.x.bottom = element_line(color = "black", linewidth = 0.5))

ggsave(file.path(output_dir, "SHAP_pdp_heartrate.png"), plot = pdp_heartrate, width = 8, height = 6, dpi = 600)

pdp_egfr <- ggplot(df_overall, aes(x = egfr, y = ITE)) +
  geom_point(color = "#9EC9E2", size = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "steelblue", fill = "#9EC9E2", linewidth = 1) +
  labs(x = "Estimated glomerular filtration rate (mL/min)", y = "Mean estimated individual treatment effect") +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 16), axis.text = element_text(size = 16), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_line(linewidth = 0.5, color = "gray80"), panel.grid.minor.y = element_line(linewidth = 0.5, color = "gray70"), axis.line.y.left = element_line(color = "black", linewidth = 0.5), axis.line.x.bottom = element_line(color = "black", linewidth = 0.5))

ggsave(file.path(output_dir, "SHAP_pdp_egfr.png"), plot = pdp_egfr, width = 8, height = 6, dpi = 600)
