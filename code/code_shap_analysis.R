rm(list = ls())

library(openxlsx)
library(kernelshap)
library(shapviz)
library(tidyverse)
library(tidyr)
library(caret)

pfun <- function(object, newdata) {
  as.numeric(predict(object, as.matrix(newdata))$predictions)
}

csf_fit <- readRDS("csf_fit.rds")
df_imputed <- read.xlsx("df_imputed.xlsx")

x_vars <- c(
  "sex", "age", "BMI", "smoking_status", "drinking_status",
  "hypertension", "diabetes", "stroke", "hypercholesterolemia",
  "prior_CAD", "prior_PAD", "prior_HF", "bradycardia",
  "conduction_disease", "copd", "asthma", "prior_BB", "MI_type",
  "Killip_class", "multivessel_disease", "revascularization", "LVEF",
  "eGFR", "hemoglobin", "LDL", "HbA1c", "sodium", "potassium",
  "asprin", "P2Y12_receptor_blocker", "ACEI_or_ARB", "diuretics",
  "CCB", "MRA", "anticoagulants", "statin", "HR", "SBP", "cTnI",
  "NT_proBNP", "WBC", "PLT", "INR", "APTT"
)

X <- df_imputed[, x_vars]

# 分类变量哑变量展开
dummy_model <- dummyVars(" ~ .", data = X, fullRank = TRUE)

# 连续变量标准化
X_preprocessed <- predict(dummy_model, newdata = X)
X_scaled <- scale(X_preprocessed)

set.seed(2026)

bg_idx <- sample(nrow(X_scaled), min(1000, nrow(X_scaled)))
bg_X <- X_scaled[bg_idx, ]
explain_idx <- sample(nrow(X_scaled), min(5000, nrow(X_scaled)))
explain_X <- X_scaled[explain_idx, ]
explain_X_orig <- X_preprocessed[explain_idx, ]

# kernelshap + shapviz
ks <- kernelshap(
  object = csf_fit,
  X = explain_X,
  bg_X = bg_X,
  pred_fun = pfun,
  seed = 2026
)
save(ks, file = "code_all/ks.RData")

shp <- shapviz(ks, X = explain_X_orig)

if (!dir.exists("paper figures and tables")) {
  dir.create("paper figures and tables")
}

shap_bar <- sv_importance(shp, kind = "bar", show_numbers = TRUE, max_display = 15) +
  theme_minimal(base_size = 14) +
  ggtitle("SHAP Feature Importance (Bar Plot)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(shap_bar)
ggsave("paper figures and tables/SHAP_Bar_Plot.png", shap_bar, dpi = 600, width = 8, height = 6)

shap_beeswarm <- sv_beeswarm(shp, max_display = 15) +
  theme_minimal(base_size = 14) +
  ggtitle("SHAP Summary (Beeswarm Plot)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(shap_beeswarm)
ggsave("paper figures and tables/SHAP_Beeswarm_Plot.png", shap_beeswarm, dpi = 600, width = 8, height = 6)

shap_df <- as.data.frame(get_shap_values(shp))
shap_long <- shap_df %>%
  pivot_longer(cols = everything(), names_to = "Feature", values_to = "SHAP_value") %>%
  group_by(Feature) %>%
  mutate(mean_abs_shap = mean(abs(SHAP_value))) %>%
  ungroup() %>%
  arrange(mean_abs_shap) %>%
  mutate(Feature = factor(Feature, levels = unique(Feature)))

shap_violin <- ggplot(shap_long, aes(x = Feature, y = SHAP_value)) +
  geom_violin(aes(fill = Feature), show.legend = FALSE, alpha = 0.6, color = "gray50", scale = "width") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(title = "SHAP Violin Plot", x = "", y = "SHAP Value (Impact on ITE)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), panel.grid.minor = element_blank())
print(shap_violin)
ggsave("paper figures and tables/SHAP_Violin_Plot.png", shap_violin, dpi = 600, width = 8, height = 6)

top_feature <- rev(levels(shap_long$Feature))[3]
p_scatter <- sv_dependence(shp, v = top_feature, color_var = "auto", alpha = 0.6, size = 2) +
  theme_minimal(base_size = 14) +
  ggtitle(paste0("SHAP Scatter/Dependence Plot (", top_feature, ")")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(p_scatter)
ggsave(paste0("paper figures and tables/SHAP_Scatter_", top_feature, ".png"), p_scatter, dpi = 600, width = 7, height = 5)

p_waterfall <- sv_waterfall(shp, row_id = 1) +
  theme_minimal(base_size = 14) +
  ggtitle("SHAP Waterfall Plot (Patient 1)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(p_waterfall)
ggsave("paper figures and tables/SHAP_Waterfall_Plot.png", p_waterfall, dpi = 600, width = 8, height = 6)

p_force <- sv_force(shp, row_id = 1) +
  theme_minimal(base_size = 14) +
  ggtitle("SHAP Force Plot (Patient 1)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(p_force)
ggsave("paper figures and tables/SHAP_Force_Plot.png", p_force, dpi = 600, width = 10, height = 4)

top10_features <- shap_long %>%
  group_by(Feature) %>%
  summarise(mean_abs_shap = mean(abs(SHAP_value), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abs_shap)) %>%
  slice_head(n = 10) %>%
  mutate(
    Scaled_Importance = round(mean_abs_shap / max(mean_abs_shap) * 100, 1),
    Feature = factor(Feature, levels = rev(Feature))
  )

varimp_plot <- ggplot(top10_features, aes(x = Scaled_Importance, y = Feature)) +
  geom_col(fill = "steelblue", width = 0.7) +
  scale_x_continuous(limits = c(0, 105), expand = c(0, 0)) +
  labs(x = "Relative importance (%)", y = "Patient characteristic") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    panel.grid = element_blank()
  )

print(varimp_plot)
write.csv(top10_features, "中间分析结果数据/top10_feature_importance.csv", row.names = FALSE)
ggsave("TOP10_Importance.png", width = 8, height = 6, dpi = 600)
