library(openxlsx)
library(grf)
library(caret)
library(dplyr)

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

feature_columns <- x_vars
X <- df_imputed[, feature_columns]
dummy_model <- dummyVars(" ~ .", data = X, fullRank = TRUE)
X_preprocessed <- predict(dummy_model, newdata = X)
X_scaled <- scale(X_preprocessed)

W <- as.numeric(df_imputed[["BB_use"]] == "Yes")
Y <- df_imputed[["survival_time"]]
D <- ifelse(df_imputed[["mace"]] == "Happened", 1, 0)

csf_fit <- causal_survival_forest(
  X = X_scaled,
  Y = Y,
  W = W,
  D = D,
  target = "RMST",
  horizon = 36,
  honesty = TRUE,
  compute.oob.predictions = TRUE,
  seed = 2026
)

ites <- predict(csf_fit, X_scaled, estimate.variance = TRUE)
df_imputed$ITE <- ites$predictions
df_imputed$ITE_se <- sqrt(ites$variance.estimates)

if (!dir.exists("中间分析结果数据")) {
  dir.create("中间分析结果数据", recursive = TRUE)
}

ite_values <- df_imputed %>%
  mutate(
    ITE_group = case_when(
      ITE < quantile(ITE, 1 / 3, na.rm = TRUE) ~ "Lower third",
      ITE > quantile(ITE, 2 / 3, na.rm = TRUE) ~ "Upper third",
      TRUE ~ "Middle third"
    ),
    ITE_group = factor(
      ITE_group,
      levels = c("Lower third", "Middle third", "Upper third")
    ),
    ITE_rank = rank(ITE, ties.method = "first")
  )

write.csv(ite_values, "中间分析结果数据/ite_values.csv", row.names = FALSE)
write.csv(ite_values, "中间分析结果数据/df_overall.csv", row.names = FALSE)
write.xlsx(ite_values, "df_overall.xlsx", sheetName = "df_overall", overwrite = TRUE)
saveRDS(csf_fit, "中间分析结果数据/csf_fit.rds")
