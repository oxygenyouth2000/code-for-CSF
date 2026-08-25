library(openxlsx)
library(tidyverse)
library(cobalt)
library(tibble)

asmd_input <- read.xlsx("df_overall.xlsx", sheet = "df_overall") %>%
  mutate(
    W = ifelse(BB_use %in% c("Yes", 1, "1"), 1, 0),
    weight = as.numeric(weight),
    ITE_group = factor(as.character(ITE_group), levels = c("Lower third", "Middle third", "Upper third"))
  )

asmd_covariates <- c(
  "sex", "age", "BMI", "smoking_status", "drinking_status",
  "hypertension", "diabetes", "stroke", "hypercholesterolemia",
  "prior_CAD", "prior_PAD", "prior_HF", "bradycardia",
  "conduction_disease", "copd", "asthma", "prior_BB", "MI_type",
  "Killip_class", "multivessel_disease", "revascularization", "LVEF",
  "eGFR", "hemoglobin", "LDL", "HbA1c", "sodium", "potassium",
  "asprin", "P2Y12_receptor_blocker", "ACEI_or_ARB", "diuretics",
  "CCB", "MRA", "anticoagulants", "statin", "HR", "SBP", "hs_cTnI",
  "NT_proBNP", "WBC", "PLT", "INR", "APTT"
)

extract_asmd_long <- function(data, tertile_label, covariates) {
  available_covariates <- covariates[covariates %in% names(data)]
  analysis_data <- data %>%
    filter(!is.na(W), !is.na(weight), is.finite(weight), weight > 0)

  if (length(available_covariates) == 0) {
    stop("No eligible covariates available for ASMD calculation.")
  }

  bal_formula <- as.formula(paste("W ~", paste(available_covariates, collapse = " + ")))
  bal_res <- cobalt::bal.tab(
    x = bal_formula,
    data = analysis_data,
    weights = analysis_data$weight,
    method = "weighting",
    estimand = "ATE",
    un = TRUE,
    s.d.denom = "pooled"
  )

  bal_res$Balance %>%
    as.data.frame() %>%
    rownames_to_column("variable") %>%
    transmute(
      variable = variable,
      tertile = tertile_label,
      n_group = nrow(analysis_data),
      Unweighted = abs(Diff.Un),
      `IPTW weighted` = abs(Diff.Adj)
    ) %>%
    pivot_longer(cols = c(Unweighted, `IPTW weighted`), names_to = "Method", values_to = "ASMD") %>%
    mutate(ASMD = round(ASMD, 3)) %>%
    arrange(variable, Method)
}

overall_asmd <- extract_asmd_long(asmd_input, "Overall", asmd_covariates)
lower <- extract_asmd_long(asmd_input %>% filter(ITE_group == "Lower third"), "Lower third", asmd_covariates)
middle <- extract_asmd_long(asmd_input %>% filter(ITE_group == "Middle third"), "Middle third", asmd_covariates)
upper <- extract_asmd_long(asmd_input %>% filter(ITE_group == "Upper third"), "Upper third", asmd_covariates)

all_long <- bind_rows(overall_asmd, lower, middle, upper)
summary_by_group <- all_long %>%
  group_by(tertile, Method) %>%
  summarise(
    n_variables = n(),
    mean_ASMD = round(mean(ASMD, na.rm = TRUE), 3),
    median_ASMD = round(median(ASMD, na.rm = TRUE), 3),
    max_ASMD = round(max(ASMD, na.rm = TRUE), 3),
    prop_below_0_1 = round(mean(ASMD < 0.1, na.rm = TRUE), 3),
    .groups = "drop"
  )

wb <- createWorkbook()
addWorksheet(wb, "overall")
writeData(wb, "overall", overall_asmd)
addWorksheet(wb, "lower_third")
writeData(wb, "lower_third", lower)
addWorksheet(wb, "middle_third")
writeData(wb, "middle_third", middle)
addWorksheet(wb, "upper_third")
writeData(wb, "upper_third", upper)
addWorksheet(wb, "all_tertiles_long")
writeData(wb, "all_tertiles_long", all_long)
addWorksheet(wb, "summary_check")
writeData(wb, "summary_check", summary_by_group)

walk(
  .x = names(wb),
  .f = ~ setColWidths(wb = wb, sheet = .x, cols = 1:10, widths = "auto")
)

saveWorkbook(wb = wb, file = "asmd_tertiles.xlsx", overwrite = TRUE)
