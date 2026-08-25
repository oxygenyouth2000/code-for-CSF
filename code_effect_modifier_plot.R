library(openxlsx)
library(survRM2)
library(dplyr)
library(tibble)
library(ggplot2)

df_overall <- read.xlsx("df_overall.xlsx", sheet = "df_overall") %>%
  mutate(
    mace = ifelse(mace %in% c("Happened", 1, "1"), 1, 0),
    BB_use = ifelse(BB_use %in% c("Yes", 1, "1"), 1, 0),
    ITE_group = factor(as.character(ITE_group), levels = c("Lower third", "Middle third", "Upper third"))
  )

df_low <- df_overall %>% filter(ITE_group == "Lower third")
df_middle <- df_overall %>% filter(ITE_group == "Middle third")
df_upper <- df_overall %>% filter(ITE_group == "Upper third")
tau_rmst <- 36

overall_result <- rmst2(time = df_overall$survival_time, status = df_overall$mace, arm = df_overall$BB_use, tau = tau_rmst)
low_result <- rmst2(time = df_low$survival_time, status = df_low$mace, arm = df_low$BB_use, tau = tau_rmst)
middle_result <- rmst2(time = df_middle$survival_time, status = df_middle$mace, arm = df_middle$BB_use, tau = tau_rmst)
upper_result <- rmst2(time = df_upper$survival_time, status = df_upper$mace, arm = df_upper$BB_use, tau = tau_rmst)

extract_rmst_diff <- function(rmst_obj, group_label) {
  rmst_vector <- rmst_obj$unadjusted.result["RMST (arm=1)-(arm=0)", c("Est.", "lower .95", "upper .95")]
  tibble(
    ITE_group = group_label,
    est = as.numeric(rmst_vector[["Est."]]),
    LCI = as.numeric(rmst_vector[["lower .95"]]),
    UCI = as.numeric(rmst_vector[["upper .95"]])
  )
}

compute_hc3_vcov <- function(lm_fit) {
  model_matrix <- model.matrix(lm_fit)
  fitted_residuals <- residuals(lm_fit)
  leverage <- hatvalues(lm_fit)
  leverage_adjusted <- pmax(1 - leverage, 1e-8)
  bread_matrix <- solve(crossprod(model_matrix))
  meat_matrix <- t(model_matrix) %*% diag((fitted_residuals / leverage_adjusted)^2) %*% model_matrix
  bread_matrix %*% meat_matrix %*% bread_matrix
}

compute_joint_wald_p <- function(beta_vector, vcov_matrix) {
  wald_statistic <- drop(t(beta_vector) %*% solve(vcov_matrix, beta_vector))
  pchisq(wald_statistic, df = length(beta_vector), lower.tail = FALSE)
}

compute_rmst_interaction <- function(data, tau_value = 36) {
  analysis_data <- data %>%
    mutate(ITE_group = factor(ITE_group, levels = c("Lower third", "Middle third", "Upper third")))

  rmst_fit <- survival::survfit(
    survival::Surv(survival_time, mace) ~ 1,
    data = analysis_data,
    model = TRUE
  )

  analysis_data <- analysis_data %>%
    mutate(
      rmst_pseudo = as.numeric(
        survival::pseudo(fit = rmst_fit, times = tau_value, type = "sojourn")
      )
    )

  interaction_fit <- lm(rmst_pseudo ~ BB_use * ITE_group, data = analysis_data)
  interaction_vcov <- compute_hc3_vcov(interaction_fit)
  interaction_terms <- grep("^BB_use:ITE_group", names(coef(interaction_fit)), value = TRUE)
  interaction_p_value <- compute_joint_wald_p(
    beta_vector = coef(interaction_fit)[interaction_terms],
    vcov_matrix = interaction_vcov[interaction_terms, interaction_terms, drop = FALSE]
  )

  list(interaction_p_value = interaction_p_value)
}

format_p_value <- function(p_value) {
  ifelse(is.na(p_value), "NA", ifelse(p_value < 0.001, "<0.001", paste0("= ", formatC(p_value, format = "f", digits = 3))))
}

RMST_result <- bind_rows(
  extract_rmst_diff(overall_result, "Overall"),
  extract_rmst_diff(low_result, "Lower third"),
  extract_rmst_diff(middle_result, "Middle third"),
  extract_rmst_diff(upper_result, "Upper third")
) %>%
  mutate(ITE_group = factor(ITE_group, levels = c("Overall", "Lower third", "Middle third", "Upper third")))

save(RMST_result, file = "中间分析结果数据/RMST_result.RData")
write.csv(RMST_result, "中间分析结果数据/RMST_result.csv", row.names = FALSE)

interaction_result <- compute_rmst_interaction(data = df_overall, tau_value = tau_rmst)
interaction_summary <- tibble(tau = tau_rmst, interaction_p_value = interaction_result$interaction_p_value)
write.csv(interaction_summary, "中间分析结果数据/RMST_interaction_p.csv", row.names = FALSE)

interaction_label <- paste0("Interaction P", format_p_value(interaction_result$interaction_p_value))

ite_forest_plot <- ggplot(RMST_result, aes(x = ITE_group, y = est)) +
  geom_errorbar(aes(ymin = LCI, ymax = UCI, color = ITE_group), width = 0.2, linewidth = 1) +
  geom_point(aes(color = ITE_group), size = 4, shape = 19) +
  scale_color_manual(values = c("Overall" = "black", "Lower third" = "steelblue", "Middle third" = "grey50", "Upper third" = "red")) +
  annotate("text", x = 3.2, y = max(RMST_result$UCI, na.rm = TRUE) * 0.9, label = "Favors benefit from BB", fontface = "italic", size = 6, color = "red") +
  annotate("text", x = 2.9, y = min(RMST_result$LCI, na.rm = TRUE) * 0.9, label = "Favors no benefit from BB", fontface = "italic", size = 6, color = "steelblue") +
  annotate("text", x = 3, y = max(RMST_result$UCI, na.rm = TRUE) * 1.15, label = interaction_label, size = 6, fontface = "italic") +
  labs(
    x = "Patients grouped by predicted individualized treatment effect",
    y = "Restricted mean survival time difference (month)",
    color = NULL
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(-3, 2), breaks = seq(-3, 2, by = 1)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.5, color = "gray80"),
    panel.grid.minor.y = element_line(linewidth = 0.5, color = "gray70"),
    legend.position = "none"
  )

print(ite_forest_plot)
ggsave("RMST_Forest.png", ite_forest_plot, width = 8, height = 6, dpi = 600)
