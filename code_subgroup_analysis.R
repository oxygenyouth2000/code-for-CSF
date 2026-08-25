library(openxlsx)
library(dplyr)
library(ggplot2)
library(rlang)

df_imputed <- read.xlsx("分析结果数据/df_imputed.xlsx")
ite_values <- read.csv("中间分析结果数据/ite_values.csv")

if ("id" %in% names(df_imputed) && "id" %in% names(ite_values)) {
  df_imputed <- df_imputed %>%
    left_join(ite_values %>% select(id, ITE), by = "id")
}

df_imputed <- df_imputed %>%
  mutate(
    LVEF_group = cut(lvef, breaks = c(-Inf, 50, Inf), labels = c("LVEF 40-50%", "LVEF >50%")),
    HR_group = cut(HR, breaks = quantile(HR, probs = c(0, 0.33, 0.67, 1), na.rm = TRUE), labels = c("low HR", "middle HR", "high HR"))
  )

subgroup_analysis <- function(data, group_var) {
  results <- data %>%
    group_by(!!sym(group_var)) %>%
    summarise(
      n = n(),
      mean_ITE = mean(ITE, na.rm = TRUE),
      se_ITE = sd(ITE, na.rm = TRUE) / sqrt(n()),
      ci_lower = mean_ITE - 1.96 * se_ITE,
      ci_upper = mean_ITE + 1.96 * se_ITE,
      .groups = "drop"
    )

  plot_obj <- ggplot(results, aes(x = !!sym(group_var), y = mean_ITE)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(x = group_var, y = "mean ITE", title = paste("Treatment effect stratified by", group_var)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  list(table = results, plot = plot_obj)
}

lvef_subgroup <- subgroup_analysis(df_imputed, "LVEF_group")
hr_subgroup <- subgroup_analysis(df_imputed, "HR_group")

print(lvef_subgroup$plot)
print(hr_subgroup$plot)

write.csv(lvef_subgroup$table, "中间分析结果数据/lvef_subgroup.csv", row.names = FALSE)
write.csv(hr_subgroup$table, "中间分析结果数据/hr_subgroup.csv", row.names = FALSE)

ggsave("LVEF_Subgroup.png", lvef_subgroup$plot, width = 8, height = 6, dpi = 600)
ggsave("HR_Subgroup.png", hr_subgroup$plot, width = 8, height = 6, dpi = 600)
