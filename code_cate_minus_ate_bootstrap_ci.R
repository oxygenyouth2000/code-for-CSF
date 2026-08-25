library(tidyverse)
library(caret)
library(grf)
library(survival)
library(readr)

# Bootstrap 95%CI for CATE minus ATE across tertiles
# 计算流程与主分析保持一致：
# 1. 从总体样本中进行有放回抽样，生成bootstrap样本；
# 2. 使用既有CSF模型对每个bootstrap样本中的全部个体重新预测ITE；
# 3. 在每个bootstrap样本内部按预测ITE重新分成三个tertiles；
# 4. 计算该bootstrap样本的Overall ATE与3个tertile的CATE；
# 5. 计算 CATE - ATE；
# 6. 采用经验分位数法（2.5%和97.5%分位数）构建95%CI。

reference_data_path <- "中间分析结果数据/df_overall.csv"
csf_model_path <- "中间分析结果数据/csf_fit_new.rds"
bootstrap_draws_output_path <- "中间分析结果数据/cate_minus_ate_bootstrap_draws.csv"
bootstrap_ci_output_path <- "中间分析结果数据/cate_minus_ate_bootstrap_ci.csv"

# 主分析已得到的点估计
# 这里固定写入主分析报告的3个点估计；bootstrap部分只负责构建其95%CI。
cate_minus_ate_point_est <- tibble(
  ITE_group = factor(
    c("Lower third", "Middle third", "Upper third"),
    levels = c("Lower third", "Middle third", "Upper third")
  ),
  est = c(-3.11, -0.33, 2.50)
)

# 分析参数
tau_rmst <- 36
n_boot_cate_minus_ate <- 99
set.seed(2026)

# CSF模型训练时使用的44个候选预测变量
# 为保证bootstrap预测输入与原模型完全一致，这里的变量顺序必须与主分析保持一致。
x_vars_boot <- c(
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

# 读取主分析总体数据，并将关键变量转成后续函数需要的编码
df_boot_reference <- read_csv(reference_data_path, show_col_types = FALSE) %>%
  mutate(
    BB_use = as.numeric(BB_use == "Yes"),
    mace = as.numeric(mace == "Happened")
  )

# 读取已拟合好的CSF模型
csf_fit_boot <- readRDS(csf_model_path)

# 重新构造与原模型训练阶段一致的 dummy 编码与标准化参数
# 这里非常关键：
# bootstrap样本中的新数据不能重新估计dummy规则或缩放参数，
# 而必须复用“主分析全样本”的编码规则和标准化参数。
X_reference_raw <- df_boot_reference[, x_vars_boot]
dummy_model_boot <- caret::dummyVars(" ~ .", data = X_reference_raw, fullRank = TRUE)
X_reference_preprocessed <- predict(dummy_model_boot, newdata = X_reference_raw)
X_reference_scaled <- scale(X_reference_preprocessed)
scale_center_boot <- attr(X_reference_scaled, "scaled:center")
scale_scale_boot <- attr(X_reference_scaled, "scaled:scale")

# 为避免极少数列标准差为0导致除以0，这里将无效标准差替换为1。
scale_scale_boot[!is.finite(scale_scale_boot) | scale_scale_boot == 0] <- 1

# 为任意新样本生成与原CSF训练阶段完全一致的输入矩阵
preprocess_csf_newdata <- function(newdata,
                                   feature_columns,
                                   dummy_model,
                                   scale_center,
                                   scale_scale) {
  # 先按原始44个候选变量取子集。
  raw_x <- newdata[, feature_columns, drop = FALSE]

  # 使用主分析中拟定的dummy规则展开分类变量。
  encoded_x <- predict(dummy_model, newdata = raw_x)

  # 使用主分析全样本的中心与标准差完成标准化。
  scaled_x <- sweep(encoded_x, 2, scale_center, FUN = "-")
  scaled_x <- sweep(scaled_x, 2, scale_scale, FUN = "/")
  storage.mode(scaled_x) <- "double"

  scaled_x
}

# 在当前分析样本内重新估计稳定化IPTW
estimate_stabilized_iptw <- function(data, covariates) {
  available_covariates <- covariates[covariates %in% names(data)]

  if (length(available_covariates) == 0) {
    stop("当前样本中没有可用于PS建模的协变量。")
  }

  analysis_data <- data %>%
    mutate(
      BB_use = as.numeric(BB_use),
      mace = as.numeric(mace)
    ) %>%
    filter(
      !is.na(BB_use),
      !is.na(survival_time),
      !is.na(mace)
    ) %>%
    filter(
      if_all(all_of(available_covariates), ~ !is.na(.x))
    )

  if (nrow(analysis_data) == 0) {
    stop("当前样本在PS建模所需变量上全部缺失，无法估计IPTW。")
  }

  if (dplyr::n_distinct(analysis_data$BB_use) < 2) {
    stop("当前样本仅包含单一治疗组，无法估计IPTW。")
  }

  model_data <- analysis_data %>%
    select(BB_use, all_of(available_covariates)) %>%
    mutate(
      across(where(is.character), as.factor)
    )

  ps_formula <- as.formula(
    paste("BB_use ~", paste(available_covariates, collapse = " + "))
  )

  ps_fit <- glm(ps_formula, data = model_data, family = binomial())
  ps <- predict(ps_fit, type = "response") %>%
    pmin(0.99) %>%
    pmax(0.01)

  p_treat <- mean(model_data$BB_use == 1)
  weight <- ifelse(
    model_data$BB_use == 1,
    p_treat / ps,
    (1 - p_treat) / (1 - ps)
  )

  q_lo <- quantile(weight, probs = 0.01, na.rm = TRUE)
  q_hi <- quantile(weight, probs = 0.99, na.rm = TRUE)
  weight <- pmin(pmax(weight, q_lo), q_hi)

  analysis_data %>%
    mutate(
      ps = as.numeric(ps),
      weight = as.numeric(weight)
    )
}

# 对加权Kaplan-Meier曲线在0到tau区间内做阶梯积分，得到RMST
compute_weighted_km_rmst <- function(fit_object, tau = 36) {
  strata_sizes <- fit_object$strata

  if (is.null(strata_sizes)) {
    strata_sizes <- length(fit_object$time)
    names(strata_sizes) <- "All"
  }

  strata_names <- names(strata_sizes)
  strata_index <- rep(strata_names, strata_sizes)
  split_time <- split(fit_object$time, strata_index)
  split_surv <- split(fit_object$surv, strata_index)

  bind_rows(lapply(strata_names, function(strata_name) {
    time_vector <- split_time[[strata_name]]
    surv_vector <- split_surv[[strata_name]]
    valid_index <- time_vector < tau

    integration_time <- c(0, time_vector[valid_index], tau)
    integration_surv <- c(1, surv_vector[valid_index])

    tibble(
      strata = strata_name,
      rmst = sum(diff(integration_time) * integration_surv)
    )
  }))
}

# 计算单个治疗组在给定tau下的加权RMST及其IJ标准误
compute_weighted_group_rmst <- function(data, arm_value, tau = 36) {
  subgroup_data <- data %>%
    filter(
      BB_use == arm_value,
      !is.na(survival_time),
      !is.na(mace),
      !is.na(weight),
      is.finite(weight),
      weight > 0
    )

  subgroup_fit <- survival::survfit(
    survival::Surv(survival_time, mace) ~ 1,
    data = subgroup_data,
    weights = weight,
    robust = TRUE,
    se.fit = TRUE,
    conf.type = "plain",
    model = TRUE
  )

  rmst_estimate <- compute_weighted_km_rmst(subgroup_fit, tau = tau) %>%
    pull(rmst)

  rmst_ij <- residuals(
    subgroup_fit,
    times = tau,
    type = "rmst",
    weighted = TRUE
  )

  tibble(
    arm = arm_value,
    rmst = rmst_estimate,
    rmst_se = sqrt(sum(as.numeric(rmst_ij)^2, na.rm = TRUE))
  )
}

# 计算某个样本（Overall或某个tertile）的加权RMST difference及95%CI
compute_weighted_rmst_difference <- function(data,
                                             group_label,
                                             tau = 36,
                                             ps_covariates = x_vars_boot) {
  analysis_data <- estimate_stabilized_iptw(data, covariates = ps_covariates) %>%
    filter(
      !is.na(survival_time),
      !is.na(mace),
      !is.na(BB_use),
      !is.na(weight),
      is.finite(weight),
      weight > 0
    )

  rmst_by_arm <- bind_rows(
    compute_weighted_group_rmst(analysis_data, arm_value = 1, tau = tau),
    compute_weighted_group_rmst(analysis_data, arm_value = 0, tau = tau)
  )

  difference_estimate <- rmst_by_arm$rmst[rmst_by_arm$arm == 1] -
    rmst_by_arm$rmst[rmst_by_arm$arm == 0]
  difference_se <- sqrt(sum(rmst_by_arm$rmst_se^2))
  z_value <- qnorm(0.975)

  tibble(
    ITE_group = group_label,
    est = as.numeric(difference_estimate),
    LCI = as.numeric(difference_estimate - z_value * difference_se),
    UCI = as.numeric(difference_estimate + z_value * difference_se)
  )
}

# 在单个bootstrap样本内完成：
# 抽样 -> 预测ITE -> 重新分tertiles -> 计算CATE-ATE
compute_bootstrap_cate_minus_ate <- function(iter_id,
                                             analysis_data,
                                             csf_fit,
                                             feature_columns,
                                             dummy_model,
                                             scale_center,
                                             scale_scale,
                                             tau_value = 36) {
  # 有放回抽样，样本量与原始总体一致。
  bootstrap_index <- sample(
    x = seq_len(nrow(analysis_data)),
    size = nrow(analysis_data),
    replace = TRUE
  )

  bootstrap_data <- analysis_data[bootstrap_index, , drop = FALSE]

  # 使用既有CSF模型对bootstrap样本逐例预测ITE。
  bootstrap_x_scaled <- preprocess_csf_newdata(
    newdata = bootstrap_data,
    feature_columns = feature_columns,
    dummy_model = dummy_model,
    scale_center = scale_center,
    scale_scale = scale_scale
  )

  bootstrap_data <- bootstrap_data %>%
    mutate(
      ITE_boot = as.numeric(
        predict(csf_fit, newdata = bootstrap_x_scaled)$predictions
      )
    )

  # 在bootstrap样本内部，按预测ITE重新确定三分位切点。
  ite_cut_1_boot <- quantile(bootstrap_data$ITE_boot, probs = 1 / 3, na.rm = TRUE)
  ite_cut_2_boot <- quantile(bootstrap_data$ITE_boot, probs = 2 / 3, na.rm = TRUE)

  bootstrap_data <- bootstrap_data %>%
    mutate(
      ITE_group = factor(
        case_when(
          ITE_boot < ite_cut_1_boot ~ "Lower third",
          ITE_boot > ite_cut_2_boot ~ "Upper third",
          TRUE ~ "Middle third"
        ),
        levels = c("Lower third", "Middle third", "Upper third")
      )
    )

  # 先在该bootstrap样本内计算Overall ATE。
  ate_boot <- compute_weighted_rmst_difference(
    data = bootstrap_data,
    group_label = "Overall",
    tau = tau_value,
    ps_covariates = feature_columns
  )

  # 再分别计算3个tertile的CATE，并减去该样本对应的ATE。
  cate_boot <- bind_rows(
    compute_weighted_rmst_difference(
      data = bootstrap_data %>% filter(ITE_group == "Lower third"),
      group_label = "Lower third",
      tau = tau_value,
      ps_covariates = feature_columns
    ),
    compute_weighted_rmst_difference(
      data = bootstrap_data %>% filter(ITE_group == "Middle third"),
      group_label = "Middle third",
      tau = tau_value,
      ps_covariates = feature_columns
    ),
    compute_weighted_rmst_difference(
      data = bootstrap_data %>% filter(ITE_group == "Upper third"),
      group_label = "Upper third",
      tau = tau_value,
      ps_covariates = feature_columns
    )
  ) %>%
    mutate(
      bootstrap_id = iter_id,
      ate_est = ate_boot$est[[1]],
      cate_minus_ate = est - ate_est
    ) %>%
    select(bootstrap_id, ITE_group, ate_est, cate_est = est, cate_minus_ate)

  cate_boot
}
-
# 正式运行99次bootstrap
bootstrap_cate_minus_ate_draws <- map_dfr(
  .x = seq_len(n_boot_cate_minus_ate),
  .f = function(iter_id) {
    message("Running bootstrap iteration: ", iter_id, " / ", n_boot_cate_minus_ate)

    compute_bootstrap_cate_minus_ate(
      iter_id = iter_id,
      analysis_data = df_boot_reference,
      csf_fit = csf_fit_boot,
      feature_columns = x_vars_boot,
      dummy_model = dummy_model_boot,
      scale_center = scale_center_boot,
      scale_scale = scale_scale_boot,
      tau_value = tau_rmst
    )
  }
)

# 采用经验分位数法（percentile bootstrap）构建95%CI
cate_minus_ate_ci <- bootstrap_cate_minus_ate_draws %>%
  group_by(ITE_group) %>%
  summarise(
    n_boot = sum(!is.na(cate_minus_ate)),
    LCI = quantile(cate_minus_ate, probs = 0.025, na.rm = TRUE),
    UCI = quantile(cate_minus_ate, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )

cate_minus_ate_result <- cate_minus_ate_point_est %>%
  left_join(cate_minus_ate_ci, by = "ITE_group")


# 输出结果
print(bootstrap_cate_minus_ate_draws)
print(cate_minus_ate_result)

write_csv(bootstrap_cate_minus_ate_draws, bootstrap_draws_output_path)
write_csv(cate_minus_ate_result, bootstrap_ci_output_path)
