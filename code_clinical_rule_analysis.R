rm(list = ls())

# ------------------------------ #
# Clinical rule analysis based on CART
# ------------------------------ #
# 本脚本基于当前研究已经生成的个体化治疗效应（ITE），
# 参考 Komura et al. (2025) 的 two-step pragmatic subgroup discovery 框架：
# 1. 第一步已由现有 CSF 模型完成，即估计每位患者的 ITE；
# 2. 第二步在高度可解释的临床变量上进行重编码，并用 CART 对 ITE 进行分组发现；
# 3. 之后在每个临床亚组或 CART 终端节点内，按照当前设定的权重策略与 RMST difference 口径，
#    估计 beta-blocker 相对不用药的 36 个月 RMST 差值及其 95%CI；
# 4. 最终输出：
#    - 临床可解释变量的重编码表；
#    - CART 树图与终端节点结果表；
#    - 五个临床变量分层后的森林图与统计表。

required_packages <- c(
  "caret",
  "dplyr",
  "ggplot2",
  "grf",
  "openxlsx",
  "purrr",
  "readr",
  "rpart",
  "rpart.plot",
  "RColorBrewer",
  "stringr",
  "survival",
  "tibble",
  "tidyr"
)

load_required_packages <- function(packages) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0) {
    stop(
      paste0(
        "以下 R 包尚未安装：",
        paste(missing_packages, collapse = ", "),
        "\n请先安装缺失包后重新运行脚本。"
      ),
      call. = FALSE
    )
  }

  invisible(
    purrr::walk(
      packages,
      ~ suppressPackageStartupMessages(
        library(.x, character.only = TRUE)
      )
    )
  )
}

load_required_packages(required_packages)
rm(required_packages)
rm(load_required_packages)

project_dir <- getwd()

# 这里把效应估计阶段的权重策略单独参数化，便于在
# 1) 当前样本内重估 IPTW；
# 2) 直接使用输入数据自带的 `weight` 列；
# 3) 完全不加权；
# 三种做法之间快速切换。
effect_weighting_strategy <- "input_weight"
input_weight_col <- "weight"

valid_weighting_strategies <- c("reestimate_iptw", "input_weight", "unweighted")

if (!effect_weighting_strategy %in% valid_weighting_strategies) {
  stop(
    paste0(
      "`effect_weighting_strategy` 必须是以下之一：",
      paste(valid_weighting_strategies, collapse = ", ")
    ),
    call. = FALSE
  )
}

# 统一设定输入与输出路径，方便后续直接复跑
input_data_path <- file.path(
  project_dir,
  "paper figures and tables_code",
  "df_imputed_new.xlsx"
)
csf_model_path <- file.path(project_dir, "中间分析结果数据", "csf_fit.rds")
output_dir <- file.path(
  project_dir,
  "analysis_outputs",
  paste0("clinical_rule_", effect_weighting_strategy)
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 统一沿用研究主分析的 RMST 截断时间点与随机种子。
tau_rmst <- 36
n_boot_rmst_ci <- 200
set.seed(2026)

# 这里明确列出 CSF 训练时使用的 44 个候选变量。
# 后续在各临床亚组内重新估计 IPTW 时，也沿用这套协变量集合，
# 以尽量与当前研究的主分析框架保持一致。
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

# 按当前研究的 SHAP 结果，固定 clinical rule 分析中使用的前五个可解释变量。
clinical_rule_vars <- c("LVEF", "age", "sex", "Killip_class", "HR")

# 当前脚本将优先用固定 `csf_fit` 重新预测 ITE
analysis_data <- openxlsx::read.xlsx(input_data_path) %>%
  as_tibble() %>%
  mutate(
    BB_use = case_when(
      BB_use %in% c("Yes", 1, "1") ~ 1,
      BB_use %in% c("No", 0, "0") ~ 0,
      TRUE ~ NA_real_
    ),
    mace = case_when(
      mace %in% c("Happened", 1, "1") ~ 1,
      mace %in% c("Not Happened", 0, "0") ~ 0,
      TRUE ~ NA_real_
    )
  )

# 检查当前分析是否具备森林图与 CART 的必需字段，避免后续中途报错。
required_columns <- c(clinical_rule_vars, "ITE", "BB_use", "survival_time", "mace")
missing_columns <- setdiff(required_columns, names(analysis_data))

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "clinical rule 分析缺少以下必需字段：",
      paste(missing_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

# ------------------------------ #
# Step 1. Create interpretable covariates
# ------------------------------ #
# 这里参考 Komura 文献的“interprable covariates”思路，
# 对 SHAP 前五个变量做临床上更容易解释的分组。
recode_interpretable_covariates <- function(data) {
  data %>%
    mutate(
      # 本研究已排除 LVEF < 40% 的患者，因此这里将 LVEF 分为 40%-49% 与 >=50%。
      LVEF_rule = factor(
        case_when(
          is.na(LVEF) ~ NA_character_,
          LVEF < 50 ~ "40%-49%",
          LVEF >= 50 ~ ">=50%"
        ),
        levels = c("40%-49%", ">=50%")
      ),
      # 年龄分层采用老年心血管研究中常见的 <65 / 65-74 / >=75 岁。
      age_rule = factor(
        case_when(
          is.na(age) ~ NA_character_,
          age < 65 ~ "<65 years",
          age >= 65 & age < 75 ~ "65-74 years",
          age >= 75 ~ ">=75 years"
        ),
        levels = c("<65 years", "65-74 years", ">=75 years")
      ),
      # 性别直接保留为二分类变量。
      sex_rule = factor(
        case_when(
          str_to_lower(as.character(sex)) == "female" ~ "Female",
          str_to_lower(as.character(sex)) == "male" ~ "Male",
          TRUE ~ NA_character_
        ),
        levels = c("Female", "Male")
      ),
      # Killip 分级为临床常用的 I / II / III-IV。
      Killip_rule = factor(
        case_when(
          as.character(Killip_class) %in% c("1", "I") ~ "Class I",
          as.character(Killip_class) %in% c("2", "II") ~ "Class II",
          as.character(Killip_class) %in% c("3", "III", "4", "IV") ~ "Class III-IV",
          TRUE ~ NA_character_
        ),
        levels = c("Class I", "Class II", "Class III-IV")
      ),
      # 心率分层采用具有直接临床解释性的 <70 / 70-89 / >=90 bpm。
      HR_rule = factor(
        case_when(
          is.na(HR) ~ NA_character_,
          HR < 70 ~ "<70 bpm",
          HR >= 70 & HR < 90 ~ "70-89 bpm",
          HR >= 90 ~ ">=90 bpm"
        ),
        levels = c("<70 bpm", "70-89 bpm", ">=90 bpm")
      )
    )
}

analysis_data <- recode_interpretable_covariates(analysis_data)

# 读取已拟合的 CSF 模型，并重建与原模型训练一致的预处理流程。
# 当前输入 `df_imputed_new.xlsx` 已确认包含完整 44 个训练协变量，
# 因此这里按主分析一致的 dummy 编码与标准化流程，用固定 `csf_fit` 重预测 ITE。
# `has_full_csf_inputs` 仍保留为防御性检查，避免未来输入源切换时脚本直接失效。
csf_fit <- NULL
dummy_model_csf <- NULL
scale_center_csf <- NULL
scale_scale_csf <- NULL
has_full_csf_inputs <- all(x_vars %in% names(analysis_data))

if (has_full_csf_inputs) {
  csf_fit <- readRDS(csf_model_path)

  X_reference_raw <- analysis_data[, x_vars, drop = FALSE]
  dummy_model_csf <- caret::dummyVars(" ~ .", data = X_reference_raw, fullRank = TRUE)
  X_reference_preprocessed <- predict(dummy_model_csf, newdata = X_reference_raw)
  X_reference_scaled <- scale(X_reference_preprocessed)
  scale_center_csf <- attr(X_reference_scaled, "scaled:center")
  scale_scale_csf <- attr(X_reference_scaled, "scaled:scale")
  scale_scale_csf[!is.finite(scale_scale_csf) | scale_scale_csf == 0] <- 1
}

preprocess_csf_newdata <- function(newdata,
                                   feature_columns,
                                   dummy_model,
                                   scale_center,
                                   scale_scale) {
  raw_x <- newdata[, feature_columns, drop = FALSE]
  encoded_x <- predict(dummy_model, newdata = raw_x)

  scaled_x <- sweep(encoded_x, 2, scale_center, FUN = "-")
  scaled_x <- sweep(scaled_x, 2, scale_scale, FUN = "/")
  storage.mode(scaled_x) <- "double"

  scaled_x
}

# 这里非常关键：
# 当前输入数据具备完整的 CSF 训练协变量，
# 因此树模型拟合、节点点估计和节点区间统一使用固定 `csf_fit` 重预测 ITE。
# 同时保留文件中原有的 `ITE` 到 `ITE_from_file`，仅作为核对或回溯备用。
analysis_data <- analysis_data %>%
  mutate(
    ITE_from_file = suppressWarnings(as.numeric(ITE)),
    ITE = suppressWarnings(as.numeric(ITE))
  )

if (has_full_csf_inputs) {
  analysis_x_scaled <- preprocess_csf_newdata(
    newdata = analysis_data,
    feature_columns = x_vars,
    dummy_model = dummy_model_csf,
    scale_center = scale_center_csf,
    scale_scale = scale_scale_csf
  )

  analysis_data <- analysis_data %>%
    mutate(
      ITE = as.numeric(
        predict(csf_fit, newdata = analysis_x_scaled)[["predictions"]]
      )
    )
}

# 生成类似文献 Table S6 的重编码表。
coding_table <- tibble(
  Covariates = c("LVEF", "Age", "Sex", "Killip class", "Heart rate"),
  `Format used in ITE prediction` = c(
    "Continuous",
    "Continuous",
    "Binary categorical",
    "Ordinal categorical",
    "Continuous"
  ),
  `Format re-coded for subgroup discovery` = c(
    "2 groups: 40%-49% or >=50%",
    "3 groups: <65 years, 65-74 years, or >=75 years",
    "2 groups: Female or Male",
    "3 groups: Class I, Class II, or Class III-IV",
    "3 groups: <70 bpm, 70-89 bpm, or >=90 bpm"
  )
)


# Step 2. Weighted RMST helpers ####
# 下面几组函数服务于“在指定样本内根据预设策略准备权重，
# 并据此计算 36 个月 RMST difference”。

estimate_stabilized_iptw <- function(data, covariates) {
  available_covariates <- covariates[covariates %in% names(data)]

  analysis_df <- data %>%
    mutate(
      BB_use = as.numeric(BB_use),
      mace = as.numeric(mace)
    ) %>%
    filter(
      !is.na(BB_use),
      !is.na(survival_time),
      !is.na(mace)
    )

  if (nrow(analysis_df) == 0) {
    stop("当前样本没有可用于估计 IPTW 的有效观测。", call. = FALSE)
  }

  if (dplyr::n_distinct(analysis_df$BB_use) < 2) {
    stop("当前样本仅包含单一治疗组，无法估计 IPTW。", call. = FALSE)
  }

  analysis_df <- analysis_df %>%
    filter(
      if_all(all_of(available_covariates), ~ !is.na(.x))
    )

  varying_covariates <- available_covariates[
    vapply(
      analysis_df[available_covariates],
      function(column) dplyr::n_distinct(column) > 1,
      logical(1)
    )
  ]

  if (length(varying_covariates) == 0) {
    p_treat <- mean(analysis_df$BB_use == 1)
    return(
      analysis_df %>%
        mutate(
          ps = p_treat,
          weight = 1
        )
    )
  }

  model_df <- analysis_df %>%
    select(BB_use, all_of(varying_covariates)) %>%
    mutate(
      across(where(is.character), as.factor)
    )

  ps_formula <- as.formula(
    paste("BB_use ~", paste(varying_covariates, collapse = " + "))
  )

  ps_fit <- suppressWarnings(
    glm(ps_formula, data = model_df, family = binomial())
  )

  ps <- predict(ps_fit, type = "response")
  ps <- pmin(pmax(ps, 0.01), 0.99)

  p_treat <- mean(model_df$BB_use == 1)
  sw <- ifelse(
    model_df$BB_use == 1,
    p_treat / ps,
    (1 - p_treat) / (1 - ps)
  )

  sw_lower <- quantile(sw, probs = 0.01, na.rm = TRUE)
  sw_upper <- quantile(sw, probs = 0.99, na.rm = TRUE)
  sw <- pmin(pmax(sw, sw_lower), sw_upper)

  analysis_df %>%
    mutate(
      ps = as.numeric(ps),
      weight = as.numeric(sw)
    )
}

prepare_effect_analysis_data <- function(data,
                                         covariates = x_vars,
                                         weighting_strategy = effect_weighting_strategy,
                                         input_weight_name = input_weight_col) {
  base_df <- data %>%
    mutate(
      BB_use = as.numeric(BB_use),
      mace = as.numeric(mace)
    ) %>%
    filter(
      !is.na(BB_use),
      !is.na(survival_time),
      !is.na(mace)
    )

  if (nrow(base_df) == 0) {
    stop("当前样本没有可用于效应估计的有效观测。", call. = FALSE)
  }

  if (weighting_strategy == "reestimate_iptw") {
    return(
      estimate_stabilized_iptw(
        data = base_df,
        covariates = covariates
      )
    )
  }

  if (weighting_strategy == "input_weight") {
    if (!input_weight_name %in% names(base_df)) {
      stop(
        paste0("输入数据中不存在权重列：", input_weight_name),
        call. = FALSE
      )
    }

    return(
      base_df %>%
        mutate(
          weight = as.numeric(.data[[input_weight_name]])
        ) %>%
        filter(
          !is.na(weight),
          is.finite(weight),
          weight > 0
        )
    )
  }

  if (weighting_strategy == "unweighted") {
    return(
      base_df %>%
        mutate(weight = 1)
    )
  }

  stop("未识别的权重策略。", call. = FALSE)
}

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

  if (nrow(subgroup_data) < 5) {
    stop("当前治疗组样本量过小，无法稳定估计加权 RMST。", call. = FALSE)
  }

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

compute_weighted_rmst_difference_point <- function(data,
                                                   subgroup_name,
                                                   tau = 36,
                                                   ps_covariates = x_vars) {
  raw_n <- nrow(data)
  raw_treated_n <- sum(data$BB_use == 1, na.rm = TRUE)
  raw_control_n <- sum(data$BB_use == 0, na.rm = TRUE)

  analysis_df <- prepare_effect_analysis_data(
    data = data,
    covariates = ps_covariates,
    weighting_strategy = effect_weighting_strategy,
    input_weight_name = input_weight_col
  ) %>%
    filter(
      !is.na(survival_time),
      !is.na(mace),
      !is.na(BB_use),
      !is.na(weight),
      is.finite(weight),
      weight > 0
    )

  if (nrow(analysis_df) < 10 || dplyr::n_distinct(analysis_df$BB_use) < 2) {
    return(
      tibble(
        subgroup = subgroup_name,
        n = raw_n,
        treated_n = raw_treated_n,
        control_n = raw_control_n,
        est = NA_real_,
        lci = NA_real_,
        uci = NA_real_
      )
    )
  }

  rmst_by_arm <- bind_rows(
    compute_weighted_group_rmst(analysis_df, arm_value = 1, tau = tau),
    compute_weighted_group_rmst(analysis_df, arm_value = 0, tau = tau)
  )

  difference_estimate <- rmst_by_arm$rmst[rmst_by_arm$arm == 1] -
    rmst_by_arm$rmst[rmst_by_arm$arm == 0]
  difference_se <- sqrt(sum(rmst_by_arm$rmst_se^2))
  z_value <- qnorm(0.975)

  tibble(
    subgroup = subgroup_name,
    n = raw_n,
    treated_n = raw_treated_n,
    control_n = raw_control_n,
    weighting_strategy = effect_weighting_strategy,
    est = as.numeric(difference_estimate),
    lci = as.numeric(difference_estimate - z_value * difference_se),
    uci = as.numeric(difference_estimate + z_value * difference_se)
  )
}

# 使用非参数 bootstrap 百分位法估计 subgroup-specific RMST difference 的 95%CI。
# 这样可以同时反映重抽样、重新估计 IPTW 以及 RMST 积分的不确定性，
# 通常会比单纯的解析式近似区间更保守，也更接近主分析的展示风格。
compute_weighted_rmst_difference <- function(data,
                                             subgroup_name,
                                             tau = 36,
                                             ps_covariates = x_vars,
                                             n_boot = n_boot_rmst_ci,
                                             conf_level = 0.95,
                                             seed_value = 2026) {
  point_estimate <- compute_weighted_rmst_difference_point(
    data = data,
    subgroup_name = subgroup_name,
    tau = tau,
    ps_covariates = ps_covariates
  )

  if (is.na(point_estimate$est[[1]]) || nrow(data) < 10) {
    return(point_estimate)
  }

  set.seed(seed_value)

  bootstrap_draws <- purrr::map_dbl(
    seq_len(n_boot),
    function(boot_id) {
      sampled_index <- sample(
        x = seq_len(nrow(data)),
        size = nrow(data),
        replace = TRUE
      )

      sampled_data <- data[sampled_index, , drop = FALSE]

      tryCatch(
        compute_weighted_rmst_difference_point(
          data = sampled_data,
          subgroup_name = subgroup_name,
          tau = tau,
          ps_covariates = ps_covariates
        )$est[[1]],
        error = function(e) NA_real_
      )
    }
  )

  valid_draws <- bootstrap_draws[is.finite(bootstrap_draws)]

  if (length(valid_draws) < max(30, ceiling(n_boot * 0.6))) {
    return(point_estimate)
  }

  alpha <- 1 - conf_level

  point_estimate %>%
    mutate(
      lci = as.numeric(quantile(valid_draws, probs = alpha / 2, na.rm = TRUE)),
      uci = as.numeric(quantile(valid_draws, probs = 1 - alpha / 2, na.rm = TRUE)),
      bootstrap_n = length(valid_draws)
    )
}

# ------------------------------ #
# Step 3. Clinical subgroup forest table
# ------------------------------ #
# 这一部分不是用 CART 终端节点，而是直接展示五个临床变量各自亚组的获益情况，
# 便于做正文/补充材料中的森林图。
subgroup_spec <- list(
  list(
    variable = "LVEF",
    rule_var = "LVEF_rule",
    levels = c("40%-49%", ">=50%")
  ),
  list(
    variable = "Age",
    rule_var = "age_rule",
    levels = c("<65 years", "65-74 years", ">=75 years")
  ),
  list(
    variable = "Sex",
    rule_var = "sex_rule",
    levels = c("Female", "Male")
  ),
  list(
    variable = "Killip class",
    rule_var = "Killip_rule",
    levels = c("Class I", "Class II", "Class III-IV")
  ),
  list(
    variable = "Heart rate",
    rule_var = "HR_rule",
    levels = c("<70 bpm", "70-89 bpm", ">=90 bpm")
  )
)

overall_result <- compute_weighted_rmst_difference(
  data = analysis_data,
  subgroup_name = "Overall population",
  tau = tau_rmst,
  ps_covariates = x_vars
) %>%
  mutate(
    variable = "Overall",
    subgroup_label = "Overall population"
  )

clinical_subgroup_results <- purrr::map_dfr(
  subgroup_spec,
  function(spec) {
    purrr::map_dfr(
      spec$levels,
      function(level_value) {
        subgroup_df <- analysis_data %>%
          filter(.data[[spec$rule_var]] == level_value)

        compute_weighted_rmst_difference(
          data = subgroup_df,
          subgroup_name = level_value,
          tau = tau_rmst,
          ps_covariates = x_vars
        ) %>%
          mutate(
            variable = spec$variable,
            subgroup_label = level_value
          )
      }
    )
  }
)

clinical_forest_table <- bind_rows(
  overall_result %>%
    select(
      variable,
      subgroup_label,
      n,
      treated_n,
      control_n,
      weighting_strategy,
      est,
      lci,
      uci
    ),
  clinical_subgroup_results %>%
    select(
      variable,
      subgroup_label,
      n,
      treated_n,
      control_n,
      weighting_strategy,
      est,
      lci,
      uci
    )
) %>%
  mutate(
    variable = factor(
      variable,
      levels = c("Overall", "LVEF", "Age", "Sex", "Killip class", "Heart rate")
    )
  ) %>%
  arrange(variable)

# 为森林图准备分面顺序。
forest_plot_data <- clinical_forest_table %>%
  mutate(
    variable = as.character(variable),
    subgroup_label = case_when(
      variable == "Overall" ~ "Overall population",
      TRUE ~ subgroup_label
    )
  ) %>%
  group_by(variable) %>%
  mutate(
    subgroup_label = factor(subgroup_label, levels = rev(unique(subgroup_label))),
    benefit_direction = case_when(
      is.na(est) ~ "Unavailable",
      est >= 0 ~ "Favors beta-blocker",
      est < 0 ~ "Favors no beta-blocker"
    )
  ) %>%
  ungroup()

clinical_forest_plot <- ggplot(
  forest_plot_data,
  aes(x = est, y = subgroup_label, xmin = lci, xmax = uci, color = benefit_direction)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.7
  ) +
  geom_errorbar(
    aes(ymin = subgroup_label, ymax = subgroup_label),
    orientation = "y",
    width = 0.18,
    linewidth = 0.75,
    na.rm = TRUE
  ) +
  geom_point(
    size = 2.8,
    na.rm = TRUE
  ) +
  facet_grid(
    variable ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_color_manual(
    values = c(
      "Favors beta-blocker" = "firebrick",
      "Favors no beta-blocker" = "steelblue",
      "Unavailable" = "gray60"
    ),
    drop = FALSE
  ) +
  labs(
    x = "36-month RMST difference (month)",
    y = NULL,
    color = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.7),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

# ------------------------------ #
# Step 4. CART model for clinical rule discovery
# ------------------------------ #
# 这里直接用重编码后的五个可解释变量拟合 CART，
# 目标值为 CSF 已经估计出的 ITE。
cart_data <- analysis_data %>%
  filter(
    !is.na(ITE),
    !is.na(LVEF_rule),
    !is.na(age_rule),
    !is.na(sex_rule),
    !is.na(Killip_rule),
    !is.na(HR_rule)
  )

cart_formula <- ITE ~ LVEF_rule + age_rule + sex_rule + Killip_rule + HR_rule

cart_model <- rpart::rpart(
  formula = cart_formula,
  data = cart_data,
  method = "anova",
  control = rpart::rpart.control(
    minbucket = ceiling(0.01 * nrow(cart_data)),
    maxdepth = 2,
    cp = 0
  )
)

# 提取每个节点对应的路径规则，
# 并据此在原始 cart_data 中回溯每个节点包含的样本，
# 以便统一计算节点层面的 observed weighted RMST difference 及其 bootstrap 95%CI。
all_node_ids <- as.integer(row.names(cart_model$frame))
all_node_paths <- path.rpart(cart_model, nodes = all_node_ids, print.it = FALSE)

subset_by_node_path <- function(data, node_path) {
  filtered_data <- data
  path_terms <- node_path[node_path != "root"]

  if (length(path_terms) == 0) {
    return(filtered_data)
  }

  for (single_term in path_terms) {
    if (!grepl("=", single_term, fixed = TRUE)) {
      next
    }

    rule_var <- trimws(sub("=.*$", "", single_term))
    rule_rhs <- trimws(sub("^[^=]+=", "", single_term))
    rule_levels <- strsplit(rule_rhs, ",", fixed = TRUE)[[1]]
    rule_levels <- trimws(rule_levels)

    filtered_data <- filtered_data %>%
      filter(as.character(.data[[rule_var]]) %in% rule_levels)
  }

  filtered_data
}

# 对每个节点直接估计 observed weighted RMST difference 及其 bootstrap 百分位区间。
# 这样树图节点与森林图、终端节点结果表会采用完全一致的效应口径，
# 不再混用“平均预测 ITE 的窄区间”与“observed RMST difference 的宽区间”。
node_effect_summary <- tibble(
  node_id = all_node_ids,
  node_path = unname(all_node_paths)
) %>%
  mutate(
    node_data = purrr::map(node_path, ~ subset_by_node_path(cart_data, .x)),
    node_sample_n = purrr::map_int(node_data, nrow),
    mean_ite = purrr::map_dbl(node_data, ~ mean(.x$ITE, na.rm = TRUE)),
    rmst_result = purrr::map2(
      node_data,
      node_id,
      ~ compute_weighted_rmst_difference(
        data = .x,
        subgroup_name = paste0("Node ", .y),
        tau = tau_rmst,
        ps_covariates = x_vars,
        n_boot = n_boot_rmst_ci,
        conf_level = 0.95,
        seed_value = 2026 + .y
      )
    )
  ) %>%
  mutate(
    rmst_est = purrr::map_dbl(rmst_result, ~ .x$est[[1]]),
    rmst_lci = purrr::map_dbl(rmst_result, ~ .x$lci[[1]]),
    rmst_uci = purrr::map_dbl(rmst_result, ~ .x$uci[[1]]),
    rmst_bootstrap_n = purrr::map_dbl(
      rmst_result,
      ~ if ("bootstrap_n" %in% names(.x)) as.numeric(.x$bootstrap_n[[1]]) else NA_real_
    )
  ) %>%
  select(
    node_id,
    node_sample_n,
    mean_ite,
    rmst_est,
    rmst_lci,
    rmst_uci,
    rmst_bootstrap_n
  )

# 将每个节点对应的效应汇总结果写回 cart_model$frame，
# 这样树图绘制时可以直接读取“当前节点”的数值，
# 避免在 node.fun 中按节点号二次查表造成错配。
cart_model$frame <- cart_model$frame %>%
  rownames_to_column(var = "node_id") %>%
  as_tibble() %>%
  mutate(node_id = as.integer(node_id)) %>%
  left_join(node_effect_summary, by = "node_id") %>%
  mutate(node_id = as.character(node_id)) %>%
  tibble::column_to_rownames(var = "node_id")

format_split_label <- function(label_text) {
  label_text <- gsub("LVEF_rule", "LVEF", label_text)
  label_text <- gsub("age_rule", "Age", label_text)
  label_text <- gsub("sex_rule", "Sex", label_text)
  label_text <- gsub("Killip_rule", "Killip class", label_text)
  label_text <- gsub("HR_rule", "Heart rate", label_text)

  if (!grepl("=", label_text, fixed = TRUE)) {
    return(label_text)
  }

  lhs <- trimws(sub("=.*$", "", label_text))
  rhs <- trimws(sub("^[^=]+=", "", label_text))

  # 将 rpart 对分类变量的“多水平 yes 侧”显示，转成更直观的临床阈值语言。
  if (lhs == "Heart rate") {
    if (rhs == "<70 bpm,70-89 bpm") {
      return("Heart rate <90 bpm")
    }
    if (rhs == "70-89 bpm,>=90 bpm") {
      return("Heart rate >=70 bpm")
    }
    if (rhs == "<70 bpm") {
      return("Heart rate <70 bpm")
    }
    if (rhs == ">=90 bpm") {
      return("Heart rate >=90 bpm")
    }
  }

  if (lhs == "Age") {
    if (rhs == "<65 years,65-74 years") {
      return("Age <75 years")
    }
    if (rhs == "65-74 years,>=75 years") {
      return("Age >=65 years")
    }
  }

  if (lhs == "LVEF") {
    if (rhs == ">=50%") {
      return("LVEF >=50%")
    }
    if (rhs == "40%-49%") {
      return("LVEF 40%-49%")
    }
  }

  paste(lhs, rhs)
}

# 这段绘图函数是对你师姐版本的本研究化改写：
# 变量名映射成当前研究实际使用的 5 个临床变量，并保持树图节点显示
# “observed weighted RMST difference”及其 bootstrap CI 与样本量。
fancy_rpart_plot <- function(model,
                             main = "",
                             sub = NULL,
                             palettes = NULL,
                             type = 2,
                             ...) {
  if (is.null(sub)) {
    sub <- ""
  }

  default_palettes <- c("Greens", "Blues", "Oranges", "Purples", "Reds", "Greys")

  if (is.null(palettes)) {
    palettes <- default_palettes
  }

  missed <- setdiff(1:6, seq_along(palettes))
  palettes <- c(palettes, default_palettes[missed])

  num_pals <- 6
  pal_size <- 5
  pals <- c(
    RColorBrewer::brewer.pal(9, palettes[1])[1:5],
    RColorBrewer::brewer.pal(9, palettes[2])[1:5],
    RColorBrewer::brewer.pal(9, palettes[3])[1:5],
    RColorBrewer::brewer.pal(9, palettes[4])[1:5],
    RColorBrewer::brewer.pal(9, palettes[5])[1:5],
    RColorBrewer::brewer.pal(9, palettes[6])[1:5]
  )

  per <- as.numeric(model$frame$yval / max(model$frame$yval))
  col_index <- round(per * (pal_size - 1)) + 1
  col_index <- abs(col_index)
  total_n <- max(model$frame$n)

  rpart.plot::prp(
    model,
    type = type,
    yesno = 2,
    node.fun = function(x, labs, digits, varlen) {
      paste0(
        "RMST difference = ", round(x$frame$rmst_est, 2),
        "\n95% bootstrap CI ", round(x$frame$rmst_lci, 2),
        " to ", round(x$frame$rmst_uci, 2),
        "\nN = ", x$frame$n,
        " (", round(x$frame$n / total_n * 100, 1), "%)"
      )
    },
    split.fun = function(x, labs, digits, varlen, faclen) {
      labs <- vapply(labs, format_split_label, character(1))
      labs
    },
    box.col = pals[col_index],
    digit = 3,
    nn = TRUE,
    varlen = 0,
    faclen = 0,
    shadow.col = 0,
    fallen.leaves = TRUE,
    branch.lty = 1,
    ...
  )

  title(main = main, sub = sub)
}

# 候选树浏览只需要快速比较不同规则结构，
# 因此这里额外提供一个轻量版树图函数，仅展示节点预测值与样本量，
# 避免为每棵候选树都重复计算较耗时的 bootstrap CI。
fancy_rpart_plot_candidate <- function(model,
                                       main = "",
                                       sub = NULL,
                                       palettes = NULL,
                                       type = 2,
                                       ...) {
  if (is.null(sub)) {
    sub <- ""
  }

  default_palettes <- c("Greens", "Blues", "Oranges", "Purples", "Reds", "Greys")

  if (is.null(palettes)) {
    palettes <- default_palettes
  }

  missed <- setdiff(1:6, seq_along(palettes))
  palettes <- c(palettes, default_palettes[missed])

  pal_size <- 5
  pals <- c(
    RColorBrewer::brewer.pal(9, palettes[1])[1:5],
    RColorBrewer::brewer.pal(9, palettes[2])[1:5],
    RColorBrewer::brewer.pal(9, palettes[3])[1:5],
    RColorBrewer::brewer.pal(9, palettes[4])[1:5],
    RColorBrewer::brewer.pal(9, palettes[5])[1:5],
    RColorBrewer::brewer.pal(9, palettes[6])[1:5]
  )

  per <- as.numeric(model$frame$yval / max(model$frame$yval))
  col_index <- round(per * (pal_size - 1)) + 1
  col_index <- abs(col_index)
  total_n <- max(model$frame$n)

  rpart.plot::prp(
    model,
    type = type,
    yesno = 2,
    node.fun = function(x, labs, digits, varlen) {
      paste0(
        "Predicted ITE = ", round(x$frame$yval, 2),
        "\nN = ", x$frame$n,
        " (", round(x$frame$n / total_n * 100, 1), "%)"
      )
    },
    split.fun = function(x, labs, digits, varlen, faclen) {
      labs <- vapply(labs, format_split_label, character(1))
      labs
    },
    box.col = pals[col_index],
    digit = 3,
    nn = TRUE,
    varlen = 0,
    faclen = 0,
    shadow.col = 0,
    fallen.leaves = TRUE,
    branch.lty = 1,
    ...
  )

  title(main = main, sub = sub)
}

cart_frame_table <- cart_model$frame %>%
  rownames_to_column(var = "node_id") %>%
  as_tibble() %>%
  mutate(node_id = as.integer(node_id))

# 提取各终端节点的 if-then 规则路径，便于输出结果表。
terminal_nodes <- cart_frame_table %>%
  filter(var == "<leaf>") %>%
  pull(node_id)
terminal_paths <- path.rpart(cart_model, nodes = terminal_nodes, print.it = FALSE)

terminal_rule_table <- tibble(
  node_id = as.integer(names(terminal_paths)),
  rule = purrr::map_chr(
    terminal_paths,
    function(path_text) {
      cleaned_path <- path_text[path_text != "root"]

      if (length(cleaned_path) == 0) {
        return("All patients")
      }

      cleaned_path <- gsub("^root", "", cleaned_path)
      cleaned_path <- vapply(cleaned_path, format_split_label, character(1))
      paste(cleaned_path, collapse = " & ")
    }
  )
)

# 在每个 CART 终端节点内估计加权 RMST difference 与 95%CI，
# 对应文献中在“subgroups suggested by CART”内再做效应估计的步骤。
cart_terminal_results <- purrr::map_dfr(
  terminal_nodes,
  function(node_value) {
    node_df <- subset_by_node_path(
      data = cart_data,
      node_path = terminal_paths[[as.character(node_value)]]
    )

    compute_weighted_rmst_difference(
      data = node_df,
      subgroup_name = paste0("Node ", node_value),
      tau = tau_rmst,
      ps_covariates = x_vars
    ) %>%
      mutate(node_id = as.integer(node_value))
  }
) %>%
  left_join(terminal_rule_table, by = "node_id") %>%
  left_join(
    node_effect_summary %>%
      select(node_id, mean_ite, rmst_est, rmst_lci, rmst_uci),
    by = "node_id"
  ) %>%
  rename(
    mean_predicted_ite = mean_ite,
    node_weighted_rmst_diff = rmst_est,
    node_weighted_rmst_diff_lci = rmst_lci,
    node_weighted_rmst_diff_uci = rmst_uci,
    weighted_rmst_diff = est,
    weighted_rmst_diff_lci = lci,
    weighted_rmst_diff_uci = uci
  ) %>%
  select(
    node_id,
    rule,
    n,
    treated_n,
    control_n,
    mean_predicted_ite,
    node_weighted_rmst_diff,
    node_weighted_rmst_diff_lci,
    node_weighted_rmst_diff_uci,
    weighted_rmst_diff,
    weighted_rmst_diff_lci,
    weighted_rmst_diff_uci
  ) %>%
  arrange(node_id)

# ------------------------------ #
# Step 5. Save outputs
# ------------------------------ #
readr::write_csv(
  coding_table,
  file.path(output_dir, "clinical_rule_coding_table.csv")
)

readr::write_csv(
  clinical_forest_table,
  file.path(output_dir, "clinical_rule_forest_table.csv")
)

readr::write_csv(
  cart_frame_table,
  file.path(output_dir, "clinical_rule_cart_frame.csv")
)

readr::write_csv(
  cart_terminal_results,
  file.path(output_dir, "clinical_rule_cart_terminal_results.csv")
)

clinical_rule_workbook <- openxlsx::createWorkbook()

openxlsx::addWorksheet(clinical_rule_workbook, "Coding_table")
openxlsx::writeData(clinical_rule_workbook, "Coding_table", coding_table)

openxlsx::addWorksheet(clinical_rule_workbook, "Forest_table")
openxlsx::writeData(clinical_rule_workbook, "Forest_table", clinical_forest_table)

openxlsx::addWorksheet(clinical_rule_workbook, "CART_frame")
openxlsx::writeData(clinical_rule_workbook, "CART_frame", cart_frame_table)

openxlsx::addWorksheet(clinical_rule_workbook, "CART_terminal")
openxlsx::writeData(clinical_rule_workbook, "CART_terminal", cart_terminal_results)

openxlsx::saveWorkbook(
  clinical_rule_workbook,
  file.path(output_dir, "clinical_rule_outputs.xlsx"),
  overwrite = TRUE
)

ggsave(
  filename = file.path(output_dir, "clinical_rule_forest_plot.png"),
  plot = clinical_forest_plot,
  width = 8.5,
  height = 10,
  dpi = 600,
  bg = "white"
)

png(
  filename = file.path(output_dir, "clinical_rule_cart_tree.png"),
  width = 2200,
  height = 1600,
  res = 220,
  bg = "white"
)
par(bg = "white", mar = c(3, 1, 3, 1))
fancy_rpart_plot(cart_model, sub = NULL)
dev.off()

# ------------------------------ #
# Step 6. Generate multiple candidate CART trees
# ------------------------------ #
# 说明：
# 对同一份数据，rpart 通常近乎确定性，只改 seed 往往不会真正改变树结构。
# 因此这里采用“不同 seed × 不同抽样比例”的重抽样方案来生成多棵候选树，
# 让 CART 在较小或中等样本比例下更有机会选中 Age / Sex / Killip class，
# 供人工挑选更有临床解释性的 rule 结构。
candidate_tree_dir <- file.path(output_dir, "candidate_trees")
dir.create(candidate_tree_dir, recursive = TRUE, showWarnings = FALSE)

candidate_tree_config <- tibble(
  candidate_id = seq_len(12),
  seed = c(2026:2031, 2032:2037),
  sample_fraction = c(1.00, 0.90, 0.80, 0.70, 0.60, 0.50, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40)
)

extract_terminal_rule_table <- function(model_object) {
  frame_df <- model_object$frame %>%
    rownames_to_column(var = "node_id") %>%
    as_tibble() %>%
    mutate(node_id = as.integer(node_id))

  terminal_node_ids <- frame_df %>%
    filter(var == "<leaf>") %>%
    pull(node_id)

  terminal_paths_local <- path.rpart(
    model_object,
    nodes = terminal_node_ids,
    print.it = FALSE
  )

  tibble(
    node_id = as.integer(names(terminal_paths_local)),
    rule = purrr::map_chr(
      terminal_paths_local,
      function(path_text) {
        cleaned_path <- path_text[path_text != "root"]

        if (length(cleaned_path) == 0) {
          return("All patients")
        }

        cleaned_path <- gsub("^root", "", cleaned_path)
        cleaned_path <- vapply(cleaned_path, format_split_label, character(1))
        paste(cleaned_path, collapse = " & ")
      }
    )
  )
}

candidate_tree_summary <- purrr::pmap_dfr(
  candidate_tree_config,
  function(candidate_id, seed, sample_fraction) {
    seed_value <- seed

    set.seed(seed_value)

    sampled_size <- max(
      ceiling(sample_fraction * nrow(cart_data)),
      ceiling(0.20 * nrow(cart_data))
    )

    sampled_index <- sample(
      x = seq_len(nrow(cart_data)),
      size = sampled_size,
      replace = TRUE
    )

    candidate_data <- cart_data[sampled_index, , drop = FALSE]

    candidate_model <- rpart::rpart(
      formula = cart_formula,
      data = candidate_data,
      method = "anova",
      control = rpart::rpart.control(
        minbucket = ceiling(0.01 * nrow(candidate_data)),
        maxdepth = 2,
        cp = 0,
        model = TRUE
      )
    )

    candidate_rules <- extract_terminal_rule_table(candidate_model)

    split_variables <- candidate_model$frame$var[
      candidate_model$frame$var != "<leaf>"
    ] %>%
      unique() %>%
      vapply(format_split_label, character(1))

    png(
      filename = file.path(
        candidate_tree_dir,
        sprintf(
          "candidate_tree_%02d_seed_%d_frac_%s.png",
          candidate_id,
          seed_value,
          formatC(sample_fraction, format = "f", digits = 2)
        )
      ),
      width = 2200,
      height = 1600,
      res = 220,
      bg = "white"
    )
    par(bg = "white", mar = c(3, 1, 3, 1))
    fancy_rpart_plot_candidate(
      candidate_model,
      main = paste0("Candidate tree ", sprintf("%02d", candidate_id)),
      sub = paste0(
        "Bootstrap sample with seed ",
        seed_value,
        "; sample fraction = ",
        formatC(sample_fraction, format = "f", digits = 2)
      )
    )
    dev.off()

    candidate_rules %>%
      mutate(
        candidate_id = candidate_id,
        seed = seed_value,
        sample_fraction = sample_fraction,
        sampled_n = nrow(candidate_data),
        split_variables = paste(split_variables, collapse = " | "),
        n_terminal_nodes = nrow(candidate_rules)
      ) %>%
      select(
        candidate_id,
        seed,
        sample_fraction,
        sampled_n,
        n_terminal_nodes,
        split_variables,
        node_id,
        rule
      )
  }
)

readr::write_csv(
  candidate_tree_summary,
  file.path(candidate_tree_dir, "candidate_tree_summary.csv")
)

# ------------------------------ #
# Step 7. Generate candidate CART trees by minbucket-maxdepth grid
# ------------------------------ #
# 这一组候选树不再依赖不同抽样比例，
# 而是固定 5 个可解释变量与全样本数据，
# 系统比较不同 minbucket 与 maxdepth 组合下的树结构。
# 这样更接近临床上“挑一棵解释性更好的树”的过程。
candidate_tree_grid_dir <- file.path(output_dir, "candidate_trees_grid")
dir.create(candidate_tree_grid_dir, recursive = TRUE, showWarnings = FALSE)

candidate_tree_grid_config <- tidyr::expand_grid(
  minbucket_fraction = c(0.005, 0.01, 0.02, 0.03),
  maxdepth = c(1L, 2L, 3L)
) %>%
  mutate(grid_id = row_number(), .before = 1)

grid_tree_summary <- purrr::pmap_dfr(
  candidate_tree_grid_config,
  function(grid_id, minbucket_fraction, maxdepth) {
    minbucket_value <- max(
      ceiling(minbucket_fraction * nrow(cart_data)),
      20L
    )

    grid_model <- rpart::rpart(
      formula = cart_formula,
      data = cart_data,
      method = "anova",
      control = rpart::rpart.control(
        minbucket = minbucket_value,
        maxdepth = maxdepth,
        cp = 0,
        model = TRUE
      )
    )

    grid_rules <- extract_terminal_rule_table(grid_model)

    split_variables <- grid_model$frame$var[
      grid_model$frame$var != "<leaf>"
    ] %>%
      unique() %>%
      vapply(format_split_label, character(1))

    png(
      filename = file.path(
        candidate_tree_grid_dir,
        sprintf(
          "grid_tree_%02d_minbucket_%s_maxdepth_%d.png",
          grid_id,
          formatC(minbucket_fraction, format = "f", digits = 3),
          maxdepth
        )
      ),
      width = 2200,
      height = 1600,
      res = 220,
      bg = "white"
    )
    par(bg = "white", mar = c(3, 1, 3, 1))
    fancy_rpart_plot_candidate(
      grid_model,
      main = paste0("Grid tree ", sprintf("%02d", grid_id)),
      sub = paste0(
        "minbucket = ",
        minbucket_value,
        " (",
        formatC(minbucket_fraction * 100, format = "f", digits = 1),
        "%), maxdepth = ",
        maxdepth
      )
    )
    dev.off()

    grid_rules %>%
      mutate(
        grid_id = grid_id,
        minbucket_fraction = minbucket_fraction,
        minbucket_n = minbucket_value,
        maxdepth = maxdepth,
        n_terminal_nodes = nrow(grid_rules),
        split_variables = paste(split_variables, collapse = " | ")
      ) %>%
      select(
        grid_id,
        minbucket_fraction,
        minbucket_n,
        maxdepth,
        n_terminal_nodes,
        split_variables,
        node_id,
        rule
      )
  }
)

readr::write_csv(
  grid_tree_summary,
  file.path(candidate_tree_grid_dir, "candidate_tree_grid_summary.csv")
)

print(coding_table)
print(clinical_forest_table)
print(cart_terminal_results)

