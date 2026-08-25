library(dplyr)
library(ggplot2)
library(MatchIt)
library(openxlsx)
library(readr)
library(stringr)
library(survival)

project_dir <- getwd()
output_dir <- file.path(project_dir, "analysis_outputs", "sensitivity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

analysis_path <- Sys.getenv("PKUH3_RESCUER_DF_OVERALL_XLSX")
if (identical(analysis_path, "")) {
  analysis_path <- file.path(project_dir, "df_overall.xlsx")
}

raw_path <- Sys.getenv("PKUH3_RESCUER_SOURCE_CSV")
if (identical(raw_path, "")) {
  stop("Environment variable PKUH3_RESCUER_SOURCE_CSV is required.")
}

read_any_data <- function(path) {
  file_ext <- tools::file_ext(path)
  if (tolower(file_ext) == "csv") {
    return(readr::read_csv(path, show_col_types = FALSE))
  }
  if (tolower(file_ext) %in% c("xlsx", "xls")) {
    return(openxlsx::read.xlsx(path, sheet = 1))
  }
  stop("Unsupported input file format.")
}

clean_names_local <- function(x) {
  x %>%
    stringr::str_trim() %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("(^_+|_+$)", "") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_to_lower()
}

find_first_existing <- function(data, candidates) {
  hits <- candidates[candidates %in% names(data)]
  if (length(hits) == 0) {
    return(NA_character_)
  }
  hits[1]
}

rename_first_existing <- function(data, target, candidates) {
  hit <- find_first_existing(data, candidates)
  if (!is.na(hit) && hit != target) {
    names(data)[names(data) == hit] <- target
  }
  data
}

recode_yes_no_binary <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  ifelse(x_chr %in% c("1", "yes", "y", "true", "happened"), 1L, 0L)
}

prepare_analysis_data <- function(data) {
  names(data) <- clean_names_local(names(data))

  data <- rename_first_existing(data, "bb_use", c("bb_use", "bbuse"))
  data <- rename_first_existing(data, "survival_time", c("survival_time", "time", "event_time"))
  data <- rename_first_existing(data, "mace", c("mace", "event", "status"))
  data <- rename_first_existing(data, "event_mace", c("event_mace", "mace_type"))
  data <- rename_first_existing(data, "ite", c("ite"))
  data <- rename_first_existing(data, "ite_group", c("ite_group", "benefit_group"))
  data <- rename_first_existing(data, "id", c("id", "patient_id"))
  data <- rename_first_existing(data, "ctni", c("ctni", "c_tni", "hs_ctni", "hs_ctni_i"))
  data <- rename_first_existing(data, "nt_probnp", c("nt_probnp", "nt_probnp_pg_ml"))
  data <- rename_first_existing(data, "hr", c("hr", "heart_rate"))
  data <- rename_first_existing(data, "sbp", c("sbp", "systolic_blood_pressure"))
  data <- rename_first_existing(data, "killip_class", c("killip_class", "killip"))
  data <- rename_first_existing(data, "mi_type", c("mi_type", "myocardial_infarction_type"))
  data <- rename_first_existing(data, "prior_hf", c("prior_hf", "prior_heart_failure"))
  data <- rename_first_existing(data, "prior_cad", c("prior_cad"))
  data <- rename_first_existing(data, "prior_pad", c("prior_pad"))
  data <- rename_first_existing(data, "prior_bb", c("prior_bb"))
  data <- rename_first_existing(data, "multivessel_disease", c("multivessel_disease", "multi_vessel", "multi_vessel_disease"))
  data <- rename_first_existing(data, "conduction_disease", c("conduction_disease", "conduction_block"))
  data <- rename_first_existing(data, "acei_or_arb", c("acei_or_arb", "acei_arb"))
  data <- rename_first_existing(data, "p2y12_receptor_blocker", c("p2y12_receptor_blocker"))
  data <- rename_first_existing(data, "asprin", c("asprin", "aspirin"))
  data <- rename_first_existing(data, "ccb", c("ccb"))
  data <- rename_first_existing(data, "mra", c("mra"))
  data <- rename_first_existing(data, "lvef", c("lvef"))
  data <- rename_first_existing(data, "egfr", c("egfr"))
  data <- rename_first_existing(data, "hba1c", c("hba1c", "hba1c_percent"))
  data <- rename_first_existing(data, "wbc", c("wbc", "wbc_count"))
  data <- rename_first_existing(data, "plt", c("plt", "platelet_count"))

  data %>%
    mutate(
      arm = recode_yes_no_binary(bb_use),
      time36 = pmin(as.numeric(survival_time), 36),
      status_mace = ifelse(recode_yes_no_binary(mace) == 1L & as.numeric(survival_time) <= 36, 1L, 0L),
      ite_group = factor(
        as.character(ite_group),
        levels = c("Lower third", "Middle third", "Upper third")
      )
    )
}

merge_raw_columns <- function(analysis_data, raw_data) {
  if (!("id" %in% names(analysis_data) && "id" %in% names(raw_data))) {
    return(analysis_data)
  }

  raw_only <- setdiff(names(raw_data), names(analysis_data))
  if (length(raw_only) == 0) {
    return(analysis_data)
  }

  analysis_data %>%
    left_join(raw_data %>% dplyr::select(id, all_of(raw_only)), by = "id")
}

x_vars <- c(
  "sex", "age", "bmi", "smoking_status", "drinking_status",
  "hypertension", "diabetes", "stroke", "hypercholesterolemia",
  "prior_cad", "prior_pad", "prior_hf", "bradycardia",
  "conduction_disease", "copd", "asthma", "prior_bb", "mi_type",
  "killip_class", "multivessel_disease", "revascularization", "lvef",
  "egfr", "hemoglobin", "ldl", "hba1c", "sodium", "potassium",
  "asprin", "p2y12_receptor_blocker", "acei_or_arb", "diuretics",
  "ccb", "mra", "anticoagulants", "statin", "hr", "sbp", "ctni",
  "nt_probnp", "wbc", "plt", "inr", "aptt"
)

add_ps_weights <- function(data, covariates, method = c("iptw", "overlap")) {
  method <- match.arg(method)
  available_covariates <- covariates[covariates %in% names(data)]

  if (length(available_covariates) == 0) {
    stop("No propensity model covariates available.")
  }

  model_data <- data %>%
    dplyr::select(arm, all_of(available_covariates)) %>%
    mutate(across(where(is.character), as.factor))

  ps_formula <- as.formula(
    paste("arm ~", paste(available_covariates, collapse = " + "))
  )

  ps_fit <- glm(ps_formula, data = model_data, family = binomial())
  ps <- predict(ps_fit, type = "response") %>%
    pmin(0.99) %>%
    pmax(0.01)

  if (method == "iptw") {
    p_treat <- mean(model_data$arm == 1)
    weight <- ifelse(model_data$arm == 1, p_treat / ps, (1 - p_treat) / (1 - ps))
    q_lo <- quantile(weight, probs = 0.01, na.rm = TRUE)
    q_hi <- quantile(weight, probs = 0.99, na.rm = TRUE)
    weight <- pmin(pmax(weight, q_lo), q_hi)
  } else {
    weight <- ifelse(model_data$arm == 1, 1 - ps, ps)
  }

  data %>%
    mutate(
      ps = ps,
      analysis_weight = as.numeric(weight)
    )
}

# 对加权Kaplan-Meier生存曲线在0到tau区间内做阶梯积分，
compute_weighted_km_rmst <- function(fit_object, tau = 36) {
  # survfit(~1) 时通常没有 strata；这里统一处理成单层结构，
  # 便于后续按相同逻辑拆分时间点与生存概率。
  strata_sizes <- fit_object$strata

  if (is.null(strata_sizes)) {
    strata_sizes <- length(fit_object$time)
    names(strata_sizes) <- "All"
  }

  # 按strata还原每一层对应的时间点与生存概率。
  strata_names <- names(strata_sizes)
  strata_index <- rep(strata_names, strata_sizes)
  split_time <- split(fit_object$time, strata_index)
  split_surv <- split(fit_object$surv, strata_index)

  bind_rows(lapply(strata_names, function(strata_name) {
    # 取出当前strata的事件时间与加权生存概率
    time_vector <- split_time[[strata_name]]
    surv_vector <- split_surv[[strata_name]]

    # RMST只在预设的36个月窗口内积分
    valid_index <- time_vector < tau
    integration_time <- c(0, time_vector[valid_index], tau)
    integration_surv <- c(1, surv_vector[valid_index])

    # 对阶梯型KM曲线做左连续积分，面积即RMST
    tibble(
      strata = strata_name,
      rmst = sum(diff(integration_time) * integration_surv)
    )
  }))
}

# 在单个治疗组内拟合加权Kaplan-Meier曲线，使用robust infinitesimal jackknife残差计算RMST标准误与95%CI
compute_group_specific_rmst <- function(data, arm_value, status_col, weight_col, tau = 36) {
  # 仅保留当前治疗组、结局信息完整且权重有效的样本。
  subgroup_data <- data %>%
    filter(
      arm == arm_value,
      !is.na(time36),
      !is.na(.data[[status_col]]),
      !is.na(.data[[weight_col]]),
      is.finite(.data[[weight_col]]),
      .data[[weight_col]] > 0
    ) %>%
    mutate(
      # 先把动态列名标准化成固定列，避免公式环境下 `.data[[...]]` 求值不稳定。
      status_value = as.integer(.data[[status_col]]),
      weight_value = as.numeric(.data[[weight_col]])
    )

  # 若当前治疗组样本量过小，则直接返回缺失结果，避免 survfit 报错。
  if (nrow(subgroup_data) < 5) {
    return(
      tibble(
        arm = arm_value,
        rmst = NA_real_,
        rmst_se = NA_real_,
        rmst_lci = NA_real_,
        rmst_uci = NA_real_
      )
    )
  }

  # 按主分析方案拟合加权 KM 曲线；
  # robust=TRUE 与 se.fit=TRUE 会调用 IJ 方差估计。
  subgroup_fit <- survival::survfit(
    survival::Surv(time36, status_value) ~ 1,
    data = subgroup_data,
    weights = weight_value,
    robust = TRUE,
    se.fit = TRUE,
    conf.type = "plain",
    model = TRUE
  )

  # 通过加权 KM 曲线积分得到该治疗组在 tau 月内的 RMST。
  rmst_estimate <- compute_weighted_km_rmst(
    fit_object = subgroup_fit,
    tau = tau
  ) %>%
    pull(rmst)

  # 直接提取 RMST 的 IJ 残差，并按平方和开根号得到标准误。
  rmst_ij <- residuals(
    subgroup_fit,
    times = tau,
    type = "rmst",
    weighted = TRUE
  )

  rmst_se <- sqrt(sum(as.numeric(rmst_ij)^2, na.rm = TRUE))
  z_975 <- qnorm(0.975)

  tibble(
    arm = arm_value,
    rmst = as.numeric(rmst_estimate),
    rmst_se = as.numeric(rmst_se),
    rmst_lci = as.numeric(rmst_estimate - z_975 * rmst_se),
    rmst_uci = as.numeric(rmst_estimate + z_975 * rmst_se)
  )
}

# 先分别估计两组的加权 RMST，
# 再以两组 RMST 之差作为 treatment effect，并沿用主分析的正态近似区间。
estimate_rmst_difference <- function(data, status_col, weight_col = NULL, tau = 36) {
  # 保留时间、结局与治疗分组完整的分析样本。
  analysis_data <- data %>%
    filter(
      !is.na(time36),
      !is.na(.data[[status_col]]),
      !is.na(arm)
    )

  # 若当前分析没有外部权重，则赋值为 1，
  # 这样 negative control 也能共用同一套 KM + IJ 估计流程。
  if (!is.null(weight_col)) {
    analysis_data <- analysis_data %>%
      filter(
        !is.na(.data[[weight_col]]),
        is.finite(.data[[weight_col]]),
        .data[[weight_col]] > 0
      )
  } else {
    analysis_data <- analysis_data %>%
      mutate(analysis_weight = 1)
    weight_col <- "analysis_weight"
  }

  # 若有效样本不足或仅有单一治疗组，则无法计算 RMST difference。
  if (nrow(analysis_data) < 10 || dplyr::n_distinct(analysis_data$arm) < 2) {
    return(tibble(estimate = NA_real_, lci = NA_real_, uci = NA_real_, p_value = NA_real_))
  }

  # 分别计算 beta-blocker 组与非 beta-blocker 组的组别 RMST 及其 IJ 标准误。
  rmst_by_arm <- bind_rows(
    compute_group_specific_rmst(
      data = analysis_data,
      arm_value = 1,
      status_col = status_col,
      weight_col = weight_col,
      tau = tau
    ),
    compute_group_specific_rmst(
      data = analysis_data,
      arm_value = 0,
      status_col = status_col,
      weight_col = weight_col,
      tau = tau
    )
  )

  # 只要任一治疗组 RMST 无法估计，就返回缺失结果。
  if (any(is.na(rmst_by_arm$rmst)) || any(is.na(rmst_by_arm$rmst_se))) {
    return(tibble(estimate = NA_real_, lci = NA_real_, uci = NA_real_, p_value = NA_real_))
  }

  # treatment effect 定义为治疗组 RMST 减去对照组 RMST。
  difference_estimate <- rmst_by_arm$rmst[rmst_by_arm$arm == 1] - rmst_by_arm$rmst[rmst_by_arm$arm == 0]

  # 两组样本互斥，因此差值标准误按两组方差之和开根号合成。
  difference_se <- sqrt(sum(rmst_by_arm$rmst_se^2, na.rm = TRUE))

  # 使用与主分析一致的 Z 统计量与 95% 置信区间。
  if (!is.finite(difference_se) || difference_se <= 0) {
    return(tibble(estimate = difference_estimate, lci = NA_real_, uci = NA_real_, p_value = NA_real_))
  }

  z_value <- difference_estimate / difference_se
  z_975 <- qnorm(0.975)

  tibble(
    estimate = as.numeric(difference_estimate),
    lci = as.numeric(difference_estimate - z_975 * difference_se),
    uci = as.numeric(difference_estimate + z_975 * difference_se),
    p_value = as.numeric(2 * pnorm(abs(z_value), lower.tail = FALSE))
  )
}

run_analysis_set <- function(data, analysis_name, status_col, weight_col = NULL) {
  cohorts <- c("Overall", "Lower third", "Middle third", "Upper third")

  bind_rows(
    lapply(cohorts, function(cohort_name) {
      cohort_data <- if (cohort_name == "Overall") {
        data
      } else {
        data %>% filter(ite_group == cohort_name)
      }

      estimate_rmst_difference(
        data = cohort_data,
        status_col = status_col,
        weight_col = weight_col,
        tau = 36
      ) %>%
        mutate(
          Analysis = analysis_name,
          Cohort = cohort_name,
          .before = 1
        )
    })
  )
}

resolve_status_column <- function(data, preferred_aliases) {
  hit <- find_first_existing(data, preferred_aliases)
  if (!is.na(hit)) {
    return(hit)
  }

  if ("event_mace" %in% names(data)) {
    event_mace_chr <- tolower(trimws(as.character(data$event_mace)))
    if (any(stringr::str_detect(event_mace_chr, "cardio"), na.rm = TRUE)) {
      data$cardiovascular_death_status <- ifelse(
        stringr::str_detect(event_mace_chr, "cardio"),
        1L,
        0L
      )
      return("cardiovascular_death_status")
    }
  }

  NA_character_
}

analysis_data <- prepare_analysis_data(read_any_data(analysis_path))

if (file.exists(raw_path)) {
  raw_data <- prepare_analysis_data(read_any_data(raw_path))
  analysis_data <- merge_raw_columns(analysis_data, raw_data)
}

analysis_data <- analysis_data %>%
  filter(!is.na(ite_group))

complete_case_data <- analysis_data %>%
  filter(
    if_all(all_of(intersect(c(x_vars, "arm", "time36", "status_mace", "ite_group"), names(analysis_data))), ~ !is.na(.x))
  ) %>%
  add_ps_weights(covariates = x_vars, method = "iptw")

ow_data <- analysis_data %>%
  filter(
    if_all(all_of(intersect(c(x_vars, "arm", "time36", "status_mace", "ite_group"), names(analysis_data))), ~ !is.na(.x))
  ) %>%
  add_ps_weights(covariates = x_vars, method = "overlap")

psm_covariates <- x_vars[x_vars %in% names(analysis_data)]
psm_formula <- as.formula(
  paste("arm ~", paste(psm_covariates, collapse = " + "))
)

psm_input <- analysis_data %>%
  filter(
    if_all(all_of(intersect(c(psm_covariates, "arm", "time36", "status_mace", "ite_group"), names(analysis_data))), ~ !is.na(.x))
  ) %>%
  mutate(across(where(is.character), as.factor))

match_object <- MatchIt::matchit(
  formula = psm_formula,
  data = psm_input,
  method = "nearest",
  ratio = 1,
  replace = FALSE
)

psm_data <- MatchIt::match.data(match_object) %>%
  mutate(analysis_weight = weights)

cvd_status_col <- resolve_status_column(
  analysis_data,
  c("cardiovascular_death_status", "cv_death", "cardiovascular_death", "death_cv")
)

sensitivity_results <- bind_rows(
  run_analysis_set(complete_case_data, "Complete-case", "status_mace", "analysis_weight"),
  run_analysis_set(ow_data, "Overlap weighting", "status_mace", "analysis_weight"),
  run_analysis_set(psm_data, "1:1 PSM", "status_mace", "analysis_weight")
)

if (!is.na(cvd_status_col)) {
  cvd_data <- analysis_data %>%
    filter(
      if_all(all_of(intersect(c(x_vars, "arm", "time36", cvd_status_col, "ite_group"), names(analysis_data))), ~ !is.na(.x))
    ) %>%
    add_ps_weights(covariates = x_vars, method = "iptw")

  sensitivity_results <- bind_rows(
    sensitivity_results,
    run_analysis_set(cvd_data, "Cardiovascular death", cvd_status_col, "analysis_weight")
  )
}

readr::write_csv(
  sensitivity_results,
  file.path(output_dir, "sensitivity_analyses_results.csv")
)


# nco analysis ####
nco_status_col <- resolve_status_column(
  analysis_data,
  c("pneumonia_status", "pneumonia_hospitalization", "pneumonia_event", "noc_status")
)

nco_time_col <- find_first_existing(
  analysis_data,
  c("pneumonia_time", "time_to_pneumonia", "noc_time")
)

if (!is.na(nco_status_col) && !is.na(nco_time_col)) {
  # 在NCO可分析样本上重新定义36个月随访时间与事件指标，避免直接复用主结局MACE的时间变量
  nco_data <- analysis_data %>%
    mutate(
      time36 = pmin(as.numeric(.data[[nco_time_col]]), 36),
      # 仅当NCO在36个月内发生时记为事件，否则按行政截尾处理
      status_nco = ifelse(
        recode_yes_no_binary(.data[[nco_status_col]]) == 1L &
          as.numeric(.data[[nco_time_col]]) <= 36,
        1L,
        0L
      )
    ) %>%
    filter(
      if_all(
        all_of(
          intersect(
            c(x_vars, "arm", "time36", "status_nco", "ite_group"),
            c(names(analysis_data), "time36", "status_nco")
          )
        ),
        ~ !is.na(.x)
      )
    ) %>%
    # 在 NCO 专属分析样本上重算 IPTW，
    # 使权重与当前可分析样本严格对应。
    add_ps_weights(covariates = x_vars, method = "iptw")

  # 输出总体及各 tertile 的加权 36 个月 RMST difference 和 95%CI。
  nco_results <- run_analysis_set(
    data = nco_data,
    analysis_name = "Negative control outcome",
    status_col = "status_nco",
    weight_col = "analysis_weight"
  )

  readr::write_csv(
    nco_results,
    file.path(output_dir, "negative_control_results.csv")
  )

  # 绘图时仅展示三个 tertiles，
  # 与正文中用于比较效应异质性的展示口径保持一致。
  nco_plot_data <- nco_results %>%
    filter(Cohort != "Overall") %>%
    mutate(
      Cohort = factor(
        Cohort,
        levels = c("Lower third", "Middle third", "Upper third")
      )
    )

  nco_plot <- ggplot(
    nco_plot_data,
    aes(x = Cohort, y = estimate, ymin = lci, ymax = uci, color = Cohort)
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 0.8
    ) +
    # 使用点估计与 95%CI 直接呈现各 tertile 的 NCO 结果。
    geom_pointrange(
      linewidth = 0.8,
      fatten = 2.8
    ) +
    scale_color_manual(
      values = c(
        "Lower third" = "steelblue",
        "Middle third" = "gray40",
        "Upper third" = "red"
      ),
      drop = FALSE
    ) +
    labs(
      x = "Patients grouped by predicted individualized treatment effect",
      y = "RMST Difference (month) for Negative Control Outcome"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray75", linewidth = 0.6),
      panel.grid.minor.y = element_blank(),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      axis.line = element_line(color = "black", linewidth = 0.8),
      panel.border = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )

  ggsave(
    filename = file.path(output_dir, "negative_control_rmst_ci.png"),
    plot = nco_plot,
    width = 7,
    height = 6,
    dpi = 600,
    bg = "white"
  )
}
