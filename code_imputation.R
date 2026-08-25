# data preprocessing ####

required_packages <- c(
  "dplyr",
  "ggplot2",
  "mice",
  "openxlsx",
  "purrr",
  "readr",
  "stringr",
  "survival",
  "tidyr"
)

# 定义依赖包加载函数
load_required_packages <- function(packages) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  # 若存在未安装包，则直接停止运行并给出清晰提示。
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

  # 若所有依赖包均已存在，则逐个加载并抑制启动提示信息。
  invisible(
    purrr::walk(
      packages,
      ~ suppressPackageStartupMessages(
        library(.x, character.only = TRUE)
      )
    )
  )
}

# 实际执行依赖包加载。
load_required_packages(required_packages)
rm(required_packages)
rm(load_required_packages)

# 1. 路径与输出目录设置
# 将当前工作目录视为项目根目录，保证脚本在 IDE 或终端中运行时路径一致。
project_dir <- getwd()

# 所有规范化后的结果统一写入一个单独的输出目录，避免污染根目录。
output_dir <- file.path(project_dir, "analysis_outputs", "mice_survival")

# 若输出目录不存在，则自动创建。
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 输入文件路径
server_csv_template <- "\\\\server\\share\\cardiology\\project\\df_original.csv"
Sys.setenv(PKUH3_RESCUER_SOURCE_CSV = "\\\\server\\share\\cardiology\\project\\df_original.csv")
input_path <- Sys.getenv(
  "PKUH3_RESCUER_SOURCE_CSV",
  unset = server_csv_template
)


# 2. 工具函数：数据读取与列名标准化
# 使用统一的列名清洗函数，将空格、连字符等统一转换为 snake_case。
clean_names_local <- function(x) {
  x %>%
    stringr::str_trim() %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("(^_+|_+$)", "") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_to_lower()
}

# 统一的数值解析函数，可将字符型数值安全转换为数值。
parse_numeric_safe <- function(x) {
  readr::parse_number(
    as.character(x),
    na = c("", "NA", "NaN", "NULL", "Missing", "missing")
  )
}

# 统一的二分类变量重编码函数，将不同写法归一到 Yes/No 两水平因子。
recode_yes_no <- function(x) {
  x_chr <- as.character(x) %>%
    stringr::str_trim() %>%
    stringr::str_to_lower()

  dplyr::case_when(
    x_chr %in% c("1", "yes", "y", "true") ~ "Yes",
    x_chr %in% c("0", "no", "n", "false") ~ "No",
    TRUE ~ NA_character_
  ) %>%
    factor(levels = c("No", "Yes"))
}

# 性别变量标准化为 male/female 两水平因子。
recode_sex <- function(x) {
  x_chr <- as.character(x) %>%
    stringr::str_trim() %>%
    stringr::str_to_lower()

  dplyr::case_when(
    x_chr %in% c("0", "male", "m") ~ "male",
    x_chr %in% c("1", "female", "f") ~ "female",
    TRUE ~ NA_character_
  ) %>%
    factor(levels = c("male", "female"))
}

# MI 类型变量标准化为 STEMI / Non-STEMI。
recode_mi_type <- function(x) {
  x_chr <- as.character(x) %>%
    stringr::str_trim() %>%
    stringr::str_to_lower()

  dplyr::case_when(
    x_chr %in% c("1", "stemi") ~ "STEMI",
    x_chr %in% c("0", "non-stemi", "nstemi", "non stemi") ~ "Non-STEMI",
    TRUE ~ NA_character_
  ) %>%
    factor(levels = c("Non-STEMI", "STEMI"))
}

# Killip 分级标准化为有序因子，便于后续 MICE 采用有序多分类插补。
recode_killip_class <- function(x) {
  x_chr <- as.character(x) %>%
    stringr::str_trim() %>%
    stringr::str_to_lower()

  dplyr::case_when(
    x_chr %in% c("1", "class 1", "i") ~ "Class 1",
    x_chr %in% c("2", "class 2", "ii") ~ "Class 2",
    x_chr %in% c("3", "class 3", "iii") ~ "Class 3",
    x_chr %in% c("4", "class 4", "iv") ~ "Class 4",
    TRUE ~ NA_character_
  ) %>%
    factor(
      levels = c("Class 1", "Class 2", "Class 3", "Class 4"),
      ordered = TRUE
    )
}

# 血运重建方式标准化为 None / PCI / CABG 三分类因子。
recode_revascularization <- function(x) {
  x_chr <- as.character(x) %>%
    stringr::str_trim() %>%
    stringr::str_to_lower()

  dplyr::case_when(
    x_chr %in% c("0", "none") ~ "None",
    x_chr %in% c("1", "pci") ~ "PCI",
    x_chr %in% c("2", "cabg") ~ "CABG",
    TRUE ~ NA_character_
  ) %>%
    factor(levels = c("None", "PCI", "CABG"))
}

# MACE 结局标准化为 Not Happened / Happened。
recode_mace <- function(x) {
  x_chr <- as.character(x) %>%
    stringr::str_trim() %>%
    stringr::str_to_lower()

  dplyr::case_when(
    x_chr %in% c("0", "not happened", "no", "none") ~ "Not Happened",
    x_chr %in% c("1", "happened", "yes", "event") ~ "Happened",
    TRUE ~ NA_character_
  ) %>%
    factor(levels = c("Not Happened", "Happened"))
}

# 按文件类型读取数据；目前同时支持 xlsx 与 csv。
read_analysis_data <- function(path) {
  file_ext <- tools::file_ext(path) %>%
    stringr::str_to_lower()

  if (file_ext == "xlsx") {
    openxlsx::read.xlsx(path, sheet = 1)
  } else if (file_ext == "csv") {
    readr::read_csv(path, show_col_types = FALSE)
  } else {
    stop("当前脚本仅支持读取 .xlsx 或 .csv 文件。", call. = FALSE)
  }
}

# 3. 工具函数：数据标准化
# 该函数负责：
# 1. 清洗列名；
# 2. 统一变量类型；
# 3. 删除旧脚本遗留的派生列，仅保留分析真正需要的变量。
standardize_analysis_data <- function(df_raw) {
  names(df_raw) <- clean_names_local(names(df_raw))

  # 定义分析中使用的变量名
  numeric_vars <- c(
    "age", "bmi", "lvef", "egfr", "hemoglobin", "ldl", "sodium", "potassium",
    "hr", "sbp", "hs_ctnt", "hs_ctni", "ck_mb", "myoglobin", "nt_probnp",
    "hba1c", "wbc", "plt", "inr", "aptt", "pt", "scr", "survival_time",
    "follow_time"
  )

  binary_vars <- c(
    "smoking_status", "drinking_status", "hypertension", "diabetes", "stroke",
    "hypercholesterolemia", "prior_hf", "prior_cad", "prior_pad", "bradycardia",
    "conduction_disease", "copd", "asthma", "prior_bb", "multivessel_disease",
    "asprin", "p2y12_receptor_blocker", "acei_or_arb", "diuretics", "ccb",
    "mra", "ivabradine", "ezetimibe", "anticoagulants", "statin"
  )

  # 通过 `any_of()` 处理数据中实际存在的列，避免因为字段略有差异而报错。
  df_clean <- df_raw %>%
    dplyr::mutate(
      dplyr::across(dplyr::any_of(numeric_vars), parse_numeric_safe)
    )

  # 对关键字符变量逐一标准化。
  if ("id" %in% names(df_clean)) {
    df_clean <- df_clean %>% dplyr::mutate(id = as.character(id))
  }

  if ("sex" %in% names(df_clean)) {
    df_clean <- df_clean %>% dplyr::mutate(sex = recode_sex(sex))
  }

  if ("bb_use" %in% names(df_clean)) {
    df_clean <- df_clean %>% dplyr::mutate(bb_use = recode_yes_no(bb_use))
  }

  if ("mi_type" %in% names(df_clean)) {
    df_clean <- df_clean %>% dplyr::mutate(mi_type = recode_mi_type(mi_type))
  }

  if ("killip_class" %in% names(df_clean)) {
    df_clean <- df_clean %>% dplyr::mutate(killip_class = recode_killip_class(killip_class))
  }

  if ("revascularization" %in% names(df_clean)) {
    df_clean <- df_clean %>% dplyr::mutate(revascularization = recode_revascularization(revascularization))
  }

  if ("mace" %in% names(df_clean)) {
    df_clean <- df_clean %>% dplyr::mutate(mace = recode_mace(mace))
  }

  # 对二分类共变量统一映射到 Yes/No。
  available_binary_vars <- intersect(binary_vars, names(df_clean))

  if (length(available_binary_vars) > 0) {
    df_clean <- df_clean %>%
      dplyr::mutate(
        dplyr::across(dplyr::all_of(available_binary_vars), recode_yes_no)
      )
  }
  df_clean
}

# 4. 工具函数：缺失概览与 MICE 设置
# 生成变量层面的缺失汇总表。
make_missing_summary <- function(df) {
  df %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::everything(),
        ~ sum(is.na(.x))
      )
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "variable",
      values_to = "missing_n"
    ) %>%
    dplyr::mutate(
      missing_pct = 100 * missing_n / nrow(df)
    ) %>%
    dplyr::arrange(dplyr::desc(missing_n))
}

# 根据变量类型构建 MICE 方法向量。
build_mice_method <- function(df, protected_columns = character()) {
  methods <- mice::make.method(df)

  # 对不需要插补的变量明确设置为空方法。
  methods[protected_columns] <- ""

  # 其余变量按类型指定方法。
  for (column_name in setdiff(names(df), protected_columns)) {
    column_data <- df[[column_name]]

    if (is.numeric(column_data)) {
      methods[column_name] <- "pmm"
    } else if (is.factor(column_data) && nlevels(column_data) == 2) {
      methods[column_name] <- "logreg"
    } else if (is.ordered(column_data)) {
      methods[column_name] <- "polr"
    } else if (is.factor(column_data) && nlevels(column_data) > 2) {
      methods[column_name] <- "polyreg"
    } else {
      methods[column_name] <- ""
    }
  }

  methods
}

# 构建预测矩阵，并确保 `id` 不作为预测变量。
build_predictor_matrix <- function(df, blocked_columns = character()) {
  predictor_matrix <- mice::make.predictorMatrix(df)

  for (column_name in blocked_columns) {
    if (column_name %in% rownames(predictor_matrix)) {
      predictor_matrix[column_name, ] <- 0
    }

    if (column_name %in% colnames(predictor_matrix)) {
      predictor_matrix[, column_name] <- 0
    }
  }

  predictor_matrix
}



# Survival analysis ####

# 估计sIPTW
compute_stabilized_iptw <- function(df, treatment_var, covariates) {
  usable_covariates <- covariates[
    vapply(
      df[covariates],
      function(x) dplyr::n_distinct(stats::na.omit(x)) > 1,
      logical(1)
    )
  ]

  if (length(usable_covariates) == 0) {
    stop("未找到可用于 PS 模型的有效协变量。", call. = FALSE)
  }

  df_ps <- df %>%
    dplyr::mutate(
      treat_indicator = as.integer(.data[[treatment_var]] == "Yes")
    )

  ps_formula <- stats::reformulate(
    termlabels = usable_covariates,
    response = "treat_indicator"
  )

  ps_model <- stats::glm(
    formula = ps_formula,
    data = df_ps,
    family = stats::binomial()
  )

  ps_hat <- stats::predict(ps_model, type = "response") %>%
    pmin(0.99) %>%
    pmax(0.01)

  treatment_probability <- mean(df_ps$treat_indicator == 1)

  iptw_raw <- dplyr::if_else(
    condition = df_ps$treat_indicator == 1,
    true = treatment_probability / ps_hat,
    false = (1 - treatment_probability) / (1 - ps_hat)
  )

  lower_cutoff <- stats::quantile(iptw_raw, probs = 0.01, na.rm = TRUE)
  upper_cutoff <- stats::quantile(iptw_raw, probs = 0.99, na.rm = TRUE)

  df_ps %>%
    dplyr::mutate(
      ps = ps_hat,
      iptw = pmin(pmax(iptw_raw, lower_cutoff), upper_cutoff)
    )
}

# 拟合加权 Cox 模型，并提取 log(HR) 及其方差。
fit_weighted_cox_single <- function(df, treatment_var, time_var, event_var) {
  df_model <- df %>%
    dplyr::mutate(
      event_indicator = as.integer(.data[[event_var]] == "Happened")
    )

  cox_formula <- stats::as.formula(
    paste0(
      "survival::Surv(",
      time_var,
      ", event_indicator) ~ ",
      treatment_var
    )
  )

  cox_fit <- survival::coxph(
    formula = cox_formula,
    data = df_model,
    weights = iptw,
    robust = TRUE,
    ties = "efron"
  )

  coefficient_name <- paste0(treatment_var, "Yes")

  if (!coefficient_name %in% names(stats::coef(cox_fit))) {
    stop("Cox 模型未能识别治疗变量的比较系数。", call. = FALSE)
  }

  log_hr <- unname(stats::coef(cox_fit)[coefficient_name])
  log_hr_var <- unname(stats::vcov(cox_fit)[coefficient_name, coefficient_name])

  tibble::tibble(
    estimate = log_hr,
    variance = log_hr_var
  )
}

# 提取加权 KM 曲线在指定时间点上的生存概率及标准误
extract_weighted_km_single <- function(df, treatment_var, time_var, event_var, time_grid) {
  df_model <- df %>%
    dplyr::mutate(
      event_indicator = as.integer(.data[[event_var]] == "Happened")
    )

  km_formula <- stats::as.formula(
    paste0(
      "survival::Surv(",
      time_var,
      ", event_indicator) ~ ",
      treatment_var
    )
  )

  km_fit <- survival::survfit(
    formula = km_formula,
    data = df_model,
    weights = iptw
  )

  km_summary <- summary(km_fit, times = time_grid, extend = TRUE)

  tibble::tibble(
    time = km_summary$time,
    strata = km_summary$strata,
    surv = km_summary$surv,
    surv_se = km_summary$std.err
  ) %>%
    dplyr::mutate(
      bb_use = stringr::str_remove(strata, paste0("^", treatment_var, "="))
    ) %>%
    dplyr::select(time, bb_use, surv, surv_se)
}

# Rubin规则合并
# pooled estimate = q_bar
# total variance = within-imputation variance + (1 + 1/m) * between-imputation variance
# 95% CI 使用对称 t 法；若自由度趋于无穷，则退化为对称 Z 法。
pool_rubin_scalar <- function(q, u, conf_level = 0.95, transform = c("identity", "exp")) {
  transform <- match.arg(transform)

  m <- length(q)
  q_bar <- mean(q)
  u_bar <- mean(u)
  b_var <- if (m > 1) stats::var(q) else 0
  total_var <- u_bar + (1 + 1 / m) * b_var
  total_se <- sqrt(total_var)

  # Rubin 规则的相对增量方差。
  relative_increase <- if (u_bar > 0) {
    ((1 + 1 / m) * b_var) / u_bar
  } else {
    Inf
  }

  # 采用 Rubin 常用自由度近似；若无法稳定计算，则退化为无穷自由度。
  degrees_of_freedom <- if (is.finite(relative_increase) && relative_increase > 0) {
    (m - 1) * (1 + 1 / relative_increase)^2
  } else {
    Inf
  }

  alpha <- 1 - conf_level
  critical_value <- if (is.finite(degrees_of_freedom)) {
    stats::qt(1 - alpha / 2, df = degrees_of_freedom)
  } else {
    stats::qnorm(1 - alpha / 2)
  }

  lower_raw <- q_bar - critical_value * total_se
  upper_raw <- q_bar + critical_value * total_se

  if (transform == "exp") {
    estimate <- exp(q_bar)
    conf_low <- exp(lower_raw)
    conf_high <- exp(upper_raw)
  } else {
    estimate <- q_bar
    conf_low <- lower_raw
    conf_high <- upper_raw
  }

  tibble::tibble(
    estimate = estimate,
    std_error = total_se,
    conf_low = conf_low,
    conf_high = conf_high,
    df = degrees_of_freedom,
    within_var = u_bar,
    between_var = b_var,
    total_var = total_var
  )
}

# 读取数据并生成缺失汇总
# 读取原始或候选数据文件
# 从院内服务器路径读入原始数据，并沿用既往脚本中的对象名 `df_original`。
df_original <- read_analysis_data(input_path)

# 统一标准化列名与变量类型。
analysis_df <- standardize_analysis_data(df_original)

# 定义用于 PS 与生存分析的基线协变量集
available_covariates <- c(
  "sex", "age", "bmi", "smoking_status", "drinking_status", "hypertension",
  "diabetes", "stroke", "hypercholesterolemia", "prior_hf", "prior_cad",
  "prior_pad", "bradycardia", "conduction_disease", "copd", "asthma",
  "prior_bb", "mi_type", "killip_class", "multivessel_disease",
  "revascularization", "lvef", "egfr", "hemoglobin", "ldl", "sodium",
  "potassium", "hba1c", "asprin", "p2y12_receptor_blocker", "acei_or_arb",
  "diuretics", "ccb", "mra", "anticoagulants", "statin", "hr", "sbp",
  "hs_ctni", "nt_probnp", "wbc", "plt", "inr", "aptt", "pt", "scr"
)

# 明确检查分析所需核心字段是否存在
required_columns <- c("bb_use", "survival_time", "mace")
missing_required_columns <- setdiff(required_columns, names(analysis_df))

if (length(missing_required_columns) > 0) {
  stop(
    paste0(
      "输入数据缺少关键分析字段：",
      paste(missing_required_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

if (length(available_covariates) < 10) {
  stop(
    paste0(
      "当前数据中可识别的基线协变量过少（仅 ",
      length(available_covariates),
      " 个），无法稳定完成 PS / IPTW 建模。"
    ),
    call. = FALSE
  )
}

# 输出原始缺失概况，便于核对缺失模式。
missing_summary_raw <- make_missing_summary(analysis_df)

readr::write_csv(
  missing_summary_raw,
  file.path(output_dir, "missing_summary_raw.csv")
)

# 输出缺失条形图，图形化查看缺失程度
missing_plot <- missing_summary_raw %>%
  dplyr::filter(missing_n > 0) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = stats::reorder(variable, missing_pct),
      y = missing_pct
    )
  ) +
  ggplot2::geom_col(fill = "steelblue") +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = NULL,
    y = "Missingness (%)"
  ) +
  ggplot2::theme_minimal(base_size = 12)

ggplot2::ggsave(
  filename = file.path(output_dir, "missingness_profile.png"),
  plot = missing_plot,
  width = 8,
  height = 10,
  dpi = 600
)

# 按研究方案执行排除与 MICE
# 计算“基线协变量层面”的个体缺失比例，用于实施 >20% 缺失患者排除。
analysis_df_filtered <- analysis_df %>%
  dplyr::mutate(
    baseline_missing_fraction = rowMeans(
      dplyr::across(
        dplyr::all_of(available_covariates),
        ~ is.na(.x)
      )
    )
  )

# 按顺序执行核心字段排除，统计 >20% 基线缺失排除人数，
analysis_df_after_core_filter <- analysis_df_filtered %>%
  dplyr::filter(
    !is.na(bb_use),
    !is.na(survival_time),
    !is.na(mace)
  )

excluded_missing_core_n <- nrow(analysis_df_filtered) - nrow(analysis_df_after_core_filter)

excluded_excess_missing_n <- sum(
  analysis_df_after_core_filter$baseline_missing_fraction > 0.20,
  na.rm = TRUE
)

# 记录顺序式队列流程
cohort_flow <- tibble::tribble(
  ~step, ~n,
  "Raw input cohort", nrow(analysis_df),
  "Excluded: missing treatment / time / event", excluded_missing_core_n,
  "Cohort after core-field exclusion", nrow(analysis_df_after_core_filter),
  "Excluded: >20% baseline covariate missingness", excluded_excess_missing_n
) %>%
  dplyr::bind_rows(
    tibble::tibble(
      step = "Final analysis cohort",
      n = nrow(analysis_df_after_core_filter) - excluded_excess_missing_n
    )
  )

# 删除不满足分析条件的观察值
analysis_df_filtered <- analysis_df_after_core_filter %>%
  dplyr::filter(
    baseline_missing_fraction <= 0.20
  )

# 保留分析必需列
mice_input_df <- analysis_df_filtered %>%
  dplyr::select(
    dplyr::any_of(c("id", "bb_use", "survival_time", "mace")),
    dplyr::all_of(available_covariates)
  )

# 生成插补前分析样本的缺失概况
missing_summary_analysis <- make_missing_summary(mice_input_df)

readr::write_csv(
  missing_summary_analysis,
  file.path(output_dir, "missing_summary_analysis_set.csv")
)

# 指定不需要插补的列：
# 1. id 仅作标识；
# 2. bb_use / survival_time / mace 为治疗与结局字段，不在本脚本中插补。
protected_columns <- intersect(
  c("id", "bb_use", "survival_time", "mace"),
  names(mice_input_df)
)

# 构建 MICE 所需的方法向量与预测矩阵
mice_methods <- build_mice_method(
  df = mice_input_df,
  protected_columns = protected_columns
)

mice_predictor_matrix <- build_predictor_matrix(
  df = mice_input_df,
  blocked_columns = c("id")
)

# 实施五重MICE插补
mice_fit <- mice::mice(
  data = mice_input_df,
  m = 5,
  method = mice_methods,
  predictorMatrix = mice_predictor_matrix,
  maxit = 20,   # `maxit = 20` 在大多数临床队列中已足以获得稳定收敛
  seed = 2026,
  printFlag = TRUE
)

# 保存完整mids对象
saveRDS(
  mice_fit,
  file = file.path(output_dir, "mice_fit_object.rds")
)

# 导出 long-format 的插补数据
imputed_long_df <- mice::complete(mice_fit, action = "long", include = FALSE) %>%
  dplyr::rename(imputation_id = .imp)

readr::write_csv(
  imputed_long_df,
  file.path(output_dir, "mice_imputed_long.csv")
)

# 每个插补数据集分别进行 IPTW + Cox + KM
# 定义 KM 曲线的输出时间网格，这里按 0–36 月、每 1 月一个时间点输出
km_time_grid <- seq(0, 36, by = 1)

# 读取 5 个独立插补完成的数据集
imputed_datasets <- mice::complete(mice_fit, action = "all")

# 在每个插补数据集中：
# 1. 估计 PS 与 sIPTW
# 2. 拟合加权 Cox
# 3. 提取加权 KM 生存概率
analysis_results <- purrr::imap(
  imputed_datasets,
  function(single_df, imp_id) {
    weighted_df <- compute_stabilized_iptw(
      df = single_df,
      treatment_var = "bb_use",
      covariates = available_covariates
    )

    cox_result <- fit_weighted_cox_single(
      df = weighted_df,
      treatment_var = "bb_use",
      time_var = "survival_time",
      event_var = "mace"
    ) %>%
      dplyr::mutate(imputation_id = as.integer(imp_id))

    km_result <- extract_weighted_km_single(
      df = weighted_df,
      treatment_var = "bb_use",
      time_var = "survival_time",
      event_var = "mace",
      time_grid = km_time_grid
    ) %>%
      dplyr::mutate(imputation_id = as.integer(imp_id))

    weight_summary <- weighted_df %>%
      dplyr::summarise(
        imputation_id = as.integer(imp_id),
        n = dplyr::n(),
        ps_min = min(ps, na.rm = TRUE),
        ps_max = max(ps, na.rm = TRUE),
        weight_min = min(iptw, na.rm = TRUE),
        weight_median = stats::median(iptw, na.rm = TRUE),
        weight_mean = mean(iptw, na.rm = TRUE),
        weight_max = max(iptw, na.rm = TRUE)
      )

    list(
      weighted_data = weighted_df %>% dplyr::mutate(imputation_id = as.integer(imp_id)),
      cox_result = cox_result,
      km_result = km_result,
      weight_summary = weight_summary
    )
  }
)

# 合并每个插补数据集的加权数据
weighted_long_df <- purrr::map_dfr(analysis_results, "weighted_data")

readr::write_csv(
  weighted_long_df,
  file.path(output_dir, "weighted_imputed_long.csv")
)

# 合并每个插补数据集的 Cox 结果
cox_results_by_imputation <- purrr::map_dfr(analysis_results, "cox_result")

readr::write_csv(
  cox_results_by_imputation,
  file.path(output_dir, "cox_results_by_imputation.csv")
)

# 合并每个插补数据集的权重诊断表
weight_diagnostics <- purrr::map_dfr(analysis_results, "weight_summary")

readr::write_csv(
  weight_diagnostics,
  file.path(output_dir, "weight_diagnostics.csv")
)

# 合并每个插补数据集的 KM 摘要
km_results_by_imputation <- purrr::map_dfr(analysis_results, "km_result")

readr::write_csv(
  km_results_by_imputation,
  file.path(output_dir, "km_results_by_imputation.csv")
)

# 使用Rubin规则合并Cox结果
# 对log(HR)及其方差应用Rubin规则，再指数化得到HR及95% CI。
pooled_cox_result <- pool_rubin_scalar(
  q = cox_results_by_imputation$estimate,
  u = cox_results_by_imputation$variance,
  conf_level = 0.95,
  transform = "exp"
) %>%
  dplyr::mutate(
    effect_measure = "Hazard ratio for beta-blocker use (Yes vs No)",
    ci_method = paste(
      "Rubin rules with symmetric t-based 95% CI on the log(HR) scale,",
      "followed by exponentiation."
    )
  ) %>%
  dplyr::select(
    effect_measure,
    estimate,
    conf_low,
    conf_high,
    std_error,
    df,
    within_var,
    between_var,
    total_var,
    ci_method
  )

readr::write_csv(
  pooled_cox_result,
  file.path(output_dir, "pooled_cox_result.csv")
)

# 使用Rubin规则合并KM生存概率
# 对每个时间点、每个治疗组的生存概率进行合并，生成基于多重插补结果的pooled KM曲线
pooled_km_curve <- km_results_by_imputation %>%
  dplyr::group_by(time, bb_use) %>%
  dplyr::summarise(
    pooled = list(
      pool_rubin_scalar(
        q = surv,
        u = surv_se^2,
        conf_level = 0.95,
        transform = "identity"
      )
    ),
    .groups = "drop"
  ) %>%
  tidyr::unnest_wider(pooled) %>%
  dplyr::mutate(
    estimate = pmin(pmax(estimate, 0), 1),
    conf_low = pmin(pmax(conf_low, 0), 1),
    conf_high = pmin(pmax(conf_high, 0), 1),
    cumulative_incidence = 100 * (1 - estimate),
    cumulative_incidence_low = 100 * (1 - conf_high),
    cumulative_incidence_high = 100 * (1 - conf_low),
    ci_method = "Rubin rules with symmetric t-based pointwise 95% CI on the survival-probability scale."
  )

readr::write_csv(
  pooled_km_curve,
  file.path(output_dir, "pooled_km_curve.csv")
)

# 绘制pooled KM累积事件曲线
pooled_km_plot <- pooled_km_curve %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = time,
      y = cumulative_incidence,
      color = bb_use
    )
  ) +
  ggplot2::geom_step(linewidth = 1) +
  ggplot2::scale_color_manual(
    values = c("No" = "steelblue", "Yes" = "firebrick")
  ) +
  ggplot2::coord_cartesian(xlim = c(0, 36)) +
  ggplot2::labs(
    x = "Follow-up (months)",
    y = "Cumulative incidence (%)",
    color = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "top"
  )

ggplot2::ggsave(
  filename = file.path(output_dir, "pooled_km_curve.png"),
  plot = pooled_km_plot,
  width = 8,
  height = 6,
  dpi = 600
)