##########################################################################################################
#Codes for "Reduced mortality burden associated with the clean air policy in China: An integrated assessment based on causal machine learning"
#Authors for codes: Yuchang Zhou, Peizheng Li, Ya Gao, Yi Guo, Yixiang Zhu, Lu Zhou, Peng Yin, Haidong Kan, Maigeng Zhou, Renjie Chen
#Correspondence to Renjie Chen, Maigeng Zhou.
###########################################################################################################
# Case cross-over analysis for Conditional Logistic Regression model
# Output: single + two-pollutant + multi-pollutant models at cummean lag0..2
##########################################################################################################

rm(list = ls()); options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(splines)
})

# ---- I/O ----
data_path <- "/d2/home/user19/LPZ/result/simulated data.rds"
out_dir   <- "/d2/home/user19/LPZ/result/simulated_result"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Global parameters ----
L_met   <- 3L
df_temp <- 6L
df_rh   <- 3L

pollutants  <- c("PM25","O3","PM2510","NO2","SO2","CO")
pollutants2 <- c("PM25","O3","PM2510","NO2","SO2","CO")

target_k     <- 2L
lag_type_tag <- "cummean"

# ✅ avoid recursive promise evaluation
UNIT_STEP_MAP <- list(
  O3     = 1.0,
  PM25   = 1.0,
  PM2510 = 1.0,
  NO2    = 1.0,
  SO2    = 1.0,
  CO     = 1.0
)

# ---- Read & basic preprocessing ----
dt <- readRDS(data_path)
setDT(dt)

# CO unit conversion if your downstream assumes ug/m^3
cols_co <- paste0("CO_lag", 0:7)
if (all(cols_co %in% names(dt))) {
  dt[, (cols_co) := lapply(.SD, function(x) x * 1000), .SDcols = cols_co]
}

# Mean imputation for numeric columns (use set() to avoid selfref warning)
num_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
for (col in num_cols) {
  x <- dt[[col]]
  if (anyNA(x)) {
    m <- mean(x, na.rm = TRUE)
    ii <- which(is.na(x))
    if (length(ii)) data.table::set(dt, i = ii, j = col, value = m)
  }
}

dt[, date    := as.IDate(date)]
dt[, case    := as.integer(case)]
dt[, holiday := as.integer(holiday)]
dt[, id      := as.integer(factor(id))]

need_cols <- function(prefix, L, DT) {
  need <- paste0(prefix, "_lag", 0:L)
  miss <- setdiff(need, names(DT))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
}

need_cols("temp", L_met, dt)
need_cols("rh",   L_met, dt)

dt[, temp_ma := rowMeans(.SD), .SDcols = paste0("temp_lag", 0:L_met)]
dt[, rh_ma   := rowMeans(.SD), .SDcols = paste0("rh_lag",   0:L_met)]

# Build PM2510_lagk from PM10_lagk and PM25_lagk if missing
make_pm2510_lags <- function(DT) {
  pm10_cols <- grep("^PM10_lag\\d+$", names(DT), value = TRUE)
  pm25_cols <- grep("^PM25_lag\\d+$", names(DT), value = TRUE)
  k10 <- as.integer(sub("^PM10_lag", "", pm10_cols))
  k25 <- as.integer(sub("^PM25_lag", "", pm25_cols))
  ks  <- sort(intersect(k10, k25))
  if (!length(ks)) return(invisible(NULL))
  for (k in ks) {
    c10  <- paste0("PM10_lag",  k)
    c25  <- paste0("PM25_lag",  k)
    cnew <- paste0("PM2510_lag", k)
    if (!cnew %in% names(DT)) DT[, (cnew) := pmax(get(c10) - get(c25), 0)]
  }
  invisible(NULL)
}
make_pm2510_lags(dt)

# Ensure cummean (lag0..k)
ensure_cummean <- function(DT, pol, k) {
  need <- paste0(pol, "_lag", 0:k)
  if (!all(need %in% names(DT))) {
    stop(sprintf("[%s] Cannot build cummean 0..%d; missing: %s",
                 pol, k, paste(setdiff(need, names(DT)), collapse = ",")))
  }
  new_name <- paste0(pol, "_ma0_", k)
  if (!new_name %in% names(DT)) {
    DT[, (new_name) := rowMeans(.SD), .SDcols = need]
  }
  new_name
}

# Fit clogit and return tidy one-row results
fit_clogit_and_tidy <- function(dtf, main_var, main_pol, controls_vars, model_type,
                                lag_type, lag_k,
                                df_temp = 6, df_rh = 3,
                                unit_step_map = UNIT_STEP_MAP) {
  
  this <- dtf[, .(case, id, holiday, temp_ma, rh_ma)]
  this[, x := dtf[[main_var]]]
  
  if (length(controls_vars)) {
    for (v in controls_vars) this[, (v) := dtf[[v]]]
  }
  
  check_cols <- c("x", "temp_ma", "rh_ma", controls_vars)
  ok <- Reduce("&", lapply(check_cols, function(cc) is.finite(this[[cc]]))) & is.finite(this$holiday)
  this <- this[ok]
  if (!nrow(this)) return(NULL)
  
  ctrl_terms <- if (length(controls_vars)) paste(controls_vars, collapse = " + ") else NULL
  rhs <- c(
    "x",
    ctrl_terms,
    sprintf("ns(temp_ma, df=%d)", df_temp),
    sprintf("ns(rh_ma, df=%d)", df_rh),
    "holiday",
    "strata(id)"
  )
  fml <- as.formula(paste("case ~", paste(rhs[nzchar(rhs)], collapse = " + ")))
  
  step_val <- if (!is.null(unit_step_map[[main_pol]])) as.numeric(unit_step_map[[main_pol]]) else 1.0
  
  fit <- try(suppressWarnings(clogit(fml, data = this, method = "efron")), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(data.frame(
      main_pollutant = main_pol,
      model_type = model_type,
      controls = if (length(controls_vars)) paste(controls_vars, collapse=",") else "",
      pollutant = main_pol, lag_type = lag_type, lag_k = lag_k,
      exp_var = main_var, n = nrow(this),
      df_temp = df_temp, df_rh = df_rh,
      beta = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_,
      unit_step = step_val,
      OR_unit = NA_real_, LCL_unit = NA_real_, UCL_unit = NA_real_,
      OR_unit_step = NA_real_, LCL_unit_step = NA_real_, UCL_unit_step = NA_real_,
      IQR = NA_real_, OR_IQR = NA_real_, LCL_IQR = NA_real_, UCL_IQR = NA_real_,
      PC_IQR = NA_real_, LCL_PC_IQR = NA_real_, UCL_PC_IQR = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  sm <- summary(fit)
  if (!"x" %in% rownames(sm$coefficients)) {
    beta <- se <- z <- p <- NA_real_
  } else {
    beta <- sm$coefficients["x","coef"]
    se   <- sm$coefficients["x","se(coef)"]
    z    <- sm$coefficients["x","z"]
    p    <- sm$coefficients["x","Pr(>|z|)"]
  }
  
  OR_unit  <- exp(beta)
  LCL_unit <- exp(beta - 1.96*se)
  UCL_unit <- exp(beta + 1.96*se)
  
  OR_unit_step  <- exp(beta * step_val)
  LCL_unit_step <- exp((beta - 1.96*se) * step_val)
  UCL_unit_step <- exp((beta + 1.96*se) * step_val)
  
  IQR_x   <- IQR(this$x, na.rm = TRUE, type = 7)
  OR_IQR  <- exp(beta * IQR_x)
  LCL_IQR <- exp((beta - 1.96*se) * IQR_x)
  UCL_IQR <- exp((beta + 1.96*se) * IQR_x)
  
  PC_IQR <- (OR_IQR - 1) * 100
  LCL_PC <- (LCL_IQR - 1) * 100
  UCL_PC <- (UCL_IQR - 1) * 100
  
  data.frame(
    main_pollutant = main_pol,
    model_type = model_type,
    controls = if (length(controls_vars)) paste(controls_vars, collapse=",") else "",
    pollutant = main_pol, lag_type = lag_type, lag_k = lag_k,
    exp_var = main_var, n = nrow(this),
    df_temp = df_temp, df_rh = df_rh,
    beta = beta, se = se, z = z, p = p,
    unit_step = step_val,
    OR_unit = OR_unit, LCL_unit = LCL_unit, UCL_unit = UCL_unit,
    OR_unit_step = OR_unit_step, LCL_unit_step = LCL_unit_step, UCL_unit_step = UCL_unit_step,
    IQR = IQR_x,
    OR_IQR = OR_IQR, LCL_IQR = LCL_IQR, UCL_IQR = UCL_IQR,
    PC_IQR = PC_IQR, LCL_PC_IQR = LCL_PC, UCL_PC_IQR = UCL_PC,
    stringsAsFactors = FALSE
  )
}

print_progress <- function(res, elapsed_sec) {
  if (is.null(res)) return(invisible(NULL))
  message(sprintf("[%s] %.3f sec | main=%s | model=%s | controls={%s} | lag=%s(%s)",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  elapsed_sec, res$main_pollutant, res$model_type, res$controls,
                  as.character(res$lag_k), res$lag_type))
  message(sprintf("  OR(IQR)=%.4f (95%%CI: %.4f–%.4f), p=%.4g | OR(step=%.3g)=%.4f | n=%d",
                  res$OR_IQR, res$LCL_IQR, res$UCL_IQR, res$p,
                  res$unit_step, res$OR_unit_step, res$n))
}

# ---- Main loop: single + two-pollutant + multi-pollutant (cummean lag0..2 only) ----
overall_t0 <- proc.time()[["elapsed"]]
results <- list(); model_count <- 0L

for (main_pol in pollutants2) {
  pol_t0 <- proc.time()[["elapsed"]]
  
  main_var <- try(ensure_cummean(dt, main_pol, target_k), silent = TRUE)
  if (inherits(main_var, "try-error")) next
  
  # 1) single-pollutant
  t0 <- proc.time()[["elapsed"]]
  res_single <- fit_clogit_and_tidy(
    dtf = dt,
    main_var = main_var, main_pol = main_pol,
    controls_vars = character(0),
    model_type = "single",
    lag_type = lag_type_tag, lag_k = target_k,
    df_temp = df_temp, df_rh = df_rh
  )
  t1 <- proc.time()[["elapsed"]]
  if (!is.null(res_single)) { results[[length(results)+1]] <- res_single; model_count <- model_count + 1L }
  print_progress(res_single, elapsed_sec = t1 - t0)
  
  # 2) two-pollutant models (add one co-pollutant at a time)
  co_pols <- setdiff(pollutants, main_pol)
  for (cp in co_pols) {
    cp_var <- try(ensure_cummean(dt, cp, target_k), silent = TRUE)
    if (inherits(cp_var, "try-error")) next
    
    t0 <- proc.time()[["elapsed"]]
    res_two <- fit_clogit_and_tidy(
      dtf = dt,
      main_var = main_var, main_pol = main_pol,
      controls_vars = cp_var,
      model_type = "two",
      lag_type = lag_type_tag, lag_k = target_k,
      df_temp = df_temp, df_rh = df_rh
    )
    t1 <- proc.time()[["elapsed"]]
    if (!is.null(res_two)) { results[[length(results)+1]] <- res_two; model_count <- model_count + 1L }
    print_progress(res_two, elapsed_sec = t1 - t0)
  }
  
  # 3) multi-pollutant
  multi_ctrl_pols <- setdiff(pollutants, main_pol)
  multi_ctrl_vars <- character(0); ok_all <- TRUE
  for (cp in multi_ctrl_pols) {
    cv <- try(ensure_cummean(dt, cp, target_k), silent = TRUE)
    if (inherits(cv, "try-error")) { ok_all <- FALSE; break }
    multi_ctrl_vars <- c(multi_ctrl_vars, cv)
  }
  
  if (ok_all) {
    t0 <- proc.time()[["elapsed"]]
    res_multi <- fit_clogit_and_tidy(
      dtf = dt,
      main_var = main_var, main_pol = main_pol,
      controls_vars = multi_ctrl_vars,
      model_type = "multi",
      lag_type = lag_type_tag, lag_k = target_k,
      df_temp = df_temp, df_rh = df_rh
    )
    t1 <- proc.time()[["elapsed"]]
    if (!is.null(res_multi)) { results[[length(results)+1]] <- res_multi; model_count <- model_count + 1L }
    print_progress(res_multi, elapsed_sec = t1 - t0)
  }
  
  pol_t1 <- proc.time()[["elapsed"]]
  message(sprintf("<< %s finished | elapsed %.3f sec >>", main_pol, pol_t1 - pol_t0))
}

if (!length(results)) stop("No results were produced. Check required *_lag0..2 columns.")

res_dt <- rbindlist(results, fill = TRUE)
setDT(res_dt)

# BH adjust within main_pollutant (includes single/two/multi together)
res_dt[, p_adj_BH := {
  ps <- p
  if (all(is.na(ps))) rep(NA_real_, .N) else p.adjust(ps, method = "BH")
}, by = .(main_pollutant)]
res_dt[, sig_BH_0.05 := as.integer(!is.na(p_adj_BH) & p_adj_BH < 0.05)]

setcolorder(res_dt, c(
  "main_pollutant","model_type","controls",
  "pollutant","lag_type","lag_k",
  "OR_unit_step","LCL_unit_step","UCL_unit_step",
  "IQR","OR_IQR","LCL_IQR","UCL_IQR",
  "exp_var","n",
  "df_temp","df_rh",
  "beta","se","z","p","p_adj_BH","sig_BH_0.05",
  "unit_step",
  "OR_unit","LCL_unit","UCL_unit",
  "PC_IQR","LCL_PC_IQR","UCL_PC_IQR"
))

outfile <- file.path(
  out_dir,
  sprintf("clogit_single_two_multi_%s_%s_Lmet%d_dfT%d_dfRH%d.csv",
          paste(pollutants2, collapse = "_"),
          "cummean_lag0_2",
          L_met, df_temp, df_rh)
)
fwrite(res_dt, outfile)

overall_t1 <- proc.time()[["elapsed"]]
message(sprintf("== Done | total elapsed %.3f sec | total models fitted: %d ==",
                overall_t1 - overall_t0, model_count))
message("Results saved to: ", outfile)
