##########################################################################################################
#Codes for "Mortality risk and burden associated with all criteria air pollutants in China: Causal machine-learning evidence"
#Authors for codes: Huimeng Liu, Jian Lei, Yunxing Jiang, Lijun Bai, et.al.
#Correspondence to Shaowei Wu, Yuewei Liu.
###########################################################################################################

##########################################################################################################
# Case cross-over analysis for AIPW model (Figure 1, Supplementary Tables 5)
###########################################################################################################

rm(list = ls()); options(stringsAsFactors = FALSE)

# ========= Reproducibility =========
SEED_MASTER <- 20251016L
set.seed(SEED_MASTER)
Sys.setenv(OMP_NUM_THREADS = "10")

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(splines)
  library(xgboost)
  library(jsonlite)
})

# ========= 0) Load data =========
dt <- readRDS("XXX/simulated data.rds")
setDT(dt)

# Create PM2.5-10 = max(PM10 - PM2.5, 0) for each lag
for (k in 0:7) {
  c10  <- sprintf("PM10_lag%d", k)
  c25  <- sprintf("PM25_lag%d", k)
  newc <- sprintf("PM2510_lag%d", k)
  stopifnot(c10 %in% names(dt), c25 %in% names(dt))
  dt[, (newc) := pmax(get(c10) - get(c25), 0)]
}

# ========= 0.2) Simple mean imputation for numeric missing values =========
data_imputed <- dt
for (col in names(data_imputed)) {
  if (is.numeric(data_imputed[[col]])) {
    data_imputed[[col]][is.na(data_imputed[[col]])] <- mean(data_imputed[[col]], na.rm = TRUE)
  }
}
rm(dt)
data <- data_imputed
rm(data_imputed)
stopifnot(exists("data"))
setDT(data)

# Standardize city field + integer city id (used as a covariate)
data[, city := trimws(city)]
data[, city_id := as.integer(factor(city, levels = sort(unique(city))))]

# ========= Global parameters =========
K_fold          <- 3L
tune_try_random <- 10L
DELTA           <- 1.0
LR_CLIP         <- c(1/50, 50)
PIT_BINS        <- 20L

VF_FLOOR_MAP <- c(PM25=0.01, PM2510=0.01, O3=0.01, NO2=0.01, SO2=0.01, CO=0.01)
VF_CAP_MAP   <- c(PM25=0.99, PM2510=0.99, O3=0.99, NO2=0.99, SO2=0.99, CO=0.99)
DELTA_MAP    <- c(PM25=1.0, PM2510=1.0, O3=1.0, NO2=1.0, SO2=1.0, CO=1.0)

tune_sample_n_id <- max(as.integer(ceiling(0.10 * uniqueN(data$id))), 2L)

# ========= Output directory + logging =========
OUT_DIR   <- "XXX/simulated_result/AIPW_RR/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
log_file  <- file.path(OUT_DIR, sprintf("runlog_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
sink(log_file, split = TRUE)

pollutants  <- c("NO2","PM25","O3","SO2","CO","PM2510")
pollutants2 <- c("SO2","CO","NO2","PM25","O3","PM2510")

# ========= Fixed lag specification (main analysis) =========
LAG_PLAN <- list(
  NO2    = list(type = "cum", k = 2),
  PM25   = list(type = "cum", k = 2),
  O3     = list(type = "cum", k = 2),
  PM2510 = list(type = "cum", k = 2),
  SO2    = list(type = "cum", k = 2),
  CO     = list(type = "cum", k = 2)
)
stopifnot(all(pollutants2 %in% names(LAG_PLAN)))

lag_tag_from_mode <- function(mode) {
  stopifnot(is.list(mode), !is.null(mode$type), !is.null(mode$k), mode$k >= 0)
  if (mode$type == "single") sprintf("lag%d", mode$k) else sprintf("lag0%d", mode$k)
}

tic <- function() as.numeric(Sys.time())
toc <- function(t0, msg="") {
  cost <- as.numeric(Sys.time()) - t0
  message(sprintf("[%.2fs] %s", cost, msg))
  invisible(cost)
}

# ========= Utility: within-stratum demeaning based on control-day mean =========
demean_by_control <- function(DT, cols) {
  ctrl_mean <- DT[, lapply(.SD, function(v){
    if (all(is.na(v[case==0]))) mean(v, na.rm = TRUE) else mean(v[case==0], na.rm = TRUE)
  }), by = id, .SDcols = cols]
  setkey(ctrl_mean, id)
  DT <- merge(DT, ctrl_mean, by="id", suffixes = c("", "_ctrlmean"))
  for (cc in cols) {
    mcol <- paste0(cc, "_ctrlmean")
    DT[(!is.na(get(cc))), (cc) := get(cc) - get(mcol)]
    DT[, (mcol) := NULL]
  }
  DT[]
}

# ========= Utility: fold assignment by matched-set id =========
make_folds <- function(ids, K, seed = SEED_MASTER + 1L) {
  K <- min(max(2L, K), length(ids))
  set.seed(seed)
  fid <- sample(rep(1:K, length.out = length(ids)))
  data.table(id = ids, fold = fid)
}

# ========= Utility: memory-safe XGBoost prediction =========
predict_xgb_in_chunks <- function(bst, Xmat, rows_idx, chunk_rows = 200000L) {
  n <- length(rows_idx)
  pred <- numeric(n)
  if (n == 0L) return(pred)
  p <- ncol(Xmat)
  if (is.finite(p) && p > 0) {
    chunk_est <- floor(3e8 / max(8 * p, 1))
    chunk_rows <- max(20000L, min(chunk_rows, chunk_est))
  }
  starts <- seq(1L, n, by = chunk_rows)
  for (s in starts) {
    e <- min(s + chunk_rows - 1L, n)
    idx <- rows_idx[s:e]
    dtmp <- xgb.DMatrix(Xmat[idx, , drop = FALSE])
    pred[s:e] <- predict(bst, dtmp)
    rm(dtmp); gc(FALSE)
  }
  pred
}

# ========= Random search (regression) for exposure model tuning =========
xgb_random_search_reg <- function(dtrain,
                                  n_try = 50L,
                                  nfold_cv = 2L,
                                  nrounds_grid = seq(200, 500, by = 50),
                                  es_rounds = 10L,
                                  seed = SEED_MASTER + 100L) {
  set.seed(seed)
  best <- NULL; best_rmse <- Inf
  for (i in seq_len(n_try)) {
    t0 <- tic()
    pars <- list(
      max_depth        = sample(3:8, 1),
      eta              = runif(1, 0.02, 0.3),
      subsample        = runif(1, 0.7, 1.0),
      colsample_bytree = runif(1, 0.7, 1.0),
      min_child_weight = sample(1:10, 1),
      tree_method      = "hist",
      max_bin          = 256,
      nthread          = 4,
      objective        = "reg:squarederror",
      eval_metric      = "rmse",
      verbosity        = 0
    )
    nrounds_cap <- sample(nrounds_grid, 1)
    cv <- xgb.cv(
      params = pars,
      data = dtrain,
      nrounds = nrounds_cap,
      nfold = nfold_cv,
      verbose = 0,
      early_stopping_rounds = es_rounds,
      showsd = TRUE,
      stratified = FALSE,
      maximize = FALSE,
      seed = seed + i
    )
    best_it <- cv$best_iteration
    rmse <- cv$evaluation_log$test_rmse_mean[best_it]
    if (is.finite(rmse) && rmse < best_rmse) {
      best_rmse <- rmse
      best <- c(pars, list(nrounds = as.integer(best_it)))
    }
    toc(t0, sprintf("XGBoost(regression) random search [%d] done (best_it=%d, rmse=%.6f)", i, best_it, rmse))
  }
  list(best = best, rmse = best_rmse)
}

# ========= Build exposure (A) under the requested lag specification =========
make_X_lag <- function(DT, pol, mode) {
  if (mode$type == "single") {
    col <- sprintf("%s_lag%d", pol, mode$k)
    stopifnot(col %in% names(DT))
    return(DT[[col]])
  } else {
    cols <- sprintf("%s_lag%d", pol, 0:mode$k)
    stopifnot(all(cols %in% names(DT)))
    return(rowMeans(DT[, ..cols], na.rm = TRUE))
  }
}

# ========= Build confounder/features (T) consistent with lag specification =========
make_T <- function(DT, pol, mode) {
  others <- setdiff(pollutants, pol)
  other_cols <- if (mode$type == "single") {
    unlist(lapply(others, function(p) sprintf("%s_lag%d", p, mode$k)))
  } else {
    others
  }
  
  if (length(other_cols)) {
    if (mode$type == "single") {
      Xothers <- DT[, ..other_cols]
    } else {
      pollutant_means <- lapply(others, function(pollutant) {
        poll_cols <- sprintf("%s_lag%d", pollutant, 0:mode$k)
        stopifnot(all(poll_cols %in% names(DT)))
        rowMeans(DT[, ..poll_cols], na.rm = TRUE)
      })
      Xothers <- as.data.table(pollutant_means)
      setnames(Xothers, others)
    }
  } else {
    Xothers <- NULL
  }
  
  temp_lag03 <- rowMeans(DT[, sprintf("temp_lag%d", 0:3), with = FALSE], na.rm = TRUE)
  rh_lag03   <- rowMeans(DT[, sprintf("rh_lag%d",   0:3), with = FALSE], na.rm = TRUE)
  
  base_cols <- c("age","gender","nation","edu","work","marriage","holiday","city_id")
  baseDT <- DT[, ..base_cols]
  
  out <- cbind(
    Xothers,
    temp_lag03 = temp_lag03,
    rh_lag03   = rh_lag03,
    baseDT
  )
  as.data.frame(out)
}

# ========= Conditional-logit probability within each stratum (softmax) =========
softmax_by_strata <- function(eta, strata_vec) {
  p <- numeric(length(eta))
  uu <- split(seq_along(eta), strata_vec)
  for (idx in uu) {
    z <- eta[idx] - max(eta[idx], na.rm = TRUE)
    ez <- exp(z); s <- sum(ez, na.rm = TRUE)
    p[idx] <- if (is.finite(s) && s > 0) ez / s else 1 / length(idx)
  }
  p[!is.finite(p)] <- mean(p[is.finite(p)], na.rm = TRUE)
  p
}

# ========= Outcome model: conditional logistic regression used in AIPW =========
fit_outcome_clogit_AIPW <- function(DT, T_df, A_vec) {
  other_names <- setdiff(colnames(T_df),
                         c("temp_lag03","rh_lag03","age","gender","nation","edu","work","marriage","holiday","city_id"))
  Aname <- "A_expo"
  fml <- as.formula(
    paste0("case ~ ", Aname,
           if (length(other_names)) paste0(" + ", paste(other_names, collapse=" + ")) else "",
           " + ns(temp_lag03, df=6)",
           " + ns(rh_lag03, df=3)",
           " + holiday + strata(id)")
  )
  DF  <- cbind(DT[, .(case, id)], T_df, A_expo = A_vec)
  fit <- clogit(fml, data = DF, method = "efron", x = TRUE)
  Xmat <- fit$x
  beta <- coef(fit)
  eta  <- drop(Xmat %*% beta)
  phat <- softmax_by_strata(eta, strata_vec = DF$id)
  beta_A <- as.numeric(beta["A_expo"])
  if (!is.finite(beta_A)) stop("Coefficient for A_expo not found (beta_A). Make sure A_expo enters linearly.")
  list(fit=fit, Xmat=Xmat, beta=beta, eta=eta, phat=phat, beta_A=beta_A)
}

# ========= Case-day-only shift helper =========
get_case_rows_first <- function(ids, case_vec) {
  dtc <- data.table(idx = which(case_vec == 1), id = ids[case_vec == 1])
  if (nrow(dtc) == 0L) stop("No case==1 days were found.")
  dtc <- dtc[!duplicated(id)]
  dtc$idx
}

# ========= AIPW with shifting exposure only on the case day within each stratum =========
aipw_shift_caseonly <- function(mA, resY, g, sigma2, A, delta, ids, case_vec,
                                beta_A, lr_clip = c(1/50, 50)) {
  idx_case <- get_case_rows_first(ids, case_vec)
  id_case  <- ids[idx_case]
  S <- length(idx_case)
  
  mult <- exp(beta_A * delta)
  m0   <- mA[idx_case]
  denom_factor <- 1 + m0 * (mult - 1)
  m_shift <- (m0 * mult) / pmax(denom_factor, 1e-300)
  
  sdv <- sqrt(pmax(sigma2[idx_case], 1e-12))
  num <- dnorm(A[idx_case] - delta, mean = g[idx_case], sd = sdv)
  den <- dnorm(A[idx_case],         mean = g[idx_case], sd = sdv)
  LR_case <- pmax(pmin(num / pmax(den, 1e-300), lr_clip[2]), lr_clip[1])
  
  mu0_hat <- pmax(m0, 1e-12)
  mu1_hat <- m_shift + LR_case * resY[idx_case]
  mu1_hat <- pmax(mu1_hat, 1e-12)
  
  rho_hat <- mu1_hat / mu0_hat
  log_rho <- log(pmax(rho_hat, 1e-300))
  
  psi_hat <- mean(log_rho, na.rm = TRUE)
  IF_s <- log_rho - psi_hat
  se_psi <- sqrt(sum(IF_s^2, na.rm = TRUE)) / S
  
  LR_full <- rep(NA_real_, length(ids)); LR_full[idx_case] <- LR_case
  IF_full <- rep(NA_real_, length(ids)); IF_full[idx_case] <- IF_s
  
  list(
    theta = psi_hat,
    IF    = IF_full,
    se_cl = se_psi,
    LR    = LR_full,
    idx_case = idx_case,
    id_case  = id_case,
    m0_case  = m0,
    m_shift_case = m_shift,
    mu1_case = mu1_hat,
    rho_case = rho_hat,
    log_rho_case = log_rho,
    mult = mult,
    denom_factor = denom_factor
  )
}

# ========= Tuning subsample by (year, city) composition =========
sample_tune_ids <- function(DT, n_id_target, seed = SEED_MASTER + 2L) {
  setDT(DT)
  keycols <- c("year", "city")
  tmp <- unique(DT[, .(id, year, city)])
  group_counts <- tmp[!is.na(year) & !is.na(city), .N, by = .(year, city)]
  total_ids <- nrow(tmp)
  n_id_target <- as.numeric(n_id_target)
  group_counts[, target_n := as.integer(round((N / total_ids) * n_id_target))]
  current_total <- sum(group_counts$target_n)
  
  if (current_total != n_id_target) {
    diff <- n_id_target - current_total
    if (diff > 0) {
      group_counts <- group_counts[order(-N)]
      for (i in 1:abs(diff)) group_counts[i, target_n := target_n + 1]
    } else {
      group_counts <- group_counts[order(N)]
      for (i in 1:abs(diff)) {
        idx <- which(group_counts$target_n > 1)[1]
        if (!is.na(idx)) group_counts[idx, target_n := target_n - 1]
      }
    }
  }
  group_counts[target_n == 0 & N > 0, target_n := 1]
  final_total <- sum(group_counts$target_n)
  
  set.seed(seed)
  if (final_total != n_id_target) {
    diff <- n_id_target - final_total
    if (diff > 0) {
      add_groups <- sample(1:nrow(group_counts), diff, replace = TRUE)
      for (grp in add_groups) group_counts[grp, target_n := target_n + 1]
    } else {
      eligible_groups <- which(group_counts$target_n > 1)
      if (length(eligible_groups) >= abs(diff)) {
        remove_groups <- sample(eligible_groups, abs(diff))
        for (grp in remove_groups) group_counts[grp, target_n := target_n - 1]
      }
    }
  }
  
  tmp2 <- tmp[group_counts, on = keycols][,
                                          .(id = id[sample(.N, min(.N, target_n))]), by = keycols
  ]
  
  result_ids <- tmp2$id
  if (length(result_ids) > n_id_target) {
    result_ids <- sample(result_ids, n_id_target)
  } else if (length(result_ids) < n_id_target) {
    remaining_ids <- setdiff(tmp$id, result_ids)
    if (length(remaining_ids) > 0) {
      need_more <- n_id_target - length(result_ids)
      additional_ids <- sample(remaining_ids, min(need_more, length(remaining_ids)))
      result_ids <- c(result_ids, additional_ids)
    }
  }
  unique(result_ids)
}

# ========= Correlation + distribution summaries (for documentation) =========
make_pollutant_matrix_for_mode <- function(DT, pols, mode) {
  if (mode$type == "single") {
    cols <- sapply(pols, function(p) sprintf("%s_lag%d", p, mode$k))
    stopifnot(all(cols %in% names(DT)))
    X <- as.data.frame(DT[, ..cols]); names(X) <- pols; return(X)
  } else {
    Xlist <- lapply(pols, function(p) {
      cns <- sprintf("%s_lag%d", p, 0:mode$k)
      stopifnot(all(cns %in% names(DT)))
      rowMeans(DT[, ..cns], na.rm = TRUE)
    })
    X <- as.data.frame(Xlist); names(X) <- pols; return(X)
  }
}

compute_and_save_corr <- function(DT, pols, mode, out_dir, lag_tag) {
  X <- make_pollutant_matrix_for_mode(DT, pols, mode)
  C <- cor(X, use = "pairwise.complete.obs", method = "pearson")
  cor_dir <- file.path(out_dir, "corr_matrices")
  dir.create(cor_dir, showWarnings = FALSE, recursive = TRUE)
  cor_df <- cbind(pollutant = rownames(C), as.data.frame(C))
  fn_csv <- file.path(cor_dir, sprintf("corr_%s.csv", lag_tag))
  data.table::fwrite(cor_df, fn_csv)
  message(sprintf("[Correlation] Written: %s", fn_csv))
  invisible(fn_csv)
}

compute_and_save_stats <- function(DT, pols, mode, out_dir, lag_tag,
                                   probs = c(.01,.05,.10,.25,.50,.75,.90,.95,.99)) {
  X <- make_pollutant_matrix_for_mode(DT, pols, mode)
  get_row <- function(v) {
    qs <- as.numeric(quantile(v, probs = probs, na.rm = TRUE))
    data.table(
      n      = sum(is.finite(v)),
      na     = sum(!is.finite(v)),
      mean   = mean(v, na.rm = TRUE),
      sd     = sd(v, na.rm = TRUE),
      median = median(v, na.rm = TRUE),
      min    = suppressWarnings(min(v, na.rm = TRUE)),
      max    = suppressWarnings(max(v, na.rm = TRUE)),
      p01 = qs[1], p05 = qs[2], p10 = qs[3], p25 = qs[4],
      p50 = qs[5], p75 = qs[6], p90 = qs[7], p95 = qs[8], p99 = qs[9],
      IQR    = IQR(v, na.rm = TRUE)
    )
  }
  statsDT <- rbindlist(lapply(colnames(X), function(col) {
    cbind(data.table(pollutant = col), get_row(X[[col]]))
  }), use.names = TRUE, fill = TRUE)
  dir_stats <- file.path(out_dir, "stats_tables")
  dir.create(dir_stats, showWarnings = FALSE, recursive = TRUE)
  fn_csv <- file.path(dir_stats, sprintf("stats_%s.csv", lag_tag))
  fwrite(statsDT, fn_csv)
  message(sprintf("[Summary stats] Written: %s", fn_csv))
  invisible(fn_csv)
}

# ========= 3) Main loop over pollutants =========
all_results <- list()
t_all0 <- tic()

for (pol in pollutants2) {
  t_pol0 <- tic()
  
  DELTA_pol <- if (!is.null(DELTA_MAP[[pol]])) as.numeric(DELTA_MAP[[pol]]) else DELTA
  vf_floor_pol <- if (!is.null(VF_FLOOR_MAP[[pol]])) as.numeric(VF_FLOOR_MAP[[pol]]) else 0.05
  vf_cap_pol   <- if (!is.null(VF_CAP_MAP[[pol]]))   as.numeric(VF_CAP_MAP[[pol]])   else 0.95
  
  mode <- LAG_PLAN[[pol]]
  stopifnot(is.list(mode), !is.null(mode$type), !is.null(mode$k), mode$k >= 0)
  lag_tag <- lag_tag_from_mode(mode)
  
  pol_dir <- file.path(OUT_DIR, pol)
  dir.create(pol_dir, showWarnings = FALSE, recursive = TRUE)
  message(sprintf(">>> Start tuning/estimation (AIPW ): %s | %s", pol, lag_tag))
  
  # --- 3.1) Hyperparameter tuning on a smaller, stratified id subset ---
  mode_tune <- mode
  t0 <- tic()
  tune_ids <- sample_tune_ids(data, tune_sample_n_id, seed = SEED_MASTER + 10L)
  DT_tune  <- copy(data[id %in% tune_ids])
  toc(t0, sprintf("Tuning subset sampled: n_id=%d, n=%d", length(unique(DT_tune$id)), nrow(DT_tune)))
  
  X_tune_raw <- make_X_lag(DT_tune, pol, mode_tune)
  T_tune_raw <- make_T(DT_tune, pol, mode_tune)
  T_tune_DT  <- as.data.table(T_tune_raw)
  
  # Demean within matched sets (control-day mean) for: co-pollutants + temp/rh + exposure
  t0 <- tic()
  if (mode_tune$type == "single") {
    others <- setdiff(pollutants, pol)
    if (pol == "PM2510") others <- setdiff(others, "PM10")
    pol_cols_in_T <- intersect(colnames(T_tune_DT), sprintf("%s_lag%d", others, mode_tune$k))
  } else {
    pol_cols_in_T <- intersect(colnames(T_tune_DT), setdiff(pollutants, pol))
    if (pol == "PM2510") pol_cols_in_T <- setdiff(pol_cols_in_T, "PM10")
  }
  
  DT_tune[, X_tune_tmp := X_tune_raw]
  tmp_in <- cbind(
    DT_tune[, .(id, case)],
    if (length(pol_cols_in_T)) T_tune_DT[, ..pol_cols_in_T] else NULL,
    T_tune_DT[, .(temp_lag03, rh_lag03)],
    X_tune_tmp = DT_tune$X_tune_tmp
  )
  dm_cols    <- c(pol_cols_in_T, "temp_lag03", "rh_lag03", "X_tune_tmp")
  tmpdm      <- demean_by_control(tmp_in, dm_cols)
  T_tune_dm  <- as.data.frame(T_tune_DT)
  if (length(pol_cols_in_T)) T_tune_dm[, pol_cols_in_T] <- as.data.frame(tmpdm[, ..pol_cols_in_T])
  T_tune_dm[, c("temp_lag03","rh_lag03")] <- as.data.frame(tmpdm[, .(temp_lag03, rh_lag03)])
  X_tune_dm <- tmpdm$X_tune_tmp
  DT_tune[, X_tune_tmp := NULL]
  toc(t0, "Tuning stage: within-stratum demeaning done (consistent with lag plan)")
  
  print("Stage-1 exposure model tuning (XGBoost regression) starts")
  t0 <- tic()
  dtrain_exp_tune <- xgb.DMatrix(data = data.matrix(as.data.frame(T_tune_dm)), label = X_tune_dm)
  xgb_tune_exp <- xgb_random_search_reg(dtrain_exp_tune, n_try = tune_try_random, seed = SEED_MASTER + 100L)
  saveRDS(xgb_tune_exp, file.path(pol_dir, paste0("xgbEXP_tune_", lag_tag, ".rds")))
  rm(dtrain_exp_tune); gc()
  toc(t0, "Stage-1 exposure model tuning finished")
  
  # --- 3.2) K-fold cross-fitting by id ---
  id2fold <- make_folds(unique(data$id), K_fold, seed = SEED_MASTER + 1L)
  data <- merge(data, id2fold, by = "id", all.x = TRUE)
  
  t_lag0 <- tic()
  message(sprintf("---- %s | %s ----", pol, lag_tag))
  
  if (identical(pol, pollutants[1])) {
    compute_and_save_corr(DT = data, pols = pollutants, mode = mode, out_dir = OUT_DIR, lag_tag = lag_tag)
    compute_and_save_stats(DT = data, pols = pollutants, mode = mode, out_dir = OUT_DIR, lag_tag = lag_tag)
  }
  
  # --- 3.2.1) Construct A and T (raw scale for outcome model) ---
  X_raw <- make_X_lag(data, pol, mode)
  T_raw <- make_T(data, pol, mode)
  T_DT  <- as.data.table(T_raw)
  
  # --- 3.2.2) Outcome model m(A,T,id) via conditional logistic regression ---
  t0 <- tic()
  out_fit <- fit_outcome_clogit_AIPW(data, T_raw, A_vec = X_raw)
  mA      <- out_fit$phat
  beta_A  <- out_fit$beta_A
  eta     <- out_fit$eta
  U       <- data$case - mA
  toc(t0, "Outcome model fitted (mA, residual U, beta_A)")
  
  # --- Diagnostic: residual association check for outcome model ---
  {
    num_cols_U <- which(sapply(T_raw, is.numeric))
    if (length(num_cols_U) > 0) {
      Tmat_U <- as.matrix(T_raw[, num_cols_U, drop = FALSE])
      cor_vec_U <- apply(Tmat_U, 2, function(col) suppressWarnings(cor(col, U, use="pairwise.complete.obs")))
      mean_abs_cor_U_T <- mean(abs(cor_vec_U[is.finite(cor_vec_U)]), na.rm = TRUE)
    } else mean_abs_cor_U_T <- NA_real_
    U_fit_ok <- is.finite(mean_abs_cor_U_T) && mean_abs_cor_U_T < 0.1
    message(sprintf("[Diag] Outcome residual AC = %.3f | (<0.10 preferred) = %s",
                    mean_abs_cor_U_T, ifelse(U_fit_ok,"TRUE","FALSE")))
  }
  
  # --- 3.2.3) Demean for exposure model g(T) and variance model sigma^2(T) ---
  t0 <- tic()
  if (mode$type == "single") {
    others2 <- setdiff(pollutants, pol)
    if (pol == "PM2510") others2 <- setdiff(others2, "PM10")
    pol_cols_in_T2 <- intersect(colnames(T_DT), sprintf("%s_lag%d", others2, mode$k))
  } else {
    pol_cols_in_T2 <- intersect(colnames(T_DT), setdiff(pollutants, pol))
    if (pol == "PM2510") pol_cols_in_T2 <- setdiff(pol_cols_in_T2, "PM10")
  }
  
  data[, X_tmp_dm := X_raw]
  tmp_in <- cbind(
    data[, .(id, case)],
    if (length(pol_cols_in_T2)) T_DT[, ..pol_cols_in_T2] else NULL,
    T_DT[, .(temp_lag03, rh_lag03)],
    X_tmp_dm = data$X_tmp_dm
  )
  dm_cols <- c(pol_cols_in_T2, "temp_lag03", "rh_lag03", "X_tmp_dm")
  tmpdm   <- demean_by_control(tmp_in, dm_cols)
  T_dm <- as.data.frame(T_DT)
  if (length(pol_cols_in_T2)) T_dm[, pol_cols_in_T2] <- as.data.frame(tmpdm[, ..pol_cols_in_T2])
  T_dm[, c("temp_lag03","rh_lag03")] <- as.data.frame(tmpdm[, .(temp_lag03, rh_lag03)])
  X_dm <- tmpdm$X_tmp_dm
  data[, X_tmp_dm := NULL]
  toc(t0, "Within-stratum demeaning finished (co-pollutants + temp/rh + exposure)")
  
  setDT(T_dm)
  MM <- data.matrix(as.data.frame(T_dm))
  stopifnot(nrow(MM) == nrow(data))
  
  # --- 3.2.4) Cross-fitted ghat(T) via XGBoost regression ---
  ehat <- rep(NA_real_, nrow(data))
  rmse_tr_vec <- c(); rmse_te_vec <- c()
  r2_tr_vec   <- c(); r2_te_vec   <- c()
  
  r2_score <- function(y, yhat) {
    sse <- sum((y - yhat)^2, na.rm = TRUE)
    sst <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
    if (!is.finite(sst) || sst <= 0) return(NA_real_)
    1 - sse/sst
  }
  rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))
  
  t0 <- tic(); print("Stage-1 exposure model (XGBoost) cross-fitting starts")
  for (k in 1:K_fold) {
    t1 <- tic()
    tr_idx <- which(data$fold != k)
    te_idx <- which(data$fold == k)
    
    set.seed(SEED_MASTER + 1000L + k)
    dtr <- xgb.DMatrix(data = MM[tr_idx, , drop = FALSE], label = X_dm[tr_idx])
    bst <- xgb.train(params = xgb_tune_exp$best, data = dtr,
                     nrounds = xgb_tune_exp$best$nrounds, verbose = 0)
    
    pred_tr <- predict_xgb_in_chunks(bst, MM, tr_idx, chunk_rows = 200000L)
    pred_te <- predict_xgb_in_chunks(bst, MM, te_idx, chunk_rows = 200000L)
    
    rmse_tr_vec <- c(rmse_tr_vec, rmse(X_dm[tr_idx], pred_tr))
    rmse_te_vec <- c(rmse_te_vec, rmse(X_dm[te_idx], pred_te))
    r2_tr_vec   <- c(r2_tr_vec,   r2_score(X_dm[tr_idx], pred_tr))
    r2_te_vec   <- c(r2_te_vec,   r2_score(X_dm[te_idx], pred_te))
    ehat[te_idx] <- pred_te
    
    rm(dtr, bst, pred_tr, pred_te); gc(FALSE)
    toc(t1, sprintf("Stage-1 exposure model fold [%d] done", k))
  }
  ghat <- ehat
  A_resid <- X_dm - ghat
  
  # --- 3.2.5) Cross-fitted variance model sigma^2(T) on log(residual^2) ---
  t0 <- tic()
  sigma2_hat <- rep(NA_real_, nrow(data))
  for (k in 1:K_fold) {
    tr_idx <- which(data$fold != k)
    te_idx <- which(data$fold == k)
    ytr_v  <- log(pmax(A_resid[tr_idx]^2, 1e-12))
    
    set.seed(SEED_MASTER + 2000L + k)
    dtr_v <- xgb.DMatrix(MM[tr_idx, , drop = FALSE], label = ytr_v)
    bst_v <- xgb.train(params = xgb_tune_exp$best,
                       data = dtr_v,
                       nrounds = round(xgb_tune_exp$best$nrounds/2), verbose = 0)
    v_te <- exp(pmin(predict_xgb_in_chunks(bst_v, MM, te_idx, chunk_rows = 200000L), 20))
    sigma2_hat[te_idx] <- v_te
    
    rm(dtr_v, bst_v, v_te, ytr_v); gc(FALSE)
  }
  toc(t0, "Variance model cross-fitting finished")
  
  # --- Variance calibration: floor/cap based on residual^2 quantiles + inflation factor ---
  res2   <- (X_dm - ghat)^2
  floor2 <- quantile(res2, vf_floor_pol, na.rm=TRUE)
  cap2   <- quantile(res2, vf_cap_pol,   na.rm=TRUE)
  sigma2_hat <- pmax(pmin(sigma2_hat, cap2), floor2)
  infl <- mean(res2, na.rm=TRUE) / mean(sigma2_hat, na.rm=TRUE)
  sigma2_hat <- pmax(sigma2_hat * infl, 1e-12)
  message(sprintf("[Variance calib] vf_floor=%.3f (floor=%.4g), vf_cap=%.3f (cap=%.4g), infl=%.3f",
                  vf_floor_pol, floor2, vf_cap_pol, cap2, infl))
  
  # --- Overfit report for exposure model ---
  rmse_tr_mean <- mean(rmse_tr_vec, na.rm = TRUE)
  rmse_te_mean <- mean(rmse_te_vec, na.rm = TRUE)
  r2_tr_mean   <- mean(r2_tr_vec,   na.rm = TRUE)
  r2_te_mean   <- mean(r2_te_vec,   na.rm = TRUE)
  rmse_change_pct <- 100 * (rmse_te_mean - rmse_tr_mean) / max(1e-12, rmse_tr_mean)
  r2_change_pct   <- 100 * (r2_te_mean   - r2_tr_mean)   / max(1e-12, abs(r2_tr_mean))
  toc(t0, sprintf(
    "Stage-1 exposure model done | RMSE: train=%.4f, test=%.4f, Δ=%.2f%% | R2: train=%.3f, test=%.3f, Δ=%.2f%%",
    rmse_tr_mean, rmse_te_mean, rmse_change_pct, r2_tr_mean, r2_te_mean, r2_change_pct
  ))
  
  # --- Diagnostic: residual association check for exposure model ---
  {
    num_cols_A <- which(sapply(T_dm, is.numeric))
    if (length(num_cols_A) > 0) {
      Tmat_A <- as.matrix(as.data.frame(T_dm)[, num_cols_A, drop = FALSE])
      cor_vec_A <- apply(Tmat_A, 2, function(col) suppressWarnings(cor(col, A_resid, use="pairwise.complete.obs")))
      mean_abs_cor_A_T <- mean(abs(cor_vec_A[is.finite(cor_vec_A)]), na.rm = TRUE)
    } else mean_abs_cor_A_T <- NA_real_
    A_fit_ok <- is.finite(mean_abs_cor_A_T) && mean_abs_cor_A_T < 0.1
    message(sprintf("[Diag] Exposure residual AC = %.3f | (<0.10 preferred) = %s",
                    mean_abs_cor_A_T, ifelse(A_fit_ok,"TRUE","FALSE")))
  }
  
  # --- 3.2.6) AIPW: case-day-only shift, reporting per-DELTA and per-IQR effects ---
  t0 <- tic()
  
  aipw1 <- aipw_shift_caseonly(
    mA = mA,
    resY = U,
    g = ghat,
    sigma2 = sigma2_hat,
    A = X_dm,
    delta = DELTA_pol,
    ids = data$id,
    case_vec = data$case,
    beta_A = beta_A,
    lr_clip = LR_CLIP
  )
  logRR_delta  <- aipw1$theta
  se_logRR_del <- aipw1$se_cl
  
  logRR_perDelta    <- logRR_delta / DELTA_pol
  se_logRR_perDelta <- se_logRR_del / DELTA_pol
  
  Zval <- logRR_perDelta / se_logRR_perDelta
  Pval <- 2 * pnorm(-abs(Zval))
  CI95_logRR_delta <- c(logRR_perDelta - 1.96 * se_logRR_perDelta,
                        logRR_perDelta + 1.96 * se_logRR_perDelta)
  
  RR_perDelta   <- exp(logRR_perDelta)
  CI95_RR_delta <- exp(CI95_logRR_delta)
  
  IQR_x <- IQR(X_raw, na.rm = TRUE)
  aipwI <- aipw_shift_caseonly(
    mA = mA,
    resY = U,
    g = ghat,
    sigma2 = sigma2_hat,
    A = X_dm,
    delta = IQR_x,
    ids = data$id,
    case_vec = data$case,
    beta_A = beta_A,
    lr_clip = LR_CLIP
  )
  IQR_logRR <- aipwI$theta
  SE_IQR    <- aipwI$se_cl
  CI95_IQR_logRR <- c(IQR_logRR - 1.96 * SE_IQR,
                      IQR_logRR + 1.96 * SE_IQR)
  
  RR_IQR      <- exp(IQR_logRR)
  CI95_RR_IQR <- exp(CI95_IQR_logRR)
  
  toc(t0, "AIPW(case-day-only shift) finished: logRR(per-DELTA) + logRR(IQR) + SE/CI/p")
  
  message(sprintf(
    paste0("[Result] %s|%s| logRR(per Δ=%.3g)=%.4e (SE_cl=%.4e, p=%.3g, CI=[%.4e, %.4e])",
           " | RR(per Δ)=%.4f [%.4f, %.4f]",
           " | IQR=%.4f | IQR_logRR=%.4e [%.4e, %.4e] | IQR_RR=%.4f [%.4f, %.4f]",
           " | AC_Y=%.3f | AC_A=%.3f | ΔRMSE=%.1f%% ΔR2=%.1f%%"),
    pol, lag_tag, DELTA_pol,
    logRR_perDelta, se_logRR_perDelta, Pval, CI95_logRR_delta[1], CI95_logRR_delta[2],
    RR_perDelta, CI95_RR_delta[1], CI95_RR_delta[2],
    IQR_x, IQR_logRR, CI95_IQR_logRR[1], CI95_IQR_logRR[2],
    RR_IQR, CI95_RR_IQR[1], CI95_RR_IQR[2],
    mean_abs_cor_U_T, mean_abs_cor_A_T, rmse_change_pct, r2_change_pct
  ))
  
  # Construct mA_pD for compatibility: normalize probabilities after shifting only the case day
  idx_case <- aipw1$idx_case
  id_case  <- aipw1$id_case
  denom_factor_case <- aipw1$denom_factor
  mult1 <- aipw1$mult
  
  pos_map <- match(data$id, id_case)
  denom_factor_all <- denom_factor_case[pos_map]
  denom_factor_all[!is.finite(denom_factor_all) | denom_factor_all <= 0] <- 1
  
  mApD <- mA / denom_factor_all
  mApD[idx_case] <- mApD[idx_case] * mult1
  
  saveRDS(list(
    id      = data$id,
    case    = data$case,
    A_raw   = X_raw,
    A_dm    = X_dm,
    ghat    = ghat,
    sigma2  = sigma2_hat,
    mA      = mA,
    mA_pD   = mApD,
    U       = U,
    LR      = aipw1$LR,
    theta_RD = logRR_delta,
    IF_RD    = aipw1$IF,
    u0       = NA_real_,
    u1       = NA_real_,
    logRR_perDelta = logRR_perDelta,
    SE_logRR       = se_logRR_perDelta,
    CI95_logRR     = CI95_logRR_delta,
    RR_perDelta    = RR_perDelta,
    CI95_RR        = CI95_RR_delta,
    IQR_x          = IQR_x,
    IQR_logRR      = IQR_logRR,
    SE_IQR         = SE_IQR,
    CI95_IQR_logRR = CI95_IQR_logRR,
    RR_IQR         = RR_IQR,
    CI95_RR_IQR    = CI95_RR_IQR
  ), file = file.path(pol_dir, sprintf("%s_%s_DRshift.rds", pol, lag_tag)))
  
  all_results[[length(all_results)+1]] <- data.table(
    pollutant      = pol,
    lag_mode       = lag_tag,
    is_cum         = as.integer(mode$type=="cum"),
    DELTA          = DELTA_pol,
    logRR_perDelta = logRR_perDelta,
    logRR_CI_L     = CI95_logRR_delta[1],
    logRR_CI_U     = CI95_logRR_delta[2],
    RR_perDelta    = RR_perDelta,
    RR_CI_L        = CI95_RR_delta[1],
    RR_CI_U        = CI95_RR_delta[2],
    IQR            = IQR_x,
    IQR_logRR      = IQR_logRR,
    IQR_logRR_CI_L = CI95_IQR_logRR[1],
    IQR_logRR_CI_U = CI95_IQR_logRR[2],
    IQR_RR         = RR_IQR,
    IQR_RR_CI_L    = CI95_RR_IQR[1],
    IQR_RR_CI_U    = CI95_RR_IQR[2],
    SE_logRR       = se_logRR_perDelta,
    Z              = Zval,
    P              = Pval,
    SE_IQR         = SE_IQR,
    AC_Y           = mean_abs_cor_U_T,
    AC_A           = mean_abs_cor_A_T,
    RMSE_train     = rmse_tr_mean,
    RMSE_test      = rmse_te_mean,
    R2_train       = r2_tr_mean,
    R2_test        = r2_te_mean,
    dRMSE_pct      = rmse_change_pct,
    dR2_pct        = r2_change_pct,
    LR_P1          = as.numeric(quantile(aipw1$LR, .01, na.rm=TRUE)),
    LR_P50         = as.numeric(quantile(aipw1$LR, .50, na.rm=TRUE)),
    LR_P95         = as.numeric(quantile(aipw1$LR, .95, na.rm=TRUE)),
    LR_P99         = as.numeric(quantile(aipw1$LR, .99, na.rm=TRUE))
  )
  
  toc(t_lag0, sprintf("%s | %s completed", pol, lag_tag))
  toc(t_pol0, sprintf("Pollutant %s finished under the fixed lag plan", pol))
  
  data[, fold := NULL]
  rm(T_DT, T_dm, T_raw, X_raw, MM); gc(FALSE)
}

# ========= 4) Save summary table =========
stopifnot(exists("all_results"), length(all_results) > 0)
resDT <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
resDT[, P_BH := p.adjust(P, method = "BH"), by = .(pollutant)]
resDT[, `:=`(
  sig_BH_05 = as.logical(P_BH < 0.05),
  sig_BH_10 = as.logical(P_BH < 0.10),
  rank_P    = frank(P,    ties.method = "min"),
  rank_BH   = frank(P_BH, ties.method = "min")
), by = .(pollutant)]
setorder(resDT, pollutant, is_cum, lag_mode)
f_csv <- file.path(OUT_DIR, "DRshift_summary.csv")
fwrite(resDT, f_csv)
message(sprintf("Summary with BH-adjusted p-values saved: %s", f_csv))

toc(t_all0, "All pollutants completed (DR/shift AIPW, case-day-only shift)")

# ========= 5) Session info + run metadata =========
SESSION_INFO <- capture.output({
  si <- sessionInfo(); print(si)
})
writeLines(SESSION_INFO, con = file.path(OUT_DIR, "SESSION_INFO.txt"))
message("SESSION_INFO.txt written (path-free content)")

run_summary <- list(
  runtime      = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  params = list(
    K_fold          = K_fold,
    tune_try_random = tune_try_random,
    DELTA_default   = DELTA,
    LR_CLIP         = as.numeric(LR_CLIP),
    PIT_BINS        = PIT_BINS,
    VF_FLOOR_MAP    = as.list(VF_FLOOR_MAP),
    VF_CAP_MAP      = as.list(VF_CAP_MAP),
    DELTA_MAP       = as.list(DELTA_MAP),
    pollutants      = pollutants,
    pollutants2     = pollutants2,
    LAG_PLAN        = LAG_PLAN,
    shift_type      = "case-day-only"
  ),
  packages = list(
    R_version   = R.version$version.string,
    platform    = R.version$platform,
    xgboost     = as.character(packageVersion("xgboost")),
    data.table  = as.character(packageVersion("data.table")),
    survival    = as.character(packageVersion("survival")),
    splines     = "base"
  ),
  outputs = list(
    summary_csv           = "DRshift_summary.csv",
    per_lag_rds_template  = "<pollutant>_<lag>_DRshift.rds"
  )
)
jsonlite::write_json(
  run_summary,
  path = file.path(OUT_DIR, "RUN_SUMMARY.json"),
  pretty = TRUE, auto_unbox = TRUE
)
message("RUN_SUMMARY.json written (filenames only, no absolute paths)")

sink()
