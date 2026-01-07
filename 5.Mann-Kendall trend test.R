##########################################################################################################
#Codes for "Reduced mortality burden associated with the clean air policy in China: An integrated assessment based on causal machine learning"
#Authors for codes: Yuchang Zhou, Peizheng Li, Ya Gao, Yi Guo, Yixiang Zhu, Lu Zhou, Peng Yin, Haidong Kan, Maigeng Zhou, Renjie Chen
#Correspondence to Renjie Chen, Maigeng Zhou.
###########################################################################################################

##########################################################################################################
# Notes on required inputs and expected outputs
#
# This script performs trend analysis on annual attributable fractions (AF%) for each pollutant and model.
# It is a post-processing step and does not fit exposure–response models. It reads annual AF tables
# produced upstream (typically one from AIPW and one from conditional logistic regression), then:
#   (1) extracts AF% time series for 2013–2019,
#   (2) computes endpoint-based declines (2013 vs 2019) with endpoint CI substitution,
#   (3) applies robust trend tests (Theil–Sen slope + Mann–Kendall test).
#
# ------------------------------
# 1) REQUIRED INPUT DATA (FROM UPSTREAM SCRIPTS)
# ------------------------------
# The script expects two annual AF wide tables (CSV), one per model, with the following structure:
#
# (A) AIPW_FILE (annual_AF_all_pollutants.csv under AIPW_DIR)
# (B) CLOGIT_FILE (annual_AF_all_pollutants.csv under CLOGIT_DIR)
#
# Each file must contain:
#   - year: calendar year as integer-like (e.g., 2013–2019). A row "Total" is allowed but will be dropped.
#   - For each pollutant p in pollutants_all (depending on INCLUDE_SO2 / INCLUDE_CO):
#       p        : annual AF% (midpoint estimate)
#       p_LCL    : annual AF% lower CI endpoint (optional; if absent will be filled as NA)
#       p_UCL    : annual AF% upper CI endpoint (optional; if absent will be filled as NA)
#
# Pollutant column names must match exactly, e.g.:
#   PM25, PM25_LCL, PM25_UCL,
#   O3, O3_LCL, O3_UCL, etc.
#
# ------------------------------
# 2) METHODS IMPLEMENTED
# ------------------------------
# For each (model, pollutant) time series over 2013–2019:
#   - Endpoint absolute reduction (percentage points, pp):
#       abs_reduction_pp = AF_2013 - AF_2019
#     with CI computed by endpoint substitution:
#       abs_LCL = AF_2013_LCL - AF_2019_LCL
#       abs_UCL = AF_2013_UCL - AF_2019_UCL
#   - Endpoint relative reduction:
#       rel_reduction = (AF_2013 - AF_2019) / AF_2013
#     reported also as percentage (×100), with endpoint-substitution CI.
#   - Robust trend estimation:
#       * Theil–Sen slope (pp per year): trend::sens.slope (fallback: zyp::zyp.sen)
#       * Mann–Kendall trend test: trend::mk.test (fallback: Kendall::MannKendall)
#
# ------------------------------
# 3) OUTPUTS
# ------------------------------
# The script writes one summary CSV:
#   OUT_CSV = annual_AF_decline_summary_2013_2019_endpointCI.csv
# containing, for each (model, pollutant):
#   - AF_2013_mid/LCL/UCL and AF_2019_mid/LCL/UCL
#   - absolute reduction (pp) mid/LCL/UCL
#   - relative reduction (ratio) mid/LCL/UCL and relative reduction (%) mid/LCL/UCL
#   - Theil–Sen slope (pp per year)
#   - Mann–Kendall Z, p-value, and significance flag (p < 0.05)
#
# ------------------------------
# 4) IMPORTANT NOTES / PITFALLS
# ------------------------------
# - If a pollutant column (p) exists but CI columns (p_LCL/p_UCL) are missing, CI fields will be NA.
# - If AF_2013 is 0 (or missing), relative reductions are set to NA to avoid division by zero.
# - Trend tests require at least 3 non-missing annual points; otherwise slope/test outputs are NA.
# - This script assumes "AF" is already expressed in percent (%) in the input tables.
###########################################################################################################
rm(list = ls()); options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(trend)  # sens.slope, mk.test
  suppressWarnings(suppressMessages(requireNamespace("zyp", quietly = TRUE)))
  suppressWarnings(suppressMessages(requireNamespace("Kendall", quietly = TRUE)))
})

# ---- Paths (HIDDEN) ----
AIPW_DIR    <- "PATH/TO/AIPW_RESULTS_DIR"
CLOGIT_DIR  <- "PATH/TO/CLOGIT_RESULTS_DIR"
AIPW_FILE   <- file.path(AIPW_DIR,   "annual_AF_all_pollutants.csv")
CLOGIT_FILE <- file.path(CLOGIT_DIR, "annual_AF_all_pollutants.csv")
OUT_CSV     <- file.path(AIPW_DIR, "annual_AF_decline_summary_2013_2019_endpointCI.csv")

# ---- Switches: include additional pollutants ----
INCLUDE_SO2 <- TRUE
INCLUDE_CO  <- TRUE

# ---- Candidate pollutants ----
pollutants_base <- c("PM25","PM2510","O3","NO2")
pollutants_all  <- pollutants_base
if (INCLUDE_SO2) pollutants_all <- c(pollutants_all, "SO2")
if (INCLUDE_CO)  pollutants_all <- c(pollutants_all, "CO")

# ---- Wide-to-long: extract mid / LCL / UCL (if LCL/UCL missing, fill NA) ----
to_long <- function(file, model_name) {
  if (!file.exists(file)) return(NULL)
  dt <- fread(file)
  if (!"year" %in% names(dt)) stop("Missing column 'year' in: ", file)
  
  pols_in <- intersect(pollutants_all, names(dt))
  if (length(pols_in) == 0L) return(NULL)
  
  lst <- lapply(pols_in, function(p) {
    lcl_col <- paste0(p, "_LCL")
    ucl_col <- paste0(p, "_UCL")
    dt[, .(
      year = year,
      pollutant = p,
      mid = as.numeric(get(p)),
      lcl = if (lcl_col %in% names(dt)) as.numeric(get(lcl_col)) else NA_real_,
      ucl = if (ucl_col %in% names(dt)) as.numeric(get(ucl_col)) else NA_real_,
      model = model_name
    )]
  })
  rbindlist(lst, use.names = TRUE, fill = TRUE)
}

# ---- Read and combine AIPW + clogit AF tables ----
dt_list <- list(
  to_long(AIPW_FILE,   "AIPW"),
  to_long(CLOGIT_FILE, "clogit")
)
af_dt <- rbindlist(Filter(Negate(is.null), dt_list), use.names = TRUE, fill = TRUE)
if (nrow(af_dt) == 0L) stop("No annual AF data loaded. Check file paths and column names.")

# ---- Clean: keep years 2013–2019 (drop 'Total' row if present) ----
af_dt <- af_dt[year != "Total"]
af_dt[, year := suppressWarnings(as.integer(as.character(year)))]
af_dt <- af_dt[!is.na(year) & !is.na(mid)]
af_dt <- af_dt[year >= 2013 & year <= 2019]

# ---- Core computation: endpoint-based decline + robust trend tests ----
calc_one <- function(sub) {
  setorder(sub, year)
  yrs  <- as.numeric(sub$year)
  vals <- as.numeric(sub$mid)
  
  # 2013 and 2019 endpoints (mid/lcl/ucl)
  a13 <- sub[year == 2013][1]
  a19 <- sub[year == 2019][1]
  
  A13m <- if (nrow(a13)) a13$mid else NA_real_
  A13L <- if (nrow(a13)) a13$lcl else NA_real_
  A13U <- if (nrow(a13)) a13$ucl else NA_real_
  
  A19m <- if (nrow(a19)) a19$mid else NA_real_
  A19L <- if (nrow(a19)) a19$lcl else NA_real_
  A19U <- if (nrow(a19)) a19$ucl else NA_real_
  
  # Absolute decline (percentage points), computed by endpoint substitution
  abs_mid <- if (is.na(A13m) | is.na(A19m)) NA_real_ else (A13m - A19m)
  abs_LCL <- if (is.na(A13L) | is.na(A19L)) NA_real_ else (A13L - A19L)
  abs_UCL <- if (is.na(A13U) | is.na(A19U)) NA_real_ else (A13U - A19U)
  
  # Relative decline (ratio), computed by endpoint substitution
  rel_mid <- if (is.na(A13m) | is.na(A19m) | A13m == 0) NA_real_ else (A13m - A19m) / A13m
  rel_LCL <- if (is.na(A13L) | is.na(A19L) | A13L == 0) NA_real_ else (A13L - A19L) / A13L
  rel_UCL <- if (is.na(A13U) | is.na(A19U) | A13U == 0) NA_real_ else (A13U - A19U) / A13U
  
  rel_pct_mid <- if (is.na(rel_mid)) NA_real_ else rel_mid * 100
  rel_pct_LCL <- if (is.na(rel_LCL)) NA_real_ else rel_LCL * 100
  rel_pct_UCL <- if (is.na(rel_UCL)) NA_real_ else rel_UCL * 100
  
  # Robust trend: Theil–Sen slope + Mann–Kendall test
  ts_slope <- NA_real_
  mk_Z <- NA_real_
  mk_p <- NA_real_
  
  idx  <- which(!is.na(yrs) & !is.na(vals))
  yrs2 <- yrs[idx]; vals2 <- vals[idx]
  
  if (length(vals2) >= 3) {
    if (var(vals2) == 0) {
      ts_slope <- 0
      mk_Z <- 0
      mk_p <- 1
    } else {
      # Theil–Sen slope (primary: trend::sens.slope; fallback: zyp::zyp.sen)
      ts <- try(trend::sens.slope(vals2 ~ yrs2), silent = TRUE)
      if (!inherits(ts, "try-error")) {
        if (!is.null(ts$estimates) && ("slope" %in% names(ts$estimates))) {
          ts_slope <- as.numeric(ts$estimates[["slope"]])
        } else if (!is.null(ts$estimates)) {
          ts_slope <- as.numeric(ts$estimates)[1]
        } else if (!is.null(ts$coefficients) && ("yrs2" %in% names(ts$coefficients))) {
          ts_slope <- as.numeric(ts$coefficients[["yrs2"]])
        }
      }
      if (is.na(ts_slope) && requireNamespace("zyp", quietly = TRUE)) {
        ts2 <- try(zyp::zyp.sen(vals2 ~ yrs2), silent = TRUE)
        if (!inherits(ts2, "try-error") && !is.null(ts2$coefficients) && ("yrs2" %in% names(ts2$coefficients))) {
          ts_slope <- as.numeric(ts2$coefficients[["yrs2"]])
        }
      }
      
      # Mann–Kendall test (primary: trend::mk.test; fallback: Kendall::MannKendall)
      mk <- try(trend::mk.test(vals2), silent = TRUE)
      if (!inherits(mk, "try-error")) {
        if (!is.null(mk$estimates) && ("Z" %in% names(mk$estimates))) {
          mk_Z <- as.numeric(mk$estimates[["Z"]])
        } else if (!is.null(mk$statistic)) {
          mk_Z <- as.numeric(mk$statistic)
        }
        mk_p <- suppressWarnings(as.numeric(mk$p.value))
      } else if (requireNamespace("Kendall", quietly = TRUE)) {
        mk2 <- try(Kendall::MannKendall(vals2), silent = TRUE)
        if (!inherits(mk2, "try-error")) {
          tau  <- as.numeric(mk2$tau[1])
          mk_p <- as.numeric(mk2$sl[1])
          mk_Z <- if (is.na(mk_p)) NA_real_ else (sign(tau) * qnorm(1 - mk_p/2))
        }
      }
    }
  }
  
  data.table(
    model      = sub$model[1],
    pollutant  = sub$pollutant[1],
    
    AF_2013_mid = A13m, AF_2013_LCL = A13L, AF_2013_UCL = A13U,
    AF_2019_mid = A19m, AF_2019_LCL = A19L, AF_2019_UCL = A19U,
    
    abs_reduction_pp_mid = abs_mid,
    abs_reduction_pp_LCL = abs_LCL,
    abs_reduction_pp_UCL = abs_UCL,
    
    rel_reduction_mid    = rel_mid,
    rel_reduction_LCL    = rel_LCL,
    rel_reduction_UCL    = rel_UCL,
    
    rel_reduction_pct_mid = rel_pct_mid,
    rel_reduction_pct_LCL = rel_pct_LCL,
    rel_reduction_pct_UCL = rel_pct_UCL,
    
    theilsen_slope_pp_per_year = ts_slope,
    MK_Z = mk_Z,
    MK_p = mk_p,
    sig_p_lt_0_05 = !is.na(mk_p) & (mk_p < 0.05)
  )
}

# ---- Apply per (model, pollutant) ----
res <- af_dt[, calc_one(.SD), by = .(model, pollutant)]

# ---- Reorder columns and sort ----
setcolorder(res, c(
  "model","pollutant",
  "AF_2013_mid","AF_2013_LCL","AF_2013_UCL",
  "AF_2019_mid","AF_2019_LCL","AF_2019_UCL",
  "abs_reduction_pp_mid","abs_reduction_pp_LCL","abs_reduction_pp_UCL",
  "rel_reduction_mid","rel_reduction_LCL","rel_reduction_UCL",
  "rel_reduction_pct_mid","rel_reduction_pct_LCL","rel_reduction_pct_UCL",
  "theilsen_slope_pp_per_year","MK_Z","MK_p","sig_p_lt_0_05"
))
setorder(res, model, pollutant)

# ---- Export ----
fwrite(res, OUT_CSV)
message("Saved: ", OUT_CSV)

# ---- Console preview ----
print(res)
