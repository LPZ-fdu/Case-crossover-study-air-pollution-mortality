##########################################################################################################
#Codes for "Reduced mortality burden associated with the clean air policy in China: An integrated assessment based on causal machine learning"
#Authors for codes: Yuchang Zhou, Peizheng Li, Ya Gao, Yi Guo, Yixiang Zhu, Lu Zhou, Peng Yin, Haidong Kan, Maigeng Zhou, Renjie Chen
#Correspondence to Renjie Chen, Maigeng Zhou.
###########################################################################################################

##########################################################################################################
# Notes on required inputs and expected outputs
# This script calculates DAILY attributable fraction (AF) and attributable number (AN) for multiple
# pollutants based on:
#   (i) daily baseline death counts for a selected ICD outcome group, and
#   (ii) an exposure–response effect size table (AME and its CI), and
#   (iii) daily pollutant concentrations with a TMREL reference level.
#
# ------------------------------
# 1) REQUIRED INPUT DATA
# ------------------------------
# (A) Case/death record data (CASE_PATH; .rds)
#   - One row per death (or per event) with at least:
#       * date      : date of death/event (character/Date)
#       * BSICD10   : ICD-10 code (character), e.g., "I63", "J44", etc.
#   - The script converts date -> Date, and selects deaths by ICD code ranges
#     (e.g., I60–I69 for cerebrovascular diseases).
#   - Output from this data: DAILY death counts N_daily (Date, N).
#
# (B) Daily pollutant data (POLL_PATH; .csv/.txt/.rds)
#   - One row per day with:
#       * date : calendar date (character/Date), unique and continuous over the study period
#       * pollutant columns matching POLLUTANTS, e.g.:
#           PM25, O3, NO2, PM2510, (optionally SO2, CO)
#   - If PM2510 is absent but PM10 and PM25 exist, it will be constructed as:
#       PM2510 = max(PM10 - PM25, 0)
#
# (C) Effect size table (RR_XLSX; Excel)
#   - One sheet specified by sheet_sel (e.g., "AIPW02") containing:
#       * pollutant  : pollutant name (must match POLLUTANTS)
#       * RR        : average marginal effect for the unit increase 
#       * RR_CI_L   : lower CI for RR
#       * RR_CI_U   : upper CI for RR

#
# ------------------------------
# 2) KEY PARAMETERS
# ------------------------------
# - POLLUTANTS: pollutant list to be evaluated (base + optional SO2/CO switches).
# - UNIT_SCALE: unit conversion factor applied to each pollutant concentration before TMREL:
#       conc_used = poll[[p]] * UNIT_SCALE[p]
# - TMREL: theoretical minimum risk exposure level (reference value) for each pollutant,
#       conc_excess = max(conc_used - TMREL[p], 0)
# - Outcome definition: controlled by ICD selection code blocks (e.g., I60–I69, J00–J99, etc.).
#
# ------------------------------
# 3) EXPECTED OUTPUTS
# ------------------------------
# The script writes three CSV files to OUT_DIR:
#
# (1) attributable_fraction_daily_all_pollutants.csv
#     - Columns:
#         Date,
#         <pollutant>, <pollutant>_LCL, <pollutant>_UCL  for each pollutant in POLLUTANTS
#     - AF is expressed in percentage (%).
#
# (2) attributable_number_daily_all_pollutants.csv
#     - Columns:
#         Date,
#         <pollutant>, <pollutant>_LCL, <pollutant>_UCL  for each pollutant
#     - AN is the number of deaths attributable to pollutant exposure above TMREL on that day.
#
# (3) daily_death_number.csv
#     - Columns:
#         Date, N
#     - N is the observed death count for the selected ICD outcome group on that day.
#
# ------------------------------
###########################################################################################################
rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(readxl)
})

# ---- Configuration (HIDDEN PATHS) ----
CASE_PATH <- "PATH/TO/CASE_DATA.rds"
POLL_PATH <- "PATH/TO/DAILY_POLLUTANTS_2013_2019.csv"
AME_XLSX  <- "PATH/TO/RESULTS.xlsx"     # sheet: pollutant, AME, AME_CI_L, AME_CI_U
OUT_DIR   <- "PATH/TO/OUTPUT_DIR"
sheet_sel <- "AIPW02"

# Optional fixed date range (not used in this script; kept for compatibility)
DATE_START <- as.Date("2013-01-01")
DATE_END   <- as.Date("2019-12-31")

# ---- Switches: include additional pollutants ----
INCLUDE_SO2 <- TRUE
INCLUDE_CO  <- TRUE

# ---- Pollutant list (PM10 -> PM2510) ----
POLLUTANTS_BASE <- c("PM25","O3","NO2","PM2510")
POLLUTANTS <- POLLUTANTS_BASE
if (INCLUDE_SO2) POLLUTANTS <- c(POLLUTANTS, "SO2")
if (INCLUDE_CO)  POLLUTANTS <- c(POLLUTANTS, "CO")

# Unit scaling (adjust if the raw units differ)
UNIT_SCALE <- c(PM25=1, O3=1, NO2=1, PM2510=1, SO2=1, CO=1)

# TMREL thresholds (must be in the same unit as concentrations after scaling)
TMREL <- c(PM25=0.0, O3=0.0, NO2=0.0, PM2510=0.0, SO2=0.0, CO=0.0)

# ---- Read case data ----
cases_all <- readRDS(CASE_PATH); setDT(cases_all)
cases_all[, Date := as.Date(date)]
cases_all[, icd_chr := toupper(trimws(BSICD10))]

tmp <- toupper(cases_all$BSICD10)
num <- suppressWarnings(as.numeric(sub("^I", "", tmp)))

# Outcome selection: I60–I69 (cerebrovascular diseases)
sel_cases <- cases_all[!is.na(num) & grepl("^I", tmp) & num >= 60.0 & num <= 69.9]

# Alternative example: J00–J99 (respiratory diseases)
# sel_cases <- cases_all[!is.na(num) & grepl("^J", tmp) & num >= 0.0 & num <= 99.9]

# I00–I52
# sel_cases <- cases_all[
#   !is.na(num) & grepl("^I", tmp) &
#     (
#       (num >= 20.0 & num <= 25.9) |
#       (num >= 10.0 & num <= 15.9) |
#       (num >= 26.0 & num <= 28.9) |
#       (num >= 30.0 & num <= 52.9) |
#       (num >= 5.0 & num <= 9.9) |
#       (num >= 0.0 & num <= 2.9)
#     ),
# ]

rm(cases_all)

# ---- Read pollutant data (CSV/RDS auto-detect) ----
read_any <- function(path){
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("rds","rda","rdata")) return(readRDS(path))
  if (ext %in% c("csv","txt")) return(fread(path))
  stop("Unsupported file type: ", path)
}

poll <- read_any(POLL_PATH); setDT(poll)
stopifnot("date" %in% names(poll))
poll[, Date := as.Date(date)]
setorder(poll, Date)
stopifnot(all(!duplicated(poll$Date)))

# Ensure PM2510 exists: if missing but PM10 and PM25 exist, construct it on the fly
if (!("PM2510" %in% names(poll))) {
  if (all(c("PM10","PM25") %in% names(poll))) {
    poll[, PM2510 := pmax(PM10 - PM25, 0)]
  } else {
    stop("Cannot create PM2510: missing PM2510 and also missing PM10 and/or PM25.")
  }
}

# Validate required pollutant columns
if (!all(POLLUTANTS %in% names(poll))) {
  stop("Missing pollutant columns: ", paste(setdiff(POLLUTANTS, names(poll)), collapse = ", "))
}

# ---- Read AME table (AME equals RR - 1) ----
ame_raw <- as.data.table(read_excel(AME_XLSX, sheet = sheet_sel))
need_cols <- c("pollutant","AME","AME_CI_L","AME_CI_U")
if (!all(need_cols %in% names(ame_raw))) {
  stop("Missing columns in AME sheet: ", paste(setdiff(need_cols, names(ame_raw)), collapse = ", "))
}

# Keep only the selected pollutants (and de-duplicate by pollutant)
ame_tab <- ame_raw[pollutant %in% POLLUTANTS, .SD[1], by = pollutant]

# ---- Daily death counts for the selected ICD outcome ----
all_dates <- data.table(Date = poll$Date)
N_daily   <- sel_cases[, .N, by = .(Date)]
N_daily   <- merge(all_dates, N_daily, by = "Date", all.x = TRUE)
N_daily[is.na(N), N := 0L]
setorder(N_daily, Date)

# ---- Compute daily AN and AF (with endpoint CI propagation via AME CI) ----
AN_dt <- data.table(Date = all_dates$Date)
AF_dt <- data.table(Date = all_dates$Date)

for (p in POLLUTANTS) {
  r <- ame_tab[pollutant == p]
  if (nrow(r) != 1L) {
    warning("AME row is missing or not unique for: ", p, " (skipping).")
    next
  }
  
  # Convert AME to RR and log(RR)
  RR   <- 1 + as.numeric(r$AME)
  RR_L <- 1 + as.numeric(r$AME_CI_L)
  RR_U <- 1 + as.numeric(r$AME_CI_U)
  
  lnRR <- log(RR)
  lnL  <- log(RR_L)
  lnU  <- log(RR_U)
  
  # Concentrations: apply unit scaling and TMREL thresholding
  conc <- as.numeric(poll[[p]]) * UNIT_SCALE[[p]]
  conc <- pmax(conc - TMREL[[p]], 0)
  
  # Attributable fraction on day i: AF = (exp(lnRR * x) - 1) / exp(lnRR * x)
  p_mid <- (exp(lnRR * conc) - 1) / exp(lnRR * conc)
  p_lo  <- (exp(lnL  * conc) - 1) / exp(lnL  * conc)
  p_hi  <- (exp(lnU  * conc) - 1) / exp(lnU  * conc)
  
  p_mid[!is.finite(p_mid)] <- 0
  p_lo [!is.finite(p_lo )] <- 0
  p_hi [!is.finite(p_hi )] <- 0
  
  # Attributable number (AN) and AF (%), using the daily death counts as baseline
  AN_mid <- p_mid * N_daily$N
  AN_lo  <- p_lo  * N_daily$N
  AN_hi  <- p_hi  * N_daily$N
  
  AF_mid <- ifelse(N_daily$N > 0, AN_mid / N_daily$N * 100, 0)
  AF_lo  <- ifelse(N_daily$N > 0, AN_lo  / N_daily$N * 100, 0)
  AF_hi  <- ifelse(N_daily$N > 0, AN_hi  / N_daily$N * 100, 0)
  
  # Write pollutant-specific columns (mid / LCL / UCL)
  AN_dt[, (p)              := AN_mid]
  AN_dt[, paste0(p,"_LCL") := AN_lo]
  AN_dt[, paste0(p,"_UCL") := AN_hi]
  
  AF_dt[, (p)              := AF_mid]
  AF_dt[, paste0(p,"_LCL") := AF_lo]
  AF_dt[, paste0(p,"_UCL") := AF_hi]
}

# ---- Export ----
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
fwrite(AF_dt,  file.path(OUT_DIR, "attributable_fraction_daily_all_pollutants.csv"))
fwrite(AN_dt,  file.path(OUT_DIR, "attributable_number_daily_all_pollutants.csv"))
fwrite(N_daily, file.path(OUT_DIR, "daily_death_number.csv"))

cat("Done. Daily AF and AN (with CI) were saved to: ", OUT_DIR, "\n")
