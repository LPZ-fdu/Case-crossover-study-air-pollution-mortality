##########################################################################################################
#Codes for "Reduced mortality burden associated with the clean air policy in China: An integrated assessment based on causal machine learning"
#Authors for codes: Yuchang Zhou, Peizheng Li, Ya Gao, Yi Guo, Yixiang Zhu, Lu Zhou, Peng Yin, Haidong Kan, Maigeng Zhou, Renjie Chen
#Correspondence to Renjie Chen, Maigeng Zhou.
###########################################################################################################

##########################################################################################################
# Notes on required inputs and expected outputs
#
# This script is a post-processing + visualization workflow for quantifying avoidable deaths under a
# slope-based counterfactual scenario and producing Figure 4.
#
# It combines:
#   (1) marginal exposure–response effects (AME: average marginal effect; AME = RR - 1 per 1 unit),
#   (2) observed daily pollutant concentrations (to estimate annual trend using Theil–Sen slope),
#   (3) observed daily cause-specific death counts (to aggregate annual deaths),
# and then estimates avoidable deaths for 2014–2019 relative to the BASE_YEAR (2013) baseline.
#
# ------------------------------
# 1) REQUIRED INPUT DATA
# ------------------------------
# (A) AME_WORKBOOK.xlsx (FILE_AME_XLSX; sheet = SHEET_AME)
#   - Required columns (case-sensitive after trimming):
#       name        : outcome/cause name (must match cause columns in daily deaths table)
#       pollutant   : pollutant name (e.g., PM25, PM2510, O3, NO2, optionally SO2, CO)
#       AME         : average marginal effect, defined as (RR - 1) per 1-unit increase in pollutant
#       AME_CI_L    : lower bound of AME
#       AME_CI_U    : upper bound of AME
#   - Interpretation:
#       beta = log(RR) = log(1 + AME) is used as the log risk ratio per 1 unit.
#   - Constraints:
#       AME, AME_CI_L, AME_CI_U must be > -1 (otherwise log1p is invalid).
#
# (B) Daily pollutant concentrations file (FILE_POLL_DAILY; CSV)
#   - Must contain a date column: "date" (any case is accepted; the script lowercases names)
#   - Must contain daily concentration columns for each pollutant in POLLUTANTS
#     (lowercase names are expected after setnames(tolower)).
#     Example required columns when include_SO2_CO = TRUE:
#       pm25, pm2510, o3, no2, so2, co
#   - Time range:
#       BASE_YEAR..YEAR_MAX (default 2013–2019). Values outside will be filtered out.
#   - Used for:
#       computing annual mean concentrations, then Theil–Sen slope (delta_x) of annual mean vs year.
#
# (C) Daily deaths table (FILE_DEATHS_XLSX; sheet = SHEET_DEATHS)
#   - Must contain one date column (automatically detected among: Date/date/day/time/日期)
#   - Other columns are treated as different causes/outcomes ("name"), each being daily death counts.
#   - The cause column names must match AME$name exactly (after trimming); otherwise those causes are ignored.
#   - Time range:
#       BASE_YEAR..YEAR_MAX (default 2013–2019). Values outside will be filtered out.
#
# ------------------------------
# 2) CORE CALCULATION (SLOPE-BASED COUNTERFACTUAL)
# ------------------------------
# Step 1: Convert AME to log(RR) per 1-unit increase:
#   beta = log(1 + AME)
#
# Step 2: Estimate annual pollutant trend (Theil–Sen slope) from annual mean concentrations:
#   delta_x = SenSlope( annual_mean ~ year )
#   Interpreted as "change in annual mean concentration per year".
#
# Step 3: For each year y = BASE_YEAR+1..YEAR_MAX, define t = y - BASE_YEAR and compute:
#   AD_y = ( exp(beta * delta_x * t) - 1 ) * N_y
# where:
#   N_y is total deaths in year y for a given cause ("name").
#
# Endpoint year BASE_YEAR is treated as the baseline anchor; the counterfactual is driven by a linear
# concentration change implied by delta_x. CI is propagated by substituting beta_L and beta_U
# (derived from AME_CI_L / AME_CI_U).
#
# ------------------------------
# 3) EXPECTED OUTPUTS
# ------------------------------
# (A) Excel workbook (OUT_XLSX) with two sheets:
#   - "avoidable_deaths_yearly":
#       avoidable deaths by (name, pollutant, year) for years BASE_YEAR+1..YEAR_MAX, with CI
#   - "avoidable_deaths_total":
#       total avoidable deaths over 2014–2019 by (name, pollutant), with CI and supporting fields
#
# (B) Figure 4 (OUT_PNG):
#   - Bar chart by cause-system (e.g., Respiratory/Cardiovascular/Cerebrovascular)
#   - Includes a "TOTAL" bar (sum across pollutants within each system)
#   - Visualization convention:
#       non-O3 pollutants plotted above zero (positive)
#       O3 plotted below zero (negative), indicating additional burden under this convention
#
# ------------------------------
# 4) IMPORTANT NOTES / PITFALLS
# ------------------------------
# - This script does NOT compute daily attributable fractions (AF) directly from daily concentrations
#   and a TMREL; instead it estimates avoidable deaths using annual concentration trends (delta_x).
# - Ensure pollutant naming is consistent across AME (pollutant), pollutant daily file (column names),
#   and plotting order. The script uppercases AME pollutant names and lowercases daily pollutant columns.
# - If any pollutant has missing annual means or insufficient years, delta_x may become NA and all
#   derived avoidable deaths for that pollutant will be NA.
# - The plotting section hard-filters outcomes to `desired_names`; adjust as needed.
###########################################################################################################
rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(writexl)
  library(lubridate)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ===================== User switches =====================
# FALSE: PM25 / PM2510 / O3 / NO2 only
# TRUE : also include SO2 and CO
include_SO2_CO <- TRUE

# ===================== Paths (HIDDEN) =====================
FILE_AME_XLSX    <- "PATH/TO/AME_WORKBOOK.xlsx"       # contains sheet AME
SHEET_AME        <- "AME"

FILE_POLL_DAILY  <- "PATH/TO/POLL_DAILY_2013_2019.csv" # daily mean concentrations
FILE_DEATHS_XLSX <- "PATH/TO/DAILY_DEATHS.xlsx"        # daily deaths (wide table)
SHEET_DEATHS     <- 1

OUT_XLSX         <- "PATH/TO/OUTPUT/AD_2013_2019.xlsx"
OUT_PNG          <- "PATH/TO/OUTPUT/Figure4.png"

BASE_YEAR <- 2013L
YEAR_MAX  <- 2019L

POLLUTANTS <- if (include_SO2_CO) {
  c("PM25","PM2510","O3","NO2","SO2","CO")
} else {
  c("PM25","PM2510","O3","NO2")
}

# ===================== Helper: Theil–Sen slope =====================
theilsen_slope <- function(year_vec, value_vec) {
  year_vec  <- as.numeric(year_vec)
  value_vec <- as.numeric(value_vec)
  ok <- is.finite(year_vec) & is.finite(value_vec)
  year_vec <- year_vec[ok]; value_vec <- value_vec[ok]
  n <- length(year_vec)
  if (n < 2) return(NA_real_)
  
  slps <- c()
  k <- 1L
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (year_vec[j] != year_vec[i]) {
        slps[k] <- (value_vec[j] - value_vec[i]) / (year_vec[j] - year_vec[i])
        k <- k + 1L
      }
    }
  }
  stats::median(slps, na.rm = TRUE)
}

# ===================== 1) Read AME and convert to log(RR) =====================
ame <- read_xlsx(FILE_AME_XLSX, sheet = SHEET_AME) %>% as.data.table()
setnames(ame, old = names(ame), new = trimws(names(ame)))

need_cols <- c("name","pollutant","AME","AME_CI_L","AME_CI_U")
stopifnot(all(need_cols %in% names(ame)))

ame[, pollutant := toupper(trimws(pollutant))]
ame <- ame[pollutant %in% toupper(POLLUTANTS)]

# Safety check: log1p requires AME > -1
if (any(is.finite(ame$AME) & ame$AME <= -1, na.rm = TRUE)) {
  stop("Found AME <= -1, which makes log1p(AME) invalid. Please check AME inputs.")
}
if (any(is.finite(ame$AME_CI_L) & ame$AME_CI_L <= -1, na.rm = TRUE)) {
  stop("Found AME_CI_L <= -1, which makes log1p(AME_CI_L) invalid. Please check AME inputs.")
}
if (any(is.finite(ame$AME_CI_U) & ame$AME_CI_U <= -1, na.rm = TRUE)) {
  stop("Found AME_CI_U <= -1, which makes log1p(AME_CI_U) invalid. Please check AME inputs.")
}

# AME (= RR - 1 per 1 unit) -> beta = log(RR) = log(1 + AME)
ame[, beta_m   := log1p(AME)]
ame[, beta_m_L := log1p(AME_CI_L)]
ame[, beta_m_U := log1p(AME_CI_U)]

# ===================== 2) Daily pollutant -> annual mean -> Theil–Sen slope =====================
poll <- fread(FILE_POLL_DAILY)
setnames(poll, tolower(names(poll)))
if (!"date" %in% names(poll)) stop("Pollutant file must contain a 'date' column.")
poll[, date := as.Date(date)]
poll[, year := year(date)]
poll <- poll[year >= BASE_YEAR & year <= YEAR_MAX]

poll_cols_lower <- tolower(POLLUTANTS)
miss_cols <- setdiff(poll_cols_lower, names(poll))
if (length(miss_cols) > 0) stop("Missing pollutant columns: ", paste(miss_cols, collapse = ", "))

poll_year <- melt(
  poll[, c("year", poll_cols_lower), with = FALSE],
  id.vars = "year", variable.name = "pollutant", value.name = "value"
)[, .(annual_mean = mean(value, na.rm = TRUE)), by = .(pollutant, year)]

slope_tab <- poll_year[, .(delta_x = theilsen_slope(year, annual_mean)), by = pollutant]
slope_tab[, pollutant := toupper(pollutant)]

# ===================== 3) Daily deaths (wide) -> annual deaths =====================
deaths_wide <- read_xlsx(FILE_DEATHS_XLSX, sheet = SHEET_DEATHS) %>% as.data.table()
setnames(deaths_wide, old = names(deaths_wide), new = trimws(names(deaths_wide)))

date_col_idx <- which(tolower(names(deaths_wide)) %in% c("date","day","time","日期"))
if (length(date_col_idx) == 0) stop("Daily deaths file must contain a date column (e.g., Date/日期).")
setnames(deaths_wide, names(deaths_wide)[date_col_idx[1]], "Date")
deaths_wide[, Date := as.Date(Date)]

deaths_long <- melt(deaths_wide, id.vars = "Date", variable.name = "name", value.name = "deaths")
deaths_long[, deaths := as.numeric(deaths)]
deaths_long <- deaths_long[is.finite(deaths)]
deaths_long[, year := year(Date)]
deaths_long <- deaths_long[year >= BASE_YEAR & year <= YEAR_MAX]

deaths_annual <- deaths_long[, .(N_deaths = sum(deaths, na.rm = TRUE)), by = .(name, year)]

# ===================== 4) Compute avoidable deaths using slope-based counterfactual =====================
years_vec <- seq(BASE_YEAR + 1L, YEAR_MAX, by = 1L)
t_map <- data.table(year = years_vec, t = years_vec - BASE_YEAR)

names_common <- intersect(unique(ame$name), unique(deaths_annual$name))
if (length(names_common) == 0) {
  stop("No overlap between AME$name and the cause names in the deaths table. Please align naming.")
}
if (length(setdiff(unique(ame$name), names_common)) > 0) {
  message("Warning: these AME$name are not found in deaths table and will be ignored:\n",
          paste(setdiff(unique(ame$name), names_common), collapse = ", "))
}

ame <- ame[name %in% names_common]
deaths_annual <- deaths_annual[name %in% names_common]

grid <- CJ(
  name = sort(unique(ame$name)),
  pollutant = sort(unique(ame$pollutant)),
  year = years_vec
)

grid <- merge(grid, ame[, .(name, pollutant, beta_m, beta_m_L, beta_m_U)],
              by = c("name","pollutant"), all.x = TRUE)
grid <- merge(grid, slope_tab, by = "pollutant", all.x = TRUE)
grid <- merge(grid, deaths_annual, by = c("name","year"), all.x = TRUE)
grid <- merge(grid, t_map, by = "year", all.x = TRUE)

grid[is.na(N_deaths), N_deaths := 0]

# Avoidable deaths (slope-based; endpoint year BASE_YEAR fixed)
# AD_y = (exp(beta * delta_x * t) - 1) * N_y
calc_piece <- function(beta, dx, t, N) {
  if (!is.finite(beta) || !is.finite(dx) || !is.finite(t) || !is.finite(N)) return(NA_real_)
  (exp(beta * dx * t) - 1) * N
}

grid[, avoid_y   := calc_piece(beta_m,   delta_x, t, N_deaths)]
grid[, avoid_y_L := calc_piece(beta_m_L, delta_x, t, N_deaths)]
grid[, avoid_y_U := calc_piece(beta_m_U, delta_x, t, N_deaths)]

res_yearly <- grid[, .(
  beta_marginal     = beta_m,
  beta_marginal_L   = beta_m_L,
  beta_marginal_U   = beta_m_U,
  delta_x_per_year  = delta_x,
  N_deaths_year     = N_deaths,
  t                = t,
  avoidable_deaths_year   = avoid_y,
  avoidable_deaths_year_L = avoid_y_L,
  avoidable_deaths_year_U = avoid_y_U
), by = .(name, pollutant, year)]

res_total <- res_yearly[, .(
  years_included           = paste(range(year, na.rm = TRUE), collapse = "–"),
  delta_x_per_year         = unique(delta_x_per_year)[1],
  beta_marginal            = unique(beta_marginal)[1],
  beta_marginal_L          = unique(beta_marginal_L)[1],
  beta_marginal_U          = unique(beta_marginal_U)[1],
  deaths_total_2014_2019   = sum(N_deaths_year, na.rm = TRUE),
  avoidable_deaths_total   = sum(avoidable_deaths_year,   na.rm = TRUE),
  avoidable_deaths_total_L = sum(avoidable_deaths_year_L, na.rm = TRUE),
  avoidable_deaths_total_U = sum(avoidable_deaths_year_U, na.rm = TRUE)
), by = .(name, pollutant)]

setorder(res_yearly, name, pollutant, year)
setorder(res_total,  name, pollutant)

write_xlsx(list(
  "avoidable_deaths_yearly" = res_yearly,
  "avoidable_deaths_total"  = res_total
), path = OUT_XLSX)
message("Done. Results written to: ", OUT_XLSX)

# ===================== 5) Plot: centered at 0 with TOTAL bar; O3 shown below 0 =====================
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(forcats)
  library(grid)
})

dt <- copy(res_total)

desired_names <- c(
  "Respiratory system",
  "Cardiovascular system",
  "Cerebrovascular system"
)
dt <- dt[name %in% desired_names]

dt[, name_f := factor(name, levels = desired_names)]
dt[, group_id := as.numeric(name_f)]

dt[, pollutant := toupper(pollutant)]
dt[, AD := as.numeric(avoidable_deaths_total)]

# Direction convention for visualization:
# - Non-O3: shown above 0 (avoidable, positive)
# - O3   : shown below 0 (additional burden, negative)
dt[, AD_plot := ifelse(pollutant == "O3", -abs(AD), abs(AD))]

# Add per-system TOTAL (net) bar
total_dt <- dt[, .(
  pollutant              = "TOTAL",
  years_included         = paste(unique(years_included), collapse = "; "),
  delta_x_per_year       = NA_real_,
  beta_marginal          = NA_real_,
  beta_marginal_L        = NA_real_,
  beta_marginal_U        = NA_real_,
  deaths_total_2014_2019 = NA_real_,
  avoidable_deaths_total = NA_real_,
  avoidable_deaths_total_L = NA_real_,
  avoidable_deaths_total_U = NA_real_,
  AD                     = sum(AD_plot, na.rm = TRUE),
  AD_plot                = sum(AD_plot, na.rm = TRUE)
), by = .(name, name_f, group_id)]

dt <- rbind(dt, total_dt, fill = TRUE)

group_gap <- 4.3
dt[, x_group := (group_id - 1) * group_gap + 1]

all_pol <- toupper(unique(dt$pollutant))
pol_order_full <- c("TOTAL","PM25","PM2510","NO2","SO2","CO","O3")
pol_levels <- pol_order_full[pol_order_full %in% all_pol]
dt[, pollutant := factor(toupper(pollutant), levels = pol_levels)]

pm_labels_expr <- c(
  TOTAL  = "\"Total\"",
  PM25   = "PM[2.5]",
  PM2510 = "PM[2.5-10]",
  NO2    = "NO[2]",
  SO2    = "SO[2]",
  CO     = "CO",
  O3     = "O[3]"
)
pm_labels_expr <- pm_labels_expr[pol_levels]

n_pol   <- length(pol_levels)
offset  <- 0.60
offsets <- seq(-(n_pol - 1) / 2, (n_pol - 1) / 2, by = 1) * offset

offset_dt <- data.table(
  pollutant = factor(pol_levels, levels = pol_levels),
  x_off = offsets
)

dt_plot <- merge(dt, offset_dt, by = "pollutant", all.x = TRUE)
dt_plot[, x_plot := x_group + x_off]

max_pos <- max(dt_plot$AD_plot[dt_plot$AD_plot > 0], na.rm = TRUE)
max_neg <- abs(min(dt_plot$AD_plot[dt_plot$AD_plot < 0], na.rm = TRUE))
if (!is.finite(max_pos)) max_pos <- 0
if (!is.finite(max_neg)) max_neg <- 0

up_limit   <- max(20000, ceiling(max_pos / 20000) * 20000)
down_limit <- max(20000, ceiling(max_neg / 20000) * 20000)

y_lim    <- c(-1.05 * down_limit, 1.05 * up_limit)
y_breaks <- seq(-down_limit, up_limit, by = 40000)

fmt_num <- label_number(big.mark = ",", accuracy = 1)

pad_up   <- 0.04 * up_limit
pad_down <- 0.04 * down_limit

dt_plot[, label_val := ifelse(
  pollutant == "O3" & AD_plot < 0,
  paste0("-", fmt_num(abs(AD))),
  fmt_num(abs(AD))
)]

dt_plot[, y_lab := ifelse(
  AD_plot >= 0,
  AD_plot + pad_up,
  AD_plot - 2.0 * pad_down
)]

fill_cols_full <- c(
  TOTAL  = "#7E6148FF",
  PM25   = "#E64B35FF",
  PM2510 = "#4DBBD5FF",
  NO2    = "#3C5488FF",
  SO2    = "#F39B7FFF",
  CO     = "#8491B4FF",
  O3     = "#00A087FF"
)
fill_cols <- fill_cols_full[pol_levels]

p_stack_center <- ggplot(dt_plot, aes(x = x_plot, y = AD_plot, fill = pollutant)) +
  geom_hline(yintercept = 0, color = "black", size = 0.6) +
  geom_col(width = 0.30, color = "white", size = 0.25, na.rm = TRUE) +
  geom_text(
    aes(x = x_plot, y = y_lab, label = label_val),
    size = 3.3,
    family = "Times New Roman",
    color = "black",
    inherit.aes = FALSE
  ) +
  scale_x_continuous(
    breaks = sort(unique(dt_plot$x_group)),
    labels = levels(dt_plot$name_f),
    expand = expansion(mult = c(0.10, 0.10))
  ) +
  scale_y_continuous(
    limits = y_lim,
    breaks = y_breaks,
    labels = fmt_num,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_fill_manual(
    values = fill_cols,
    breaks = pol_levels,
    labels = function(x) parse(text = unname(pm_labels_expr[as.character(x)])),
    name = "Pollutant"
  ) +
  labs(
    x = "Cause-specific mortality",
    y = "Avoidable number of deaths",
    title = NULL
  ) +
  theme_bw(base_size = 14, base_family = "Times New Roman") +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black", angle = 0, hjust = 0.5, vjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "top",
    legend.key.height = unit(0.7, "lines"),
    legend.key.width  = unit(1.3, "lines"),
    legend.box = "horizontal",
    legend.box.just = "center",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave(OUT_PNG, p_stack_center, width = 8, height = 5, dpi = 300)
message("Done. Figure saved to: ", OUT_PNG)
