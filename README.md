# Case-crossover-study-air-pollution-mortality
Example Codes for “Reduced mortality burden associated with the clean air policy in China: An integrated assessment based on causal machine-learning”
Overview

This repository provides example codes and simulated (fake) data for reproducing the core analytical workflow of our study, “Reduced mortality burden associated with the clean air policy in China: An integrated assessment based on causal machine-learning.”

Note: The individual-level mortality data used in the real analysis are subject to strict governmental regulations and cannot be publicly released due to privacy and confidentiality requirements. Therefore, we provide a synthetic dataset that mimics the structure of the real time-stratified case-crossover data, as well as example output tables. As a result, all numerical results produced with the example data are for demonstration only and are not expected to match the results reported in the paper.

Software details

R, version 4.0.3 (recommended; newer versions should also work in most cases)


Repository structure (suggested)

Instructions to use data

simulated_data.rds
Simulated time-stratified case-crossover (TS-CCO) dataset with lagged exposures and covariates.


Instructions to use code

Necessary packages are provided in the codes, please first install the packages.
The shared codes are used to generate the main analysis of the paper, codes for drawing pictures are not provided.

1.AIPW.R
Implements the doubly robust AIPW estimator under the time-stratified case-crossover design to estimate pollutant-specific acute mortality effects under correlated multipollutant exposure, and exports the effect estimates (with uncertainty) for downstream attributable-fraction and burden calculations.

2.clogit.R
Fits time-stratified case-crossover conditional logistic regression models as conventional comparators, including single-pollutant, two-pollutant, and multi-pollutant specifications, and outputs harmonized effect estimates (per unit / per IQR, with CIs and FDR-adjusted p values).


3.Daily attributable fractions.R
Uses daily death counts, daily pollutant concentrations, and model-based effect estimates (AME/logRR) with specified reference levels (TMREL) to calculate daily attributable fraction (AF) and attributable number (AN) for each pollutant (including endpoint-based CI propagation), and exports daily AF/AN time series.

4.Annual attributable fractions.R
Aggregates the daily AN (and its CI) to the annual level, then computes annual AF% (and CI) by dividing annual AN by annual total deaths, producing a wide table of annual AF for all pollutants including a “Total” summary row.

5.Mann-Kendall trend test.R
Combines annual AF results from AIPW and clogit, quantifies endpoint-based changes from 2013 to 2019 (absolute/relative reductions with endpoint CIs), and performs robust trend assessment using Theil–Sen slopes and Mann–Kendall tests, exporting a summary table for trend reporting.

6.Total avoided number.R
Estimates the total avoidable (or additional) deaths attributable to temporal changes in pollutant concentrations over 2014–2019, by linking marginal effects (AME/logRR) with concentration trends (e.g., Theil–Sen slopes) and annual deaths, and outputs cause-by-pollutant totals (with CIs) for the integrated health impact summary.


Simulated data are meant only to demonstrate the workflow.

Contact

For questions about code structure, variable definitions, or how to adapt the workflow to other outcomes (e.g., respiratory, cardiovascular, cerebrovascular mortality), please contact the corresponding author or open an issue in the repository.
