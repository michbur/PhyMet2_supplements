# Supplementary code and data

Files description:
*data/seq_dat.RData*: sequences used in training (after data filtering).
*data/condition_data.csv*: culture conditions of strains used in training (after data filtering).
*train_predictor_multi.R* and *train_predictor_multi_log.R*: nested cross-validation of predictors (normal and log-transformed data).
*jackknife.R*: jackknife evaluation of chosen predictors.

# Session info

**R version 3.4.3 (2017-11-30)**

**Platform:** x86_64-pc-linux-gnu (64-bit) 

**locale:**
_LC_CTYPE=pl_PL.UTF-8_, _LC_NUMERIC=C_, _LC_TIME=pl_PL.UTF-8_, _LC_COLLATE=pl_PL.UTF-8_, _LC_MONETARY=pl_PL.UTF-8_, _LC_MESSAGES=pl_PL.UTF-8_, _LC_PAPER=pl_PL.UTF-8_, _LC_NAME=C_, _LC_ADDRESS=C_, _LC_TELEPHONE=C_, _LC_MEASUREMENT=pl_PL.UTF-8_ and _LC_IDENTIFICATION=C_

**attached base packages:** 
_stats_, _graphics_, _grDevices_, _utils_, _datasets_, _methods_ and _base_

**other attached packages:** 
pander(v.0.6.1)

**loaded via a namespace (and not attached):** 
_compiler(v.3.4.3)_, _tools(v.3.4.3)_, _yaml(v.2.1.15)_, _Rcpp(v.0.12.14)_ and _digest(v.0.6.12)_
