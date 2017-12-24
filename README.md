# Supplementary code and data

The code below was used to train [http://www.smorfland.uni.wroc.pl/shiny/mgp/](MethanoGram): predictor of methanograms' culture conditions.

Files description:
*data/seq_dat.RData*: sequences used in the training and evaluation of MethanoGram.
*data/condition_data.csv*: culture conditions of strains used in the training and evaluation of MethanoGram.
*train_predictor_multi.R* and *train_predictor_multi_log.R*: nested cross-validation of MethanoGram (normal and log-transformed data). Results are respectively: *ngram_benchmark_full.RData* and *ngram_benchmark_full_log.RData*.
*results/best_pars.RData* - hyperparameters of the final MethanoGram iterations.
*jackknife.R*: jackknife evaluation of the final MethanoGram iterations.

# Session info

**R version 3.4.3 (2017-11-30)**

**Platform:** x86_64-pc-linux-gnu (64-bit) 

**locale:**
_LC_CTYPE=pl_PL.UTF-8_, _LC_NUMERIC=C_, _LC_TIME=pl_PL.UTF-8_, _LC_COLLATE=pl_PL.UTF-8_, _LC_MONETARY=pl_PL.UTF-8_, _LC_MESSAGES=pl_PL.UTF-8_, _LC_PAPER=pl_PL.UTF-8_, _LC_NAME=C_, _LC_ADDRESS=C_, _LC_TELEPHONE=C_, _LC_MEASUREMENT=pl_PL.UTF-8_ and _LC_IDENTIFICATION=C_

**attached base packages:** 
_grid_, _stats_, _graphics_, _grDevices_, _utils_, _datasets_, _methods_ and _base_

**other attached packages:** 
_gridExtra(v.2.3)_, _reshape2(v.1.4.3)_, _ggplot2(v.2.2.1)_, _mlr(v.2.11)_, _ParamHelpers(v.1.10)_, _pbapply(v.1.3-3)_, _seqinr(v.3.4-5)_, _biogram(v.1.5)_, _slam(v.0.1-40)_ and _dplyr(v.0.7.4)_

**loaded via a namespace (and not attached):** 
_parallelMap(v.1.3)_, _gmp(v.0.5-13.1)_, _Rcpp(v.0.12.14)_, _plyr(v.1.8.4)_, _compiler(v.3.4.3)_, _bindr(v.0.1)_, _tools(v.3.4.3)_, _partitions(v.1.9-19)_, _digest(v.0.6.12)_, _bit(v.1.1-12)_, _lattice(v.0.20-35)_, _tibble(v.1.3.4)_, _checkmate(v.1.8.5)_, _gtable(v.0.2.0)_, _pkgconfig(v.2.0.1)_, _rlang(v.0.1.4)_, _Matrix(v.1.2-11)_, _yaml(v.2.1.15)_, _parallel(v.3.4.3)_, _polynom(v.1.3-9)_, _bindrcpp(v.0.2)_, _stringr(v.1.2.0)_, _combinat(v.0.0-8)_, _ade4(v.1.7-8)_, _glue(v.1.2.0)_, _data.table(v.1.10.4-3)_, _R6(v.2.2.2)_, _survival(v.2.41-3)_, _pander(v.0.6.1)_, _magrittr(v.1.5)_, _splines(v.3.4.3)_, _BBmisc(v.1.11)_, _backports(v.1.1.1)_, _scales(v.0.5.0.9000)_, _assertthat(v.0.2.0)_, _colorspace(v.1.3-2)_, _stringi(v.1.1.6)_, _entropy(v.1.2.1)_, _lazyeval(v.0.2.1)_ and _munsell(v.0.4.3)_
