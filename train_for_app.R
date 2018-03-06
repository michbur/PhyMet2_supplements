library(seqinr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(mlr)

source("./functions/validate_seqs.R")

raw_dat <- read.csv("./data/dump_17-07-23.csv") 

mcra_seqs <- list(read.fasta("./raw_data/McrA_1.txt")) %>% 
  unlist(recursive = FALSE) %>% 
  lapply(function(i) {
    res <- i[i != "-" & i != " "]
    mostattributes(res) <- attributes(i)
    res
  }) %>% 
  validate_seqs

rna_seqs <- list(read.fasta("./raw_data/RNA_1.txt")) %>% 
  unlist(recursive = FALSE) %>% 
  lapply(function(i) {
    res <- i[i != "-" & i != " "]
    mostattributes(res) <- attributes(i)
    res
  }) %>% 
  validate_rna



conditions_dat <- raw_dat[c("Name", 
                            "Growth.doubling.time..h.", "Growth.rate", 
                            "Min..growth.temp.", "Max..growth.temp.", 
                            "Min..optimal.growth.temp.", "Max..optimal.growth.temp.",  
                            "Min..growth.NaCl", "Max..growth.NaCl", 
                            "Min..optimal.growth.NaCl", "Max..optimal.growth.NaCl", 
                            "Min..growth.pH", "Max..growth.pH", 
                            "Min..optimal.growth.pH", "Max..optimal.growth.pH")] %>% 
  rename(growth_doubl = Growth.doubling.time..h.,
         growth_rate = Growth.rate,
         min_gt = Min..growth.temp.,
         max_gt = Max..growth.temp.,
         min_ogt = Min..optimal.growth.temp.,
         max_ogt = Max..optimal.growth.temp.,
         min_gn = Min..growth.NaCl,
         max_gn = Max..growth.NaCl,
         min_ogn = Min..optimal.growth.NaCl,
         max_ogn = Max..optimal.growth.NaCl,
         min_gp = Min..growth.pH,
         max_gp = Max..growth.pH,
         min_ogp = Min..optimal.growth.pH,
         max_ogp = Max..optimal.growth.pH) %>% 
  mutate(mean_gt = (min_gt + max_gt)/2,
         mean_ogt = (min_ogt + max_ogt)/2,
         mean_gn = (min_gn + max_gn)/2,
         mean_ogn = (min_ogn + max_ogn)/2,
         mean_gp = (min_gp + max_gp)/2,
         mean_ogp = (min_ogp + max_ogp)/2) %>% 
  select(Name, growth_doubl, growth_rate, mean_ogt, mean_ogn, mean_ogp) %>% 
  na.omit %>% 
  filter(Name != "Methanoculleus sediminis") 

both_mcra_rna <- intersect(unique(rownames(rna_seqs)), unique(rownames(mcra_seqs)))
both_mcra_conditions <- intersect(as.character(conditions_dat[["Name"]]), unique(rownames(mcra_seqs)))
both_rna_conditions <- intersect(as.character(conditions_dat[["Name"]]), unique(rownames(rna_seqs)))
all_three <- intersect(as.character(conditions_dat[["Name"]]), both_mcra_rna)

training_data <- expand.grid(type1 = c("mcra", "rna"), len1 = 1L:6, 
                             type2 = c("mcra", "rna"), len2 = 1L:6) %>% 
  split(1L:nrow(.)) %>% 
  do.call(rbind, .) %>% 
  mutate(dat_list = as.character(rownames(.)))

load("./results/best_pars.RData")
load("./data/seq_dat.RData")

set.seed(15390)
pred_list <- list(rna = lapply(best_pars_rna[["task.id"]], function(ith_condition) {
  train_pars <- best_pars_rna[best_pars_rna[["task.id"]] == ith_condition, ]
  
  ngram_dat <- cbind(select(seq_dat[[train_pars[["type1"]]]][[train_pars[["len1"]]]], -species),
                     seq_dat[[train_pars[["type2"]]]][[train_pars[["len2"]]]])
  
  names_dat <- conditions_dat[, c("Name", ith_condition)] %>% 
    inner_join(ngram_dat, by = c("Name" = "species")) 
  
  dat <- select(names_dat, -Name)
  
  if(train_pars[["log_val"]]) 
    dat[[ith_condition]] <- log(dat[[ith_condition]])
  
  predict_ngrams <- makeRegrTask(id = ith_condition, 
                                 data = dat, 
                                 target = ith_condition)
  
  rf_pars <- list(num.trees = as.numeric(as.character(train_pars[["num.trees"]])),
                  min.node.size = as.numeric(as.character(train_pars[["min.node.size"]])))
  
  learnerRF <- makeFilterWrapper(learner = makeLearner("regr.ranger", par.vals = rf_pars), 
                                 fw.method = "linear.correlation", 
                                 fw.perc = as.numeric(as.character(train_pars[["fw.perc"]])))
  
  train(learnerRF, predict_ngrams)
}))

save(pred_list, file = "./app/pred_list.RData")
