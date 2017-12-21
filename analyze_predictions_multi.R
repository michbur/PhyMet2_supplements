library(seqinr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(biogram)

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

load("./results/ngram_benchmark_full.RData")
benchmark_raw <- benchmark_res %>% 
  unlist(recursive = FALSE) 

load("./results/ngram_benchmark_full_log.RData")
benchmark_raw_log <- benchmark_res %>% 
  unlist(recursive = FALSE) 

res_sum <- rbind(
  benchmark_raw_log[sapply(benchmark_raw_log, function(i) class(i)[1]) != "try-error"] %>% 
  do.call(rbind, .) %>% 
  mutate(mean_error = sqrt(exp(mean_error)),
         sd_error = sqrt(exp(sd_error)),
         log_val = TRUE),
  benchmark_raw[sapply(benchmark_raw, function(i) class(i)[1]) != "try-error"] %>% 
    do.call(rbind, .) %>% 
    mutate(mean_error = sqrt(mean_error),
           sd_error = sqrt(sd_error),
           log_val = FALSE)) %>% 
  inner_join(filter(conditions_dat, Name %in% all_three) %>% 
               melt(variable.name = "task.id") %>% 
               group_by(task.id) %>% 
               summarise(median_value = median(value))) %>% 
  mutate(frac_mean = mean_error/median_value,
         frac_sd = sd_error/median_value,
         dat_list = as.character(dat_list)) %>% 
  inner_join(training_data) 


best_pars <- group_by(res_sum, task.id) %>% 
  filter(mean_error == min(mean_error)) %>% 
  select(-dat_list) 

best_pars_mcra <- group_by(res_sum, task.id) %>% 
  filter(type1 == "mcra", type2 == "mcra") %>% 
  filter(mean_error == min(mean_error))

best_pars_rna <- group_by(res_sum, task.id) %>% 
  filter(type1 == "rna", type2 == "rna") %>% 
  filter(mean_error == min(mean_error))

save(best_pars, best_pars_mcra, best_pars_rna, file = "./results/best_pars.RData")
