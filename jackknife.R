library(seqinr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(mlr)

configureMlr(show.info = FALSE)

conditions_dat <- read.csv("./data/condition_data.csv")

load("./data/seq_dat.RData")

training_data <- expand.grid(type1 = c("mcra", "rna"), len1 = 1L:6, 
                             type2 = c("mcra", "rna"), len2 = 1L:6) %>% 
  split(1L:nrow(.)) %>% 
  do.call(rbind, .) %>% 
  mutate(dat_list = as.character(rownames(.)))

load("./results/best_pars.RData")

# train with known random seed to make random forests reproducible -------------------
set.seed(15390)
jackknife_res <- lapply(best_pars[["task.id"]], function(ith_condition) {
  train_pars <- best_pars[best_pars[["task.id"]] == ith_condition, ]
  
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
  
  raw_jk <- resample(learnerRF, predict_ngrams, resampling = makeResampleDesc("CV", iters = nrow(dat)),
                     keep.pred = TRUE)
  
  res <- select(names_dat, Name) %>% 
    mutate(id = 1L:nrow(.)) %>% 
    inner_join(data.frame(getRRPredictions(raw_jk))) %>% 
    select(Name, truth, response)
  
  if(train_pars[["log_val"]]) {
    res[["response"]] <- exp(res[["response"]])
    res[["truth"]] <- exp(res[["truth"]])
  }
  
  mutate(res, condition = ith_condition)
})

jackknife_res_rna <- lapply(best_pars_rna[["task.id"]], function(ith_condition) {
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
  
  raw_jk <- resample(learnerRF, predict_ngrams, resampling = makeResampleDesc("CV", iters = nrow(dat)),
                     keep.pred = TRUE)
  
  res <- select(names_dat, Name) %>% 
    mutate(id = 1L:nrow(.)) %>% 
    inner_join(data.frame(getRRPredictions(raw_jk))) %>% 
    select(Name, truth, response)
  
  if(train_pars[["log_val"]]) {
    res[["response"]] <- exp(res[["response"]])
    res[["truth"]] <- exp(res[["truth"]])
  }
  
  mutate(res, condition = ith_condition)
})

jackknife_res_mcra <- lapply(best_pars_mcra[["task.id"]], function(ith_condition) {
  train_pars <- best_pars_mcra[best_pars_mcra[["task.id"]] == ith_condition, ]
  
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
  
  raw_jk <- resample(learnerRF, predict_ngrams, resampling = makeResampleDesc("CV", iters = nrow(dat)),
                     keep.pred = TRUE)
  
  res <- select(names_dat, Name) %>% 
    mutate(id = 1L:nrow(.)) %>% 
    inner_join(data.frame(getRRPredictions(raw_jk))) %>% 
    select(Name, truth, response)
  
  if(train_pars[["log_val"]]) {
    res[["response"]] <- exp(res[["response"]])
    res[["truth"]] <- exp(res[["truth"]])
  }
  
  mutate(res, condition = ith_condition)
})


do.call(rbind, jackknife_res) %>% 
  write.csv("final_jackknife.csv", row.names = FALSE)

median_pred <- function(x) {
  sapply(1L:length(x), function(i) sqrt((x[i] - median(x[-i]))^2))
}

res <- rbind(read.csv("final_jackknife.csv") %>% 
               mutate(error = sqrt((truth - response)^2)) %>% 
               group_by(condition) %>% 
               summarise(error = mean(error)) %>% 
               mutate(type = "MethanoGram"),
             do.call(rbind, jackknife_res_mcra) %>% 
               mutate(error = sqrt((truth - response)^2)) %>% 
               group_by(condition) %>% 
               summarise(error = mean(error)) %>% 
               mutate(type = "MethanoGram (only mcrA)"),
             do.call(rbind, jackknife_res_rna) %>% 
               mutate(error = sqrt((truth - response)^2)) %>% 
               group_by(condition) %>% 
               summarise(error = mean(error)) %>% 
               mutate(type = "MethanoGram (only 16S rRNA)"),
             data.frame(condition = c("growth_doubl", "growth_rate", "mean_ogn", "mean_ogp", "mean_ogt"),
                        error = c(mean(median_pred(filter(read.csv("final_jackknife.csv"), 
                                                          condition == "growth_doubl")[["truth"]])),
                                  mean(median_pred(filter(read.csv("final_jackknife.csv"), 
                                                          condition == "growth_rate")[["truth"]])),
                                  mean(median_pred(filter(read.csv("final_jackknife.csv"), 
                                                          condition == "mean_ogn")[["truth"]])),
                                  mean(median_pred(filter(read.csv("final_jackknife.csv"), 
                                                          condition == "mean_ogp")[["truth"]])),
                                  mean(median_pred(filter(read.csv("final_jackknife.csv"), 
                                                          condition == "mean_ogt")[["truth"]]))),
                        type = "null model"),
             data.frame(condition = c("growth_doubl", "growth_rate", "mean_ogn", "mean_ogp", "mean_ogt"),
                        error = c(27.19, 0.35, 0.21, 0.47, 8.89),
                        type = "MethanoGram (previous model)")) %>% 
  left_join(read.csv("./data/full_names.csv"), by = c("condition" = "task.id")) %>% 
  filter(condition != "growth_rate", type != "MethanoGram (previous model)") %>% 
  droplevels() %>% 
  filter()

# read.csv("final_jackknife.csv") %>% 
#   mutate(error = sqrt((truth - response)^2)) %>% 
#   group_by(condition) %>% 
#   summarise(error_mean = mean(error)) %>% 
#   left_join(read.csv("./data/full_names.csv"), by = c("condition" = "task.id")) %>% 
#   select(nice, error_mean) %>% 
#   slice(-2) %>% 
#   write.csv("summary_jackknife.csv", row.names = FALSE)

appender <- function(string) {
  sapply(string, function(i)
    switch(i,
           `Growth doubling time [h]` = TeX("Growth doubling time $\\[h\\]$" ),
           `Optimal growth NaCl` = TeX("Optimal growth NaCl $\\[mol/dm$^3\\]$"),
           `Optimal growth pH` = TeX("Optimal growth pH"),
           `Optimal growth temp.` = TeX("Optimal growth temp.$\\[^{\\degree} C \\]$"))
  )
}

library(latex2exp)

p <- inner_join(filter(res, type %in% c("MethanoGram (only 16S rRNA)", "null model")) %>%
                  droplevels() %>% 
                  mutate(type = ifelse(type == "MethanoGram (only 16S rRNA)", "MethanoGram", type)),
                summarise_if(conditions_dat, is.numeric, .funs = c("sd")) %>% 
                  melt(variable.name = "condition",
                       value.name = "sd")) %>% 
  mutate(norm = error/sd,
         nice_label = paste0(round(error, 2), "\n(", round(norm, 2), ")")) %>% 
  ggplot(aes(x = type, y = error, fill = type, label = nice_label)) +
  geom_col() +
  scale_y_continuous("Mean error", expand = c(0.15, 0)) +
  geom_text(size = 3, aes(y = error)) +
  scale_x_discrete("") +
  scale_fill_manual("", values = c("cornflowerblue", "chocolate1")) +
  facet_wrap(~ nice, scales = "free_y", nrow = 1, 
             labeller = as_labeller(appender, default = label_parsed)) +
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom")

cairo_pdf(filename = "methanogram_jackknife.pdf", height = 3.5,
          width = 6.5, pointsize = 12, fallback_resolution = 500)
p
dev.off()


library(latex2exp)

rbind(mutate(best_pars, type = "16S rRNA and mcrA"),
      mutate(best_pars_rna, type = "Only 16S rRNA"),
      mutate(best_pars_mcra, type = "Only mcrA")) %>% 
  left_join(read.csv("./data/full_names.csv"), by = c("task.id" = "task.id")) %>% 
  filter(task.id != "growth_rate") %>% 
  droplevels %>% 
  ggplot(aes(x = type, y = mean_error, fill = type, label = round(mean_error, 4))) +
  geom_col() +
  scale_y_continuous(expand = c(0.15, 0)) +
  geom_text(size = 3, aes(y = 1.05 * mean_error)) +
  scale_x_discrete("") +
  scale_fill_discrete() +
  #coord_flip() +
  #scale_fill_manual("", values = c("cornflowerblue", "chocolate1")) +
  facet_wrap(~ nice, scales = "free_y", nrow = 1, 
             labeller = as_labeller(appender, default = label_parsed)) +
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom")

library(xtable)
rbind(mutate(best_pars, type = "Best predictor"),
      mutate(best_pars_rna, type = "Only 16S rRNA"),
      mutate(best_pars_mcra, type = "Only mcrA")) %>% 
  left_join(read.csv("./data/full_names.csv"), by = c("task.id" = "task.id")) %>% 
  filter(task.id != "growth_rate") %>% 
  droplevels %>% 
  select(task.id, type, mean_error) %>% 
  dcast(task.id ~ type) %>% 
  xtable(digits = 4)


do.call(rbind, jackknife_res_rna) %>% 
  mutate(error = sqrt((truth - response)^2)) %>% 
  group_by(condition) %>% 
  mutate(rank = rank(error)) %>% 
  filter(condition != "growth_rate") %>% 
  filter(rank < 10) %>% 
  ungroup %>% 
  select(Name) %>% 
  table %>% 
  sort
  
inner_join(filter(res, type == "MethanoGram (only 16S rRNA)") %>%
             droplevels() %>% 
             data.frame() %>% 
             select(condition, nice, error),
           summarise_if(conditions_dat, is.numeric, .funs = c("sd")) %>% 
             melt(variable.name = "condition",
                  value.name = "sd")) %>% 
  mutate(norm = error/sd) %>% 
  select(nice, error, norm)
