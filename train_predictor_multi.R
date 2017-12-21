library(dplyr)
library(biogram)
library(seqinr)
library(pbapply)
library(mlr)

configureMlr(show.info = FALSE)

conditions_dat <- read.csv("./data/condition_data.csv")

load("./data/seq_dat.RData")

training_data <- expand.grid(type1 = c("mcra", "rna"), len1 = 1L:6, 
                             type2 = c("mcra", "rna"), len2 = 1L:6) %>% 
  split(1L:nrow(.)) 

ngram_dat_list <- lapply(training_data, function(i) {
  if(i[["type1"]] == i[["type2"]] & i[["len1"]] == i[["len2"]]) {
    seq_dat[[i[["type1"]]]][[i[["len1"]]]]
  } else {
    cbind(select(seq_dat[[i[["type1"]]]][[i[["len1"]]]], -species),
          seq_dat[[i[["type2"]]]][[i[["len2"]]]])
  }
})

library("parallelMap")

options(parallelMap.default.mode = "multicore",
        parallelMap.default.cpus = 4,
        parallelMap.default.show.info = FALSE)

parallelStartSocket(4)

benchmark_res <- pblapply(c("growth_doubl", "growth_rate", "mean_ogt", 
                            "mean_ogn", "mean_ogp"), function(ith_condition)
                              lapply(1L:length(ngram_dat_list), function(ith_ngram_dat_list) {
                                try({
                                  ngram_dat <- ngram_dat_list[[ith_ngram_dat_list]]
                                  
                                  dat <- conditions_dat[, c("Name", ith_condition)] %>% 
                                    inner_join(ngram_dat, by = c("Name" = "species")) %>% 
                                    select(-Name)
                                  
                                  predict_ngrams <- makeRegrTask(id = ith_condition, 
                                                                 data = dat, 
                                                                 target = ith_condition)
                                  
                                  learnerRF <- makeFilterWrapper(learner = makeLearner("regr.ranger"), 
                                                                 fw.method = "linear.correlation")
                                  
                                  learner_pars <- makeParamSet(
                                    makeDiscreteParam("num.trees", values = c(250, 500, 750, 1000)),
                                    makeDiscreteParam("min.node.size", values = c(3, 5, 7)),
                                    makeDiscreteParam("fw.perc", values = c(0.1, 0.25, 0.5))
                                  )
                                  
                                  set.seed(1410)
                                  
                                  inner <- makeResampleDesc("CV", iters = 5)
                                  outer <- makeResampleDesc("CV", iters = 3)
                                  learnerRF_tuned <- makeTuneWrapper(learnerRF, 
                                                                     resampling = inner, 
                                                                     par.set = learner_pars, 
                                                                     control = makeTuneControlGrid())
                                  
                                  nested_cv <- resample(learnerRF_tuned, predict_ngrams, outer, extract = getTuneResult)
                                  
                                  nested_res <- getNestedTuneResultsOptPathDf(nested_cv) 
                                  
                                  group_by(nested_res, num.trees, min.node.size, fw.perc) %>% 
                                    summarise(mean_error = mean(mse.test.mean),
                                              sd_error = sd(mse.test.mean)) %>% 
                                    mutate(task.id = ith_condition) %>% 
                                    mutate(dat_list = ith_ngram_dat_list)
                                }, silent = TRUE)
                              })
)

parallelStop()

save(benchmark_res, file = "./results/ngram_benchmark_full.RData")
