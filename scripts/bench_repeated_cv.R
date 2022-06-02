######################################################
# Simple Benchmark of 4 untuned survival learners on #
# the miRNA and mRNA tasks using repeated CV         #
######################################################
library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
library(dplyr)
library(stringr)
library(ggplot2)

# Logging ----
if (FALSE) { # logger doesn't work with parallelization enabled
  if (!dir.exists('log')) dir.create('log')

  logger = lgr::get_logger("mlr3")
  logger$set_threshold('info')
  filename = tempfile(pattern = 'bench_repeated_cv_', tmpdir = 'log', fileext = ".log")
  logger$add_appender(lgr::AppenderFile$new(filename), name = "logfile")
}

# Tasks ----
# first run `scripts/prepare_tasks.R`
tasks = readRDS(file = 'data/tasks.rds')
tasks = list(tasks$mRNA, tasks$miRNA) # use only mRNA and miRNA tasks

# Learners ----
# with default hyperparameters (no tuning)
glmnet_lrn  = lrn('surv.glmnet', id = 'CoxLasso', standardize = FALSE, lambda = 0.01, alpha = 1)
xgboost_lrn = lrn('surv.xgboost', id = 'XGBoost Survival Learner')
rpart_lrn   = lrn('surv.rpart', id = 'Survival Tree')
ranger_lrn  = lrn('surv.ranger', id = 'Survival Forest')

learners = list(glmnet_lrn, rpart_lrn, ranger_lrn, xgboost_lrn)

# Benchmark (Repeated-CV) ----
resampling = rsmp('repeated_cv', folds = 5, repeats = 10)
design = benchmark_grid(tasks, learners, resampling)

# parallelization + reproducibility
future::plan("multisession")
set.seed(42)

bm_res = benchmark(design)
#saveRDS(object = bm_res, file = 'results/bm_res_repeated_cv.rds')
