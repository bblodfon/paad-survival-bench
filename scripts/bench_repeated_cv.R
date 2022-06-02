#'####################################################
# Simple Benchmark of 4 untuned survival learners on #
# the miRNA and mRNA tasks using repeated CV         #
#'####################################################
library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
library(tictoc)

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

# Less logging
lgr::get_logger("bbotk")$set_threshold("warn")
lgr::get_logger("mlr3")$set_threshold("warn")

# parallelization + reproducibility
future::plan("multisession")
set.seed(42)

tic()
bm_res = benchmark(design, store_models = TRUE, store_backends = FALSE)
toc()

saveRDS(object = bm_res, file = 'results/bm_res_repeated_cv.rds')
