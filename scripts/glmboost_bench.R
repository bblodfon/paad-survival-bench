#'#############################################
#' Benchmark mboost::glmboost() using mRNA data
#' Dataset: mRNA, with simple train/test split
#'          + bootstrap results on test set
#' Tuning: Bayesian Optimization
#'#############################################
library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
library(mlr3mbo)
library(progressr)
library(tictoc)
source('scripts/boot_mlr3.R')
source('scripts/helpers.R')

res_path = 'results/glmboost/'
if (!dir.exists(res_path)) {
  dir.create(res_path)
}

# Reproducibility
set.seed(42)

# Progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers('progress')

# Logging (less)
lgr::get_logger('bbotk')$set_threshold('warn')
lgr::get_logger('mlr3')$set_threshold('warn')

# Global variables ----
n_folds = 5 # CV folds
n_evals = 250 # number of hyperparameter configurations to search
harrell_cindex = msr('surv.cindex')

# Task mRNA ----
## ~10000 features
tasks = readRDS(file = 'data/tasks.rds')
task_mRNA = tasks$mRNA

## ~70%/30% train/test split
train_indx = sample(seq_len(task_mRNA$nrow), 100)
test_indx  = setdiff(seq_len(task_mRNA$nrow), train_indx)

# Glmboost ----
## CoxPH ----
glmboost_lrn = lrn('surv.glmboost', family = 'coxph', center = FALSE,
  mstop = to_tune(p_int(50, 500)),
  nu = to_tune(p_dbl(1e-03, 1, logscale = TRUE)),
  fallback = lrn('surv.kaplan')
)

glmboost_at = AutoTuner$new(
  learner = glmboost_lrn,
  resampling = rsmp('cv', folds = n_folds),
  measure = harrell_cindex,
  terminator = trm('evals', n_evals = n_evals),
  tuner = tnr('mbo')
)

print('Train Start: Glmboost (CoxPH family, mRNA)')
tic()
glmboost_at$train(task_mRNA, train_indx)
toc()
print('Train Finished')

print('Performance on test set of hpc that BO chose')
glmboost_at$predict(task_mRNA, test_indx)$score(harrell_cindex)

print('Performance on test set of best hpc')
best_hpc = glmboost_at$archive$best()$x_domain[[1L]]
glmboost_lrn$param_set$values = mlr3misc::insert_named(glmboost_lrn$param_set$values, best_hpc)
glmboost_lrn$train(task_mRNA, train_indx)
glmboost_lrn$predict(task_mRNA, test_indx)$score(harrell_cindex)

print('Calculate bootstrap CIs for C-index on test set (best hpc)')
tic()
boot_cox = boot_ci(data = task_mRNA$data(rows = test_indx), nrsmps = 1000,
  learner = glmboost_lrn, num_threads = 1, parallel = 'no')
toc()

print('Calculate C-indexes for all hpcs (test + train set)')
tic()
hpc_res = get_cindex_all_hps(at = glmboost_at, task = task_mRNA, train_indx = train_indx)
toc()

print('Saving results...')
saveRDS(object = list(
  train_indx = train_indx,
  test_indx  = test_indx,
  glmboost_lrn = glmboost_lrn, # trained learner with best hpc from BO
  boot_cox = boot_cox, # bootstrap test set results of learner with best hpc
  hpc_res = hpc_res),
  file = paste0(res_path, 'glmboost_CoxPH_mRNA.rds')
)

## AFT ----
glmboost_lrn = lrn('surv.glmboost', center = FALSE,
  family = to_tune(p_fct(c('weibull', 'loglog', 'lognormal', 'gehan'))),
  mstop = to_tune(p_int(50, 500)),
  nu = to_tune(p_dbl(1e-03, 1, logscale = TRUE)),
  fallback = lrn('surv.kaplan')
)

glmboost_at = AutoTuner$new(
  learner = glmboost_lrn,
  resampling = rsmp('cv', folds = n_folds),
  measure = harrell_cindex,
  terminator = trm('evals', n_evals = n_evals),
  tuner = tnr('mbo')
)

print('Train Start: Glmboost (AFT family, mRNA)')
tic()
glmboost_at$train(task_mRNA, train_indx)
toc()
print('Train Finished')

print('Performance on test set of hpc that BO chose')
glmboost_at$predict(task_mRNA, test_indx)$score(harrell_cindex)

print('Performance on test set of best hpc')
best_hpc = glmboost_at$archive$best()$x_domain[[1L]]
glmboost_lrn$param_set$values = mlr3misc::insert_named(glmboost_lrn$param_set$values, best_hpc)
glmboost_lrn$train(task_mRNA, train_indx)
glmboost_lrn$predict(task_mRNA, test_indx)$score(harrell_cindex)

print('Calculate bootstrap CIs for C-index on test set (best hpc)')
tic()
boot_aft = boot_ci(data = task_mRNA$data(rows = test_indx), nrsmps = 1000,
  learner = glmboost_lrn, num_threads = 1, parallel = 'no')
toc()

print('Calculate C-indexes for all hpcs (test + train set)')
tic()
hpc_res = get_cindex_all_hps(at = glmboost_at, task = task_mRNA, train_indx = train_indx)
toc()

print('Saving results...')
saveRDS(object = list(
  train_indx = train_indx,
  test_indx  = test_indx,
  glmboost_lrn = glmboost_lrn, # trained learner with best hpc from BO
  boot_aft = boot_aft, # bootstrap test set results of learner with best hpc
  hpc_res = hpc_res),
  file = paste0(res_path, 'glmboost_AFT_mRNA.rds')
)

## Cindex (Uno) family ----
glmboost_lrn = lrn('surv.glmboost', center = FALSE, family = 'cindex',
  sigma = to_tune(p_dbl(0.1, 0.5)), # smoothing hp for sigmoid
  mstop = to_tune(p_int(50, 500)),
  nu = to_tune(p_dbl(1e-03, 1, logscale = TRUE)),
  fallback = lrn('surv.kaplan')
)

glmboost_at = AutoTuner$new(
  learner = glmboost_lrn,
  resampling = rsmp('cv', folds = n_folds),
  measure = harrell_cindex,
  terminator = trm('evals', n_evals = n_evals),
  tuner = tnr('mbo')
)

print('Train Start: Glmboost (C-index family, mRNA)')
tic()
glmboost_at$train(task_mRNA, train_indx)
toc()
print('Train Finished')

print('Performance on test set of hpc that BO chose')
glmboost_at$predict(task_mRNA, test_indx)$score(harrell_cindex)

print('Performance on test set of best hpc')
best_hpc = glmboost_at$archive$best()$x_domain[[1L]]
glmboost_lrn$param_set$values = mlr3misc::insert_named(glmboost_lrn$param_set$values, best_hpc)
glmboost_lrn$train(task_mRNA, train_indx)
glmboost_lrn$predict(task_mRNA, test_indx)$score(harrell_cindex)

print('Calculate bootstrap CIs for C-index on test set (best hpc)')
tic()
boot_cindex = boot_ci(data = task_mRNA$data(rows = test_indx), nrsmps = 1000,
  learner = glmboost_lrn, num_threads = 1, parallel = 'no')
toc()

print('Calculate C-indexes for all hpcs (test + train set)')
tic()
hpc_res = get_cindex_all_hps(at = glmboost_at, task = task_mRNA, train_indx = train_indx)
toc()

print('Saving results...')
saveRDS(object = list(
  train_indx = train_indx,
  test_indx  = test_indx,
  glmboost_lrn = glmboost_lrn, # trained learner with best hpc from BO
  boot_cindex = boot_cindex, # bootstrap test set results of learner with best hpc
  hpc_res = hpc_res),
  file = paste0(res_path, 'glmboost_Cindex_mRNA.rds')
)
