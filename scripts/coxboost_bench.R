#'#############################################
#' Benchmark
#' - CoxBoost using mRNA only vs mRNA + clinical
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

res_path = 'results/coxboost/'
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

# Tasks ----
## mRNA (~10000 features) ----
tasks = readRDS(file = 'data/tasks.rds')
task_mRNA = tasks$mRNA

## ~70%/30% train/test split
train_indx = sample(seq_len(task_mRNA$nrow), 100)
test_indx  = setdiff(seq_len(task_mRNA$nrow), train_indx)

## Clinical (8 features) ----
task_clinical = readRDS(file = 'data/task_clinical.rds')

## Clinical + mRNA ----
# make combined task with clinical first and then mRNA features
cmb_task = task_clinical$clone(deep = TRUE)
cmb_task$cbind(task_mRNA$data(cols = task_mRNA$feature_names))

# Benchmarks ----
## CoxBoost (only mRNA) ----
coxboost_lrn = lrn('surv.coxboost',
  standardize = FALSE, # data already standardized
  return.score = FALSE, # don't need this in the output
  fallback = lrn('surv.kaplan'),
  stepno = to_tune(p_int(50, 500)),
  penalty = to_tune(p_int(10, 1000, logscale = TRUE)), # leave at default => 9 * sum(status == 1)?
  stepsize.factor = to_tune(p_dbl(1e-01, 10, logscale = TRUE)), # leave at default => 1?
  criterion = to_tune(p_fct(c('pscore', 'hpscore'))) # penalized scores, `hpscore` is faster to calculate
)

coxboost_at = AutoTuner$new(
  learner = coxboost_lrn,
  resampling = rsmp('cv', folds = n_folds),
  measure = harrell_cindex,
  terminator = trm('evals', n_evals = n_evals),
  tuner = tnr('mbo')
)

print('Train Start: CoxBoost (mRNA only)')
tic()
coxboost_at$train(task_mRNA, train_indx)
toc()
print('Train Finished')

print('Performance on test set of hpc that BO chose')
coxboost_at$predict(task_mRNA, test_indx)$score(harrell_cindex)

print('Performance on test set of best hpc')
best_hpc = coxboost_at$archive$best()$x_domain[[1L]]
coxboost_lrn$param_set$values = mlr3misc::insert_named(coxboost_lrn$param_set$values, best_hpc)
coxboost_lrn$train(task_mRNA, train_indx)
coxboost_lrn$predict(task_mRNA, test_indx)$score(harrell_cindex)

print('Calculate bootstrap CIs for C-index on test set (best hpc)')
tic()
boot_cox = boot_ci(data = task_mRNA$data(rows = test_indx), nrsmps = 1000,
  learner = coxboost_lrn, num_threads = 1, parallel = 'no')
toc()

print('Calculate C-indexes for all hpcs (test + train set)')
tic()
hpc_res = get_cindex_all_hps(at = coxboost_at, task = task_mRNA, train_indx = train_indx)
toc()

print('Saving results...')
saveRDS(object = list(
  train_indx = train_indx,
  test_indx  = test_indx,
  coxboost_lrn = coxboost_lrn, # trained learner with best hpc from BO
  boot_cox = boot_cox, # bootstrap test set results of learner with best hpc
  hpc_res  = hpc_res),
  file = paste0(res_path, 'coxboost_mRNA_test.rds')
)

## CoxBoost (clinical + mRNA) ----
coxboost_lrn = lrn('surv.coxboost',
  standardize = FALSE, # data already standardized
  return.score = FALSE, # don't need this in the output
  fallback = lrn('surv.kaplan'),
  stepno = to_tune(p_int(50, 500)),
  penalty = to_tune(p_int(10, 1000, logscale = TRUE)), # leave at default => 9 * sum(status == 1)?
  stepsize.factor = to_tune(p_dbl(1e-01, 10, logscale = TRUE)), # leave at default => 1?
  criterion = to_tune(p_fct(c('pscore', 'hpscore'))), # penalized scores, `hpscore` is faster to calculate
  unpen.index = 1:length(task_clinical$feature_names) # add clinical features as mandatory
)

coxboost_at = AutoTuner$new(
  learner = coxboost_lrn,
  resampling = rsmp('cv', folds = n_folds),
  measure = harrell_cindex,
  terminator = trm('evals', n_evals = n_evals),
  tuner = tnr('mbo')
)

print('Train Start: CoxBoost (mRNA + clinical)')
tic()
coxboost_at$train(cmb_task, train_indx)
toc()
print('Train Finished')

print('Performance on test set of hpc that BO chose')
coxboost_at$predict(cmb_task, test_indx)$score(harrell_cindex)

print('Performance on test set of best hpc')
best_hpc = coxboost_at$archive$best()$x_domain[[1L]]
coxboost_lrn$param_set$values = mlr3misc::insert_named(coxboost_lrn$param_set$values, best_hpc)
coxboost_lrn$train(cmb_task, train_indx)
coxboost_lrn$predict(cmb_task, test_indx)$score(harrell_cindex)

print('Calculate bootstrap CIs for C-index on test set (best hpc)')
tic()
boot_cox = boot_ci(data = cmb_task$data(rows = test_indx), nrsmps = 1000,
  learner = coxboost_lrn, num_threads = 1, parallel = 'no')
toc()

print('Calculate C-indexes for all hpcs (test + train set)')
tic()
hpc_res = get_cindex_all_hps(at = coxboost_at, task = cmb_task, train_indx = train_indx)
toc()

print('Saving results...')
saveRDS(object = list(
  train_indx = train_indx,
  test_indx  = test_indx,
  coxboost_lrn = coxboost_lrn, # trained learner with best hpc from BO
  boot_cox = boot_cox, # bootstrap test set results of learner with best hpc
  hpc_res  = hpc_res),
  file = paste0(res_path, 'coxboost_clinical_mRNA_test.rds')
)
