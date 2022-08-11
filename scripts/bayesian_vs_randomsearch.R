#'##################################################
# Compare Bayesian Optimization vs Random Search   #
# strategies for hyperparameter tuning (mRNA data) #
#'##################################################
library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
library(mlr3mbo)
library(progressr)
library(ggplot2)
library(boot)

res_path = 'results/bayes_vs_randsearch/'
if (!dir.exists(res_path)) {
  dir.create(res_path)
}

# for reproducibility
set.seed(42)

# Add progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers('progress')

# Less logging
lgr::get_logger('bbotk')$set_threshold('warn')
lgr::get_logger('mlr3')$set_threshold('warn')

# Global variables ----
nfolds = 5
rsmp_cv = rsmp('cv', folds = nfolds)
harrell_cindex = msr('surv.cindex')
n_evals = 250 # number of hyperparameter configurations to search
eval_trm = trm('evals', n_evals = n_evals)
rand_tnr = tnr('random_search', batch_size = 1) # no parallelization
bayes_tnr = tnr('mbo')
num_threads = 20 # implicit parallelization for RFs

# Bootstrap functions ----

#' `data` is a data.table/data.frame object with the test data
#' and has to have the same structure (features and target columns)
#' as the task that was used to train the learner
boot_fun = function(data, index, learner, measure) {
  learner$predict_newdata(data[index])$score(measure)
}

get_bootCIs = function(data, learner, measure = mlr3::msr('surv.cindex'),
  num_threads = 1, nrsmps = 1000) {
  boot_res = boot::boot(data, statistic = boot_fun, R = nrsmps,
    parallel = 'multicore', ncpus = num_threads, learner = learner, measure = measure)

  bootci_res = boot::boot.ci(boot_res, type = c('basic', 'norm', 'perc', 'bca'))

  return(list(boot_res = boot_res, bootci_res = bootci_res))
}

# Tasks (1) ----
# mRNA (~10000 features)
tasks = readRDS(file = 'data/tasks.rds')
task_mRNA = tasks$mRNA

# ~70%/30% train/test split
train_indx = sample(seq_len(task_mRNA$nrow), 100)
test_indx  = setdiff(seq_len(task_mRNA$nrow), train_indx)
#intersect(train_indx, test_indx)

# Learners (2) ----
## CoxNet (2 HPs) ----
coxnet_lrn = lrn('surv.glmnet', id = 'CoxNet',
  standardize = FALSE, maxit = 10^3,
  lambda = to_tune(p_dbl(1e-04, 1, logscale = TRUE)),
  alpha = to_tune(0, 1)) # from Ridge to Lasso penalty

coxnet_bayes_at = AutoTuner$new(
  learner = coxnet_lrn,
  resampling = rsmp_cv,
  measure = harrell_cindex,
  terminator = eval_trm,
  tuner = bayes_tnr)

coxnet_rand_at = AutoTuner$new(
  learner = coxnet_lrn,
  resampling = rsmp_cv,
  measure = harrell_cindex,
  terminator = eval_trm,
  tuner = rand_tnr)

## Survival Forests (4 HPs) ----
ranger_lrn = lrn('surv.ranger', verbose = FALSE, id = 'SurvivalForest',
  num.trees = to_tune(100, 1500),
  mtry.ratio = to_tune(0.01, 0.1), # e.g. choose between ~100 and 1000 mRNA features
  min.node.size = to_tune(1, 20),
  splitrule = to_tune(c('logrank', 'maxstat', 'C')),
  num.threads = num_threads)

ranger_bayes_at = AutoTuner$new(
  learner = ranger_lrn,
  resampling = rsmp_cv,
  measure = harrell_cindex,
  terminator = eval_trm,
  tuner = bayes_tnr)

ranger_rand_at = AutoTuner$new(
  learner = ranger_lrn,
  resampling = rsmp_cv,
  measure = harrell_cindex,
  terminator = eval_trm,
  tuner = rand_tnr)

# Benchmark ----
## CoxNet ----
print('Start: CoxNet (BO)')
coxnet_bayes_at$train(task_mRNA, row_ids = train_indx)
coxnet_bayes_at$
  predict(task_mRNA, row_ids = test_indx)$
  score(harrell_cindex)
print('Finish: CoxNet (BO)')

print('Start: CoxNet (RS)')
coxnet_rand_at$train(task_mRNA, row_ids = train_indx)
coxnet_rand_at$
  predict(task_mRNA, row_ids = test_indx)$
  score(harrell_cindex)
print('Finish: CoxNet (RS)')

print('Save tuning results')
saveRDS(object = list(
  coxnet_bayes_at = coxnet_bayes_at,
  coxnet_rand_at  = coxnet_rand_at),
  file = paste0(res_path, 'coxnet_tuning_res.rds'))

print('Calculate bootstrap CIs for C-index on test set')
bootci_coxnet_bayes = get_bootCIs(data = task_mRNA$data()[test_indx], nrsmps = 1000,
  learner = coxnet_bayes_at$learner, num_threads = num_threads)
bootci_coxnet_rand  = get_bootCIs(data = task_mRNA$data()[test_indx], nrsmps = 1000,
  learner = coxnet_rand_at$learner, num_threads = num_threads)

saveRDS(object = list(
  bootci_coxnet_bayes = bootci_coxnet_bayes,
  bootci_coxnet_rand  = bootci_coxnet_rand),
  file = paste0(res_path, 'coxnet_bootci.rds'))

## Survival Forests ----
print('Start: Survival Forest (BO)')
ranger_bayes_at$train(task_mRNA, row_ids = train_indx)
ranger_bayes_at$
  predict(task_mRNA, row_ids = test_indx)$
  score(harrell_cindex)
print('Finish: Survival Forest (BO)')

print('Start: Survival Forest (RS)')
ranger_rand_at$train(task_mRNA, row_ids = train_indx)
ranger_rand_at$
  predict(task_mRNA, row_ids = test_indx)$
  score(harrell_cindex)
print('Finish: Survival Forest (RS)')

print('Save tuning results')
saveRDS(object = list(
  ranger_bayes_at = ranger_bayes_at,
  ranger_rand_at  = ranger_rand_at),
  file = paste0(res_path, 'ranger_tuning_res.rds'))

print('Calculate bootstrap CIs for C-index on test set')
bootci_ranger_bayes = get_bootCIs(data = task_mRNA$data()[test_indx], nrsmps = 1000,
  learner = ranger_bayes_at$learner, num_threads = num_threads)
bootci_ranger_rand  = get_bootCIs(data = task_mRNA$data()[test_indx], nrsmps = 1000,
  learner = ranger_rand_at$learner, num_threads = num_threads)
saveRDS(object = list(
  bootci_ranger_bayes = bootci_ranger_bayes,
  bootci_ranger_rand  = bootci_ranger_rand),
  file = paste0(res_path, 'ranger_bootci.rds'))
