#'###############################################
# Benchmark xgboost on mRNA task using a simple #
# train/test split and Bayesian Optimization    #
#'###############################################
library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
library(mlr3mbo)
library(progressr)
library(tictoc)
source('scripts/boot_mlr3.R')

res_path = 'results/xgboost/'
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
num_threads = 50 # implicit parallelization
harrell_cindex = msr('surv.cindex')

# mRNA Task ----
# ~10000 features
tasks = readRDS(file = 'data/tasks.rds')
task_mRNA = tasks$mRNA

# ~70%/30% train/test split
train_indx = sample(seq_len(task_mRNA$nrow), 100)
test_indx  = setdiff(seq_len(task_mRNA$nrow), train_indx)

# Learners (4) ----
## XGBoost (survival-cox, 3 HPs) ----
learner1 = lrn('surv.xgboost',
  nthread = num_threads, booster = 'gbtree', early_stopping_rounds = 10,
  nrounds = to_tune(50, 5000),
  eta = to_tune(p_dbl(1e-04, 1, logscale = TRUE)),
  max_depth = to_tune(2, 10))

xgboost_at1 = AutoTuner$new(
  learner = learner1,
  resampling = rsmp('cv', folds = n_folds),
  measure = harrell_cindex,
  terminator = trm('evals', n_evals = n_evals),
  tuner = tnr('mbo')
)

print('Train Start: XGBoost (survival-cox, 3 HPs)')
tic()
xgboost_at1$train(task_mRNA, row_ids = train_indx)
toc()
print('Train Finished')

xgboost_at1$
  predict(task_mRNA, row_ids = test_indx)$
  score(harrell_cindex)

print('Calculate bootstrap CIs for C-index on test set')
boot_ci1 = boot_ci(data = task_mRNA$data()[test_indx], nrsmps = 1000,
  learner = xgboost_at1$learner, num_threads = 1)

print('Saving results...')
saveRDS(object = list(
  xgboost_at1 = xgboost_at1,
  boot_ci1 = boot_ci1,
  test_indx = test_indx,
  objective = 'survival-cox',
  hp_num = 3),
  file = paste0(res_path, 'xgboost1_res.rds'))

## XGBoost (survival-cox, 6 HPs) ----
learner2 = lrn('surv.xgboost',
  nthread = num_threads, booster = 'gbtree', early_stopping_rounds = 10,
  nrounds = to_tune(50, 5000),
  eta = to_tune(p_dbl(1e-04, 1, logscale = TRUE)),
  max_depth = to_tune(2, 10),
  min_child_weight = to_tune(1, 100, logscale = TRUE),
  alpha  = to_tune(1e-03, 10, logscale = TRUE),
  lambda = to_tune(1e-03, 10, logscale = TRUE))

xgboost_at2 = AutoTuner$new(
  learner = learner2,
  resampling = rsmp('cv', folds = n_folds),
  measure = harrell_cindex,
  terminator = trm('evals', n_evals = n_evals),
  tuner = tnr('mbo')
)

print('Train Start: XGBoost (survival-cox, 6 HPs)')
tic()
xgboost_at2$train(task_mRNA, row_ids = train_indx)
toc()
print('Train Finished')

xgboost_at2$
  predict(task_mRNA, row_ids = test_indx)$
  score(harrell_cindex)

print('Calculate bootstrap CIs for C-index on test set')
boot_ci2 = boot_ci(data = task_mRNA$data()[test_indx], nrsmps = 1000,
  learner = xgboost_at2$learner, num_threads = 1)

print('Saving results...')
saveRDS(object = list(
  xgboost_at2 = xgboost_at2,
  boot_ci2 = boot_ci2,
  objective = 'survival-cox',
  hp_num = 6),
  file = paste0(res_path, 'xgboost2_res.rds'))

## XGBoost (survival-aft, 5 HPs) ----
learner3 = lrn('surv.xgboost',
  nthread = num_threads, booster = 'gbtree', early_stopping_rounds = 10,
  nrounds = to_tune(50, 5000),
  eta = to_tune(p_dbl(1e-04, 1, logscale = TRUE)),
  max_depth = to_tune(2, 10),
  objective = 'survival:aft',
  aft_loss_distribution = to_tune(c('normal', 'logistic', 'extreme')),
  aft_loss_distribution_scale = to_tune(0.5, 2.0))

xgboost_at3 = AutoTuner$new(
  learner = learner3,
  resampling = rsmp('cv', folds = n_folds),
  measure = harrell_cindex,
  terminator = trm('evals', n_evals = n_evals),
  tuner = tnr('mbo')
)

print('Train Start: XGBoost (survival-aft, 5 HPs)')
tic()
xgboost_at3$train(task_mRNA, row_ids = train_indx)
toc()
print('Train Finished')

xgboost_at3$
  predict(task_mRNA, row_ids = test_indx)$
  score(harrell_cindex)

print('Calculate bootstrap CIs for C-index on test set')
boot_ci3 = boot_ci(data = task_mRNA$data()[test_indx], nrsmps = 1000,
  learner = xgboost_at3$learner, num_threads = 1)

print('Saving results...')
saveRDS(object = list(
  xgboost_at3 = xgboost_at3,
  boot_ci3 = boot_ci3,
  objective = 'survival-aft',
  hp_num = 5),
  file = paste0(res_path, 'xgboost3_res.rds'))

## XGBoost (survival-aft, 8 HPs) ----
learner4 = lrn('surv.xgboost',
  nthread = num_threads, booster = 'gbtree', early_stopping_rounds = 10,
  nrounds = to_tune(50, 5000),
  eta = to_tune(p_dbl(1e-04, 1, logscale = TRUE)),
  max_depth = to_tune(2, 10),
  min_child_weight = to_tune(1, 100, logscale = TRUE),
  alpha  = to_tune(1e-03, 10, logscale = TRUE),
  lambda = to_tune(1e-03, 10, logscale = TRUE),
  objective = 'survival:aft',
  aft_loss_distribution = to_tune(c('normal', 'logistic', 'extreme')),
  aft_loss_distribution_scale = to_tune(0.5, 2.0))

xgboost_at4 = AutoTuner$new(
  learner = learner4,
  resampling = rsmp('cv', folds = n_folds),
  measure = harrell_cindex,
  terminator = trm('evals', n_evals = n_evals),
  tuner = tnr('mbo')
)

print('Train Start: XGBoost (survival-aft, 8 HPs)')
tic()
xgboost_at4$train(task_mRNA, row_ids = train_indx)
toc()
print('Train Finished')

xgboost_at4$
  predict(task_mRNA, row_ids = test_indx)$
  score(harrell_cindex)

print('Calculate bootstrap CIs for C-index on test set')
boot_ci4 = boot_ci(data = task_mRNA$data()[test_indx], nrsmps = 1000,
  learner = xgboost_at4$learner, num_threads = 1)

print('Saving results...')
saveRDS(object = list(
  xgboost_at4 = xgboost_at4,
  boot_ci4 = boot_ci4,
  objective = 'survival-aft',
  hp_num = 8),
  file = paste0(res_path, 'xgboost4_res.rds'))
