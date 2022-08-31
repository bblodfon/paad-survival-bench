#'###############################################
# Benchmark xgboost on mRNA task using a simple #
# train/test split and Bayesian Optimization v2 #
#'###############################################
# Additions compared to `bench_xgboost.R`:
# 1. Tune `early_stopping_rounds` and `gamma`
# 2. Get best hpc (hyperparam config) + learner out of BO
# (previously it was best hpc as calculated by the surrogate)
# 3. Get C-index on train and test sets using every hpc
# (test for overfitting)

library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
library(mlr3mbo)
library(progressr)
library(tictoc)
suppressMessages(library(dplyr))
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
nfolds = 5
rsmp_cv = rsmp('cv', folds = nfolds)
harrell_cindex = msr('surv.cindex')
n_evals = 250 # number of hyperparameter configurations to search
eval_trm = trm('evals', n_evals = n_evals)
bayes_tnr = tnr('mbo')
num_threads = 50 # implicit parallelization

# mRNA Task ----
# ~10000 features
tasks = readRDS(file = 'data/tasks.rds')
task_mRNA = tasks$mRNA

# ~70%/30% train/test split
train_indx = sample(seq_len(task_mRNA$nrow), 100)
test_indx  = setdiff(seq_len(task_mRNA$nrow), train_indx)

# XGBoost Learners (2) ----
## Base learners ----
xgboost_cox = lrn('surv.xgboost', nthread = num_threads, booster = 'gbtree',
  fallback = lrn('surv.kaplan'))
xgboost_aft = lrn('surv.xgboost', nthread = num_threads, booster = 'gbtree',
  fallback = lrn('surv.kaplan'), objective = 'survival:aft')

## Parameter spaces ----
ps_cox = ps(
  nrounds = p_int(50, 5000),
  eta = p_dbl(1e-04, 1, logscale = TRUE),
  max_depth = p_int(2, 10),
  min_child_weight = p_dbl(1, 100, logscale = TRUE),
  alpha  = p_dbl(1e-03, 10, logscale = TRUE),
  lambda = p_dbl(1e-03, 10, logscale = TRUE),
  # tune gamma as well
  gamma  = p_dbl(0, 5),
  # early stopping => 10% of nrounds
  .extra_trafo = function(x, param_set) {
    x$early_stopping_rounds = as.integer(ceiling(0.1 * x$nrounds))
    x
  }
)

ps_aft = ps_cox$clone(deep = TRUE)
ps_aft$add(ps(
  aft_loss_distribution = p_fct(c('normal', 'logistic', 'extreme')),
  aft_loss_distribution_scale = p_dbl(0.5, 2.0))
)

# Train & Test benchmark ----
## XGBoost Cox ----
xgboost_at_cox = AutoTuner$new(
  learner = xgboost_cox,
  resampling = rsmp_cv,
  measure = harrell_cindex,
  search_space = ps_cox,
  terminator = eval_trm,
  tuner = bayes_tnr
)

print('Train Start: XGBoost (survival-cox, 7 HPs)')
tic()
xgboost_at_cox$train(task_mRNA, train_indx)
toc()
print('Train Finished')

print('Performance on test set of hpc that BO chose')
xgboost_at_cox$predict(task_mRNA, test_indx)$score(harrell_cindex)

# Get best C-index hpc and train learner with it (different from above)
best_hpc = xgboost_at_cox$archive$best()$x_domain[[1L]]
xgboost_cox$param_set$values =
  mlr3misc::insert_named(xgboost_cox$param_set$values, best_hpc)

xgboost_cox$train(task_mRNA, train_indx)

print('Performance on test set of best hpc')
xgboost_cox$predict(task_mRNA, test_indx)$score(harrell_cindex)

print('Calculate bootstrap CIs for C-index on test set (best hpc)')
tic()
boot_cox = boot_ci(data = task_mRNA$data()[test_indx], nrsmps = 1000,
  learner = xgboost_cox, num_threads = 1, parallel = 'no')
toc()

print('Calculate C-indexes for all hpcs (test + train set)')
# get all hpcs tested by BO
hpc_list = xgboost_at_cox$archive$data$x_domain
stopifnot(length(hpc_list) == n_evals)

# empty base xgboost learner
base_lrn = xgboost_cox$clone(deep = TRUE)
base_lrn$reset()

res = list()
for (i in 1:length(hpc_list)) {
  hpc = hpc_list[[i]]
  base_lrn$param_set$values =
    mlr3misc::insert_named(base_lrn$param_set$values, hpc)

  base_lrn$train(task_mRNA, train_indx)
  res[[i]] = dplyr::tibble(
    index = i,
    train_cindex = base_lrn$predict(task_mRNA, train_indx)$score(),
    test_cindex  = base_lrn$predict(task_mRNA, test_indx) $score()
  )

  base_lrn$reset()
}

hpc_res = dplyr::bind_rows(res)
# add the stored average CV C-indexes on the train data
hpc_res = hpc_res %>%
  mutate(traincv_cindex = xgboost_at_cox$archive$data$surv.cindex)

print('Saving results...')
saveRDS(object = list(
  train_indx = train_indx,
  test_indx  = test_indx,
  #xgboost_at_cox = xgboost_at_cox, # save space
  xgboost_cox = xgboost_cox, # trained learner with best hpc from BO
  boot_cox = boot_cox, # bootstrap results of best hpc learner on test set
  hpc_list = hpc_list,
  hpc_res  = hpc_res,
  ps_cox   = ps_cox,
  objective = 'survival-cox',
  hp_num = 7),
  file = paste0(res_path, 'xgboost_cox.rds'))

## XGBoost AFT ----
xgboost_at_aft = AutoTuner$new(
  learner = xgboost_aft,
  resampling = rsmp_cv,
  measure = harrell_cindex,
  search_space = ps_aft,
  terminator = eval_trm,
  tuner = bayes_tnr
)

print('Train Start: XGBoost (survival-aft, 9 HPs)')
tic()
xgboost_at_aft$train(task_mRNA, train_indx)
toc()
print('Train Finished')

print('Performance on test set of hpc that BO chose')
xgboost_at_aft$predict(task_mRNA, test_indx)$score(harrell_cindex)

# Get best C-index hpc and train learner with it (different from above)
best_hpc = xgboost_at_aft$archive$best()$x_domain[[1L]]
xgboost_aft$param_set$values =
  mlr3misc::insert_named(xgboost_aft$param_set$values, best_hpc)

xgboost_aft$train(task_mRNA, train_indx)

print('Performance on test set of best hpc')
xgboost_aft$predict(task_mRNA, test_indx)$score(harrell_cindex)

print('Calculate bootstrap CIs for C-index on test set (best hpc)')
tic()
boot_aft = boot_ci(data = task_mRNA$data()[test_indx], nrsmps = 1000,
  learner = xgboost_aft, num_threads = 1, parallel = 'no')
toc()

print('Calculate C-indexes for all hpcs (test + train set)')
# get all hpcs tested by BO
hpc_list = xgboost_at_aft$archive$data$x_domain
stopifnot(length(hpc_list) == n_evals)

# empty base xgboost learner
base_lrn = xgboost_aft$clone(deep = TRUE)
base_lrn$reset()

res = list()
for (i in 1:length(hpc_list)) {
  hpc = hpc_list[[i]]
  base_lrn$param_set$values =
    mlr3misc::insert_named(base_lrn$param_set$values, hpc)

  base_lrn$train(task_mRNA, train_indx)
  res[[i]] = dplyr::tibble(
    index = i,
    train_cindex = base_lrn$predict(task_mRNA, train_indx)$score(),
    test_cindex  = base_lrn$predict(task_mRNA, test_indx) $score()
  )

  base_lrn$reset()
}

hpc_res = dplyr::bind_rows(res)
# add the stored average CV C-indexes on the train data
hpc_res = hpc_res %>%
  mutate(traincv_cindex = xgboost_at_aft$archive$data$surv.cindex)

print('Saving results...')
saveRDS(object = list(
  train_indx = train_indx,
  test_indx  = test_indx,
  #xgboost_at_aft = xgboost_at_aft, # save space
  xgboost_aft = xgboost_aft, # trained learner with best hpc from BO
  boot_aft = boot_aft, # bootstrap results of best hpc learner on test set
  hpc_list = hpc_list,
  hpc_res  = hpc_res,
  ps_aft   = ps_aft,
  objective = 'survival-aft',
  hp_num = 9),
  file = paste0(res_path, 'xgboost_aft.rds'))
