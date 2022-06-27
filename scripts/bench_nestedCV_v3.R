# Testing several parallelization configurations for 3 survival learners
# (CoxNet, Survival Tree, Survival Forest), using only the mRNA dataset
library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
suppressMessages(library(dplyr))
library(progressr)
library(tictoc)

# Add progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers("progress")

# Tasks ----
# first run `scripts/prepare_tasks.R`
tasks = readRDS(file = 'data/tasks.rds')
mRNA_task = tasks$mRNA # (145, 10467)

# Global variables ----
#baseline_estimator = 'nelson' # Nelson-Aalen estimator
#distr_form = 'aft' # Accelerated failure time
out_folds = 4
in_folds  = 5
out_rsmp  = rsmp('cv', folds = out_folds)
in_rsmp   = rsmp('cv', folds = in_folds)

# for reproducibility: all learners train, tune and test on the same subtasks
set.seed(42)
out_rsmp$instantiate(mRNA_task)

harrell_cindex = msr('surv.cindex')
n_evals = 10 # number of hyperparameter configurations searched
eval_trm = trm("evals", n_evals = n_evals)
rand_tnr = tnr('random_search', batch_size = 1) # tuner strategy

# Less logging
lgr::get_logger("bbotk")$set_threshold("warn")
lgr::get_logger("mlr3")$set_threshold("warn")

# CoxNet ----
coxnet_lrn = lrn('surv.glmnet', standardize = FALSE, maxit = 10^3,
  lambda = to_tune(p_dbl(1e-03, 1, logscale = TRUE)),
  alpha = to_tune(0, 1)) # from Ridge to Lasso penalty

coxnet_at = AutoTuner$new(
  learner = coxnet_lrn,
  resampling = in_rsmp,
  measure = harrell_cindex,
  terminator = eval_trm,
  tuner = rand_tnr)

# tic()
# coxnet_rr = resample(mRNA_task, coxnet_at, out_rsmp)
# toc()

# Survival Tree learner ----
rpart_lrn = lrn('surv.rpart', id = 'SurvivalTree',
  minsplit = to_tune(1, 20),
  cp = to_tune(1e-04, 0.1, logscale = TRUE))

rpart_at = AutoTuner$new(
  learner = rpart_lrn,
  resampling = in_rsmp,
  measure = harrell_cindex,
  terminator = eval_trm,
  tuner = rand_tnr)

# tic()
# rpart_rr = resample(mRNA_task, rpart_at, out_rsmp)
# toc()

# check if benchmark behaves the same in terms of CPU utilization/parallelization (it does!!!)
# design = benchmark_grid(tasks = list(mRNA_task), learners = list(rpart_at), out_rsmp)
# tic()
# benchmark(design)
# toc()

# Survival Forest learner (logrank) ----
ranger_lrn = lrn('surv.ranger', verbose = FALSE, id = 'SurvivalForest',
  num.trees = to_tune(100, 1500),
  mtry.ratio = to_tune(0.01, 0.1), # choose between ~100 and 1000 mRNA features
  min.node.size = to_tune(1, 20),
  num.threads = 4,
  splitrule = 'logrank')

ranger_at = AutoTuner$new(
  learner = ranger_lrn,
  resampling = in_rsmp,
  measure = harrell_cindex,
  terminator = eval_trm,
  tuner = rand_tnr)

print('RSF (logrank)')
tic()
ranger_rr_logrank = resample(mRNA_task, ranger_at, out_rsmp)
toc()

# Survival Forest learner (maxstat) ----
ranger_lrn = lrn('surv.ranger', verbose = FALSE, id = 'SurvivalForest',
  num.trees = to_tune(100, 1500),
  mtry.ratio = to_tune(0.01, 0.1), # choose between ~100 and 1000 mRNA features
  min.node.size = to_tune(1, 20),
  num.threads = 30,
  splitrule = 'maxstat')

ranger_at = AutoTuner$new(
  learner = ranger_lrn,
  resampling = in_rsmp,
  measure = harrell_cindex,
  terminator = eval_trm,
  tuner = rand_tnr)

#future::plan(list("sequential", "multisession"))
print('RSF (maxstat)')
tic()
ranger_rr_maxstat = resample(mRNA_task, ranger_at, out_rsmp)
toc()

saveRDS(list(logrank = ranger_rr_logrank, maxstat = ranger_rr_maxstat), file = 'results/rsf_seq_implicit_parallelization.rds')
