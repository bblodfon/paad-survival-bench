#'######################################################
# Benchmark CoxNet, Survival Tree and Survival Forests #
# Same as `bench_nestedCV_v2`, but now parallelized    #
#'######################################################
library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
library(mlr3misc)
suppressMessages(library(dplyr))
library(progressr)
library(tictoc)

# for reproducibility
set.seed(42)

# Add progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers('progress')

# Less logging
lgr::get_logger('bbotk')$set_threshold('warn')
lgr::get_logger('mlr3')$set_threshold('warn')

# Tasks (3) ----
# first run `scripts/prepare_tasks.R`
tasks = readRDS(file = 'data/tasks.rds')
tasks = list(tasks$mRNA, tasks$CNA, tasks$`mRNA-CNA`)

# Global variables ----
#baseline_estimator = 'nelson' # Nelson-Aalen estimator
#distr_form = 'aft' # Accelerated failure time
out_folds = 4
in_folds  = 5
out_rsmp  = rsmp('cv', folds = out_folds)
in_rsmp   = rsmp('cv', folds = in_folds)

harrell_cindex = msr('surv.cindex')
n_evals = 100 # number of hyperparameter configurations to search
eval_trm = trm('evals', n_evals = n_evals)
batch_size_parallel = 6 # `batch_size_parallel` x `in_folds` = #CPUs
rand_tnr_seq   = tnr('random_search', batch_size = 1) # tuner strategy
rand_tnr_multi = tnr('random_search', batch_size = batch_size_parallel)
total_nthreads = 30 # implicit parallelization

# Learners (3) ----
## CoxNet ----
coxnet_lrn = lrn('surv.glmnet', id = 'CoxNet',
  standardize = FALSE, maxit = 10^3,
  lambda = to_tune(p_dbl(1e-03, 1, logscale = TRUE)),
  alpha = to_tune(0, 1)) # from Ridge to Lasso penalty

coxnet_at = AutoTuner$new(
  learner = coxnet_lrn,
  resampling = in_rsmp,
  measure = harrell_cindex,
  terminator = eval_trm,
  tuner = rand_tnr_multi)

## Survival Tree ----
rpart_lrn = lrn('surv.rpart', id = 'SurvivalTree',
  minsplit = to_tune(1, 20),
  cp = to_tune(1e-04, 0.1, logscale = TRUE))

rpart_at = AutoTuner$new(
  learner = rpart_lrn,
  resampling = in_rsmp,
  measure = harrell_cindex,
  terminator = eval_trm,
  tuner = rand_tnr_multi)

## Survival Forest ----
ranger_lrn = lrn('surv.ranger', verbose = FALSE, id = 'SurvivalForest',
  num.trees = to_tune(100, 1500),
  mtry.ratio = to_tune(0.01, 0.1), # e.g. choose between ~100 and 1000 mRNA features
  min.node.size = to_tune(1, 20),
  # include 'C' splitrule for maybe a bit better performance but much more compute time
  splitrule = to_tune(c('logrank', 'maxstat')),
  num.threads = total_nthreads)

ranger_at = AutoTuner$new(
  learner = ranger_lrn,
  resampling = in_rsmp,
  measure = harrell_cindex,
  terminator = eval_trm,
  tuner = rand_tnr_seq)

learners = list(coxnet_at, rpart_at, ranger_at)

# Benchmark ----
grid = benchmark_grid(tasks, learners, out_rsmp) # same resamplings

# subgrids
coxnet_grid = grid[mlr3misc::map_chr(grid$learner, `[[`, 'id') == 'CoxNet.tuned']
rpart_grid  = grid[mlr3misc::map_chr(grid$learner, `[[`, 'id') == 'SurvivalTree.tuned']
ranger_grid = grid[mlr3misc::map_chr(grid$learner, `[[`, 'id') == 'SurvivalForest.tuned']

# CoxNet and Survival Tree use explicit (inner fold/tuning) parallelization
future::plan(list('sequential', 'multisession'))
#future::plan(list('sequential', future::tweak('multisession', workers = total_nthreads)))

print('Start: CoxNet')
tic()
coxnet_bm = benchmark(coxnet_grid, store_models = FALSE, store_backends = FALSE)
toc()
print('Finish: CoxNet')

print('Start: Survival Tree')
tic()
rpart_bm = benchmark(rpart_grid, store_models = FALSE, store_backends = FALSE)
toc()
print('Finish: Survival Tree')

# Survival Forest uses implicit parallelization
future::plan('sequential')

print('Start: Survival Forest')
tic()
ranger_bm = benchmark(ranger_grid, store_models = FALSE, store_backends = FALSE)
toc()
print('Finish: Survival Forest')

# combine benchmark results
bm_res = coxnet_bm$combine(rpart_bm)$combine(ranger_bm)

# save prediction performance results
perf_res = bm_res$score(harrell_cindex) %>%
  as_tibble() %>%
  select(task_id, learner_id, surv.cindex)
saveRDS(perf_res, file = 'results/perf_res_nestedCV_v4.rds')
