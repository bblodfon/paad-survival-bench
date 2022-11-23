#'##################################################
#' Benchmark
#' - Use some ML survival learners
#' - Task powerset (all combinations), after FS has
#' been applied (so less features in total)
#' - Simple train/test split + bootstrap C-index
#' on test set (based on mRNA data split benchmarks)
#' - Bayesian Optimization for hyperparameter tuning
#'##################################################
library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
library(mlr3mbo)
library(progressr)
library(tictoc)
source('scripts/helpers.R')

res_path = 'results/powerset_bench/'
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
config = list(
  nfolds = 5, # for (repeated) CV resampling
  repeats = 1, # for (repeated) CV resampling
  nevals = 64, # for Bayesian Optimization: number of hpcs to search
  nthreads = 16 # implicit learner parallelization (RSFs, xgboost)
)

# Resampling for tuning the train set
resampling = rsmp('repeated_cv', repeats = config$repeats, folds = config$nfolds)

# Measures of performance
harrell_c = msr('surv.cindex')
harrell_c$label = 'HarrellC'
uno_c = msr('surv.cindex', weight_meth = 'G2')
uno_c$label = 'UnoC'
ibrier = msr('surv.graf')
ibrier$label = 'IBrier'
measure = harrell_c # used for training
test_measures = list(harrell_c, uno_c) # used for bootstrap on the test set

config$tune_measure = measure$label
config$test_measures = sapply(test_measures, function(x) x$label)

# learner ids (will be included in the benchmark result)
lrn_ids = c('coxnet', 'rsf_cindex', 'rsf_logrank', 'coxboost', 'xgboost_cox')

# Tasks ----
#' See `prepare_tasks_after_fs.R`
tasks = readRDS(file = 'data/tasks_fs.rds')
task_ids = lapply(tasks, function(task) task$id) %>% unname() %>% unlist()

## ~70%/30% train/test split
data_split = readRDS(file = 'data/data_split.rds')
train_indx = data_split$train_indx # 100
test_indx  = data_split$test_indx # 45

# Learners & Tuning Parameters ----
lrn_tbl = get_lrns_and_ps(nthreads = config$nthreads, ids = lrn_ids)

# Benchmark ----
bench_grid = data.table::CJ(task_id = task_ids, lrn_id = lrn_tbl$id, sorted = FALSE)

## for testing
#bench_grid = bench_grid[c(1,123)]

n_benchmarks = nrow(bench_grid)

res_list = list()
for(row_id in 1:n_benchmarks) {
  message('\nBenchmark ', row_id, '/', n_benchmarks)
  message('#################')

  task_id = bench_grid[row_id]$task_id
  lrn_id  = bench_grid[row_id]$lrn_id

  task = tasks[[task_id]]
  learner = lrn_tbl[id == lrn_id]$learner[[1L]]
  search_space = lrn_tbl[id == lrn_id]$param_set[[1L]]

  # CoxBoost: don't penalize clinical features in case they are included
  # Note: clinical variables should be first in a multimodal dataset
  if (grepl(pattern = '^Clinical', task_id) && (lrn_id == 'coxboost')) {
    learner$param_set$values$unpen.index = 1:length(tasks$Clinical$feature_names)
  }

  # TODO
  #' Set `early_stopping_set = test/train?` for `xgboost_{cox/aft}_reg` learner
  #' when `mlr3extralearners` is updated
  #' distr prediction is needed pending on train measure or more general?

  res = run_at(learner, task, train_indx, test_indx, resampling, measure,
    config$nevals, search_space, all_hpcs_perf = TRUE, boot_test = TRUE,
    test_measures, nrsmps = 1000, nthreads = 1)

  res_list[[row_id]] = tibble::tibble(
    task_id = task_id,
    lrn_id  = lrn_id,
    results = list(res)
  )
}

bench_res = dplyr::bind_rows(res_list)

# Save results ----
cur_time = format(Sys.time(), "%d%m%Y_%H%M%S")
saveRDS(list(config = config, bench_res = bench_res),
  file = paste0(res_path, 'bench_res_', cur_time, '.rds'))
