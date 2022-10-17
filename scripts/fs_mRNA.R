#'###################################################
#' Apply an ensemble Feature Selection (FS) strategy
#' on the mRNA dataset, using multiple learners and
#' RFE and GA wrapper-based FS methods
#'###################################################
library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
library(progressr)
library(tictoc)
source('scripts/helpers.R')

res_path = 'results/fs/mRNA'
if (!dir.exists(res_path)) {
  dir.create(res_path, recursive = TRUE)
}

# Reproducibility
set.seed(42)

# Progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers('progress')

# Less logging
lgr::get_logger('bbotk')$set_threshold('warn')
lgr::get_logger('mlr3')$set_threshold('warn')

# Global variables ----
config = list(
  num_threads = 30, # implicit parallelization for RSFs
  num_trees = 250,
  mtry_ratio = 0.1, # try 10% of the features in the RSF splits
  min_node_size = 3, # default for survival (RSF)
  repeats = 100, # how many times to run a wrapper method using a specific learner
  wrapper_methods = c('rfe', 'ga'),
  rfe_n_features = 2, # number of features to stop RFE
  rfe_feature_fraction = 0.85, # %features to keep in each iteration of RFE
  ga_iters = 100, # GA iterations
  ga_zeroToOneRatio = 200, # sparse features subsets in GA
  ga_popSize = 1000 # initial population of feature subsets
)

# mRNA Task ----
task_mRNA = readRDS(file = 'data/tasks.rds')$mRNA
mRNA_indx = readRDS(file = 'data/mRNA_indx.rds') # same train/test split as in other benchmarks
task = task_mRNA$clone()$filter(mRNA_indx$train_indx) # mRNA task to use for FS

# Learners ----
## RSFs ----
rsf_cindex = lrn('surv.ranger', verbose = FALSE, id = 'rsf_cindex',
  num.trees = config$num_trees, mtry.ratio = config$mtry_ratio,
  min.node.size = config$min_node_size,
  num.threads = config$num_threads, importance = 'permutation',
  splitrule = 'C') # C-index

rsf_logrank = lrn('surv.ranger', verbose = FALSE, id = 'rsf_logrank',
  num.trees = config$num_trees, mtry.ratio = config$mtry_ratio,
  min.node.size = config$min_node_size,
  num.threads = config$num_threads, importance = 'permutation',
  splitrule = 'logrank')

rsf_maxstat = lrn('surv.ranger', verbose = FALSE, id = 'rsf_maxstat',
  num.trees = config$num_trees, mtry.ratio = config$mtry_ratio,
  min.node.size = config$min_node_size,
  num.threads = config$num_threads, importance = 'permutation',
  splitrule = 'maxstat')

## CoxLasso ----
# with CV tuning of lambda
coxlasso = lrn('surv.cv_glmnet', alpha = 1, s = 'lambda.min', # more coefficients
  id = 'cox_lasso', nfolds = 5, type.measure = 'C', standardize = FALSE,
  maxit = 10^3, fallback = lrn('surv.kaplan'))

## CoxBoost ----
coxboost = lrn('surv.coxboost', id = 'cox_boost',
  standardize = FALSE, # data already standardized
  return.score = FALSE, # don't need this in the output
  fallback = lrn('surv.kaplan'), stepno = 100, criterion = 'hpscore')

learners = list(rsf_cindex, rsf_logrank, rsf_maxstat, coxlasso, coxboost)

# Ensemble FS ----
res_list = list()
index = 1
for (method in config$wrapper_methods) {
  for (learner in learners) {
    message('### ', method, ' - ', learner$id, ' ###')

    tic()
    res = run_wrapper_fs(learner = learner, task = task, method = method,
      repeats = config$repeats,
      rfe_n_features = config$rfe_n_features,
      rfe_feature_fraction = config$rfe_feature_fraction,
      ga_iters = config$ga_iters,
      ga_zeroToOneRatio = config$ga_zeroToOneRatio,
      ga_popSize = config$ga_popSize)
    toc()

    if (!is.null(res)) {
      res_list[[index]] = res
      index = index + 1

      # save intermediate output
      saveRDS(res, file = paste0(res_path, '/', method, '_', learner$id, '.rds'))
    }
  }
}

total_res = dplyr::bind_rows(res_list)

# save all results + configuration
saveRDS(object = list(res = total_res, config = config),
        file   = paste0(res_path, '/fs_results.rds'))

if (FALSE) {
  # Some test runs
  rsf_maxstat = lrn('surv.ranger', verbose = FALSE, id = 'rsf_maxstat',
    num.trees = 50,
    mtry.ratio = 0.01,
    min.node.size = 3,
    splitrule = 'maxstat',
    num.threads = num_threads,
    importance = 'permutation')

  ## GA-RSF-maxstat
  res1 = run_wrapper_fs(learner = rsf_maxstat, task = task,
    method = 'ga', repeats = 2, ga_iters = 5)
  ## GA-CoxLasso
  res2 = run_wrapper_fs(learner = coxlasso, task = task,
    method = 'ga', repeats = 2, ga_iters = 5)
  res22 = run_wrapper_fs(learner = coxboost, task = task,
    method = 'ga', repeats = 2, ga_iters = 5)

  ## RFE-CoxLasso
  res3 = run_wrapper_fs(learner = coxlasso, task = task,
    method = 'rfe', repeats = 2, rfe_n_features = 2, rfe_feature_fraction = 0.05)
  ## RFE-RSF-maxstat
  res4 = run_wrapper_fs(learner = rsf_maxstat, task = task,
    method = 'rfe', repeats = 2, rfe_n_features = 2, rfe_feature_fraction = 0.05)
  ## RFE-CoxBoost
  res5 = run_wrapper_fs(learner = coxboost, task = task,
    method = 'rfe', repeats = 2, rfe_n_features = 2, rfe_feature_fraction = 0.05)

  # after FS has been done
  train_task = task_mRNA$clone()$select(res1$sel_features[[2]])$filter(mRNA_indx$train_indx)
  test_task  = task_mRNA$clone()$select(res1$sel_features[[2]])$filter(mRNA_indx$test_indx)

  ## CoxLasso
  coxlasso$reset()
  coxlasso$train(train_task)
  coxlasso$predict(test_task)$score()

  ## RSF
  rsf_cindex$reset()
  rsf_cindex$train(train_task)
  rsf_cindex$model

  rsf_cindex$predict(test_task)$score()
  rsf_cindex$predict(test_task)$score(uno_cindex, task = train_task,
    train_set = 1:train_task$nrow)
}
