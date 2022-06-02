library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
suppressMessages(library(dplyr))
library(readr)
library(tictoc)

# Tasks ----
# first run `scripts/prepare_tasks.R`
tasks = readRDS(file = 'data/tasks.rds')
tasks = list(tasks$mRNA, tasks$miRNA, tasks$Methyl, tasks$CNA)

# Set variables for AutoTuner ----
rand_search_tuner = tnr('random_search')
inner_folds = 5
inner_resampling = rsmp('cv', folds = inner_folds)
outer_folds = 4
outer_resampling = rsmp('cv', folds = outer_folds) # outer folds (3-10)
harrell_cindex = msr('surv.cindex')
n_evals = 100 # number of hyperparameter configurations searched
eval_trm = trm("evals", n_evals = n_evals)

# glmnet learner ----
# mlr3tuningspaces::lts("classif.glmnet.rbv2") # 1000 too much
min_lambda = 0.001
max_lambda = 10

glmnet_lrn = lrn('surv.glmnet', standardize = FALSE,
  lambda = to_tune(p_dbl(min_lambda, max_lambda, logscale = TRUE)),
  alpha = to_tune(0, 1))

glmnet_at = AutoTuner$new(
  learner = glmnet_lrn,
  resampling = inner_resampling,
  measure = harrell_cindex,
  terminator = eval_trm,
  tuner = rand_search_tuner)

# Benchmark (Nested-CV) ----
design = benchmark_grid(tasks = tasks, glmnet_at, outer_resampling)

# Runs the outer loop in parallel and the inner loop sequentially
future::plan(list("multisession", "sequential"))
#future::plan("multisession")
set.seed(42)

tic()
bm_res = benchmark(design, store_models = TRUE, store_backends = FALSE)
toc()

saveRDS(bm_res, file = 'results/bm_res_glmnet_nested_cv.rds')
