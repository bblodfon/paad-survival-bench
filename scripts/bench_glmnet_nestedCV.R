library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
suppressMessages(library(dplyr))
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

glmnet_lrn = lrn('surv.glmnet', standardize = FALSE, id = 'CoxNet',
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

# Less logging
lgr::get_logger("bbotk")$set_threshold("warn")
lgr::get_logger("mlr3")$set_threshold("warn")

# Runs the outer loop in parallel and the inner loop sequentially
future::plan(list("multisession", "sequential"))
#future::plan("multisession")
set.seed(42)

tic()
# `store_model = TRUE` to get the inner tuning archives!
bm_res = benchmark(design, store_models = TRUE, store_backends = FALSE)
toc()

# save performance results
perf_res = bm_res$score() %>%
  as_tibble() %>%
  select(task_id, learner_id, surv.cindex)
perf_res$learner_id = 'CoxNet'
saveRDS(perf_res, file = 'results/perf_res_glmnet_nestedCV.rds')

# save inner tuning archives
hpo_res = extract_inner_tuning_archives(bm_res) %>%
  as_tibble() %>%
  select(-c(resample_result, resampling_id))
saveRDS(object = hpo_res, file = 'results/hpo_res_glmnet_nestedCV.rds')
