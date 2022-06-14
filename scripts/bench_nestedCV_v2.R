#'############################################################
# Benchmark tuned CoxNet, Survival Tree and Survival Forests #
# on 3 tasks (2 single, 1 combined) using nested CV          #
#'############################################################
library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
suppressMessages(library(dplyr))
library(tictoc)

# Tasks ----
# first run `scripts/prepare_tasks.R`
tasks = readRDS(file = 'data/tasks.rds')
tasks = list(tasks$mRNA, tasks$CNA, tasks$`mRNA-CNA`)

# Global variables ----
baseline_estimator = 'nelson' # Nelson-Aalen estimator
distr_form = 'aft' # Accelerated failure time
rand_search_tuner = tnr('random_search')
inner_folds = 5
inner_resampling = rsmp('cv', folds = inner_folds)
outer_folds = 4
outer_resampling = rsmp('cv', folds = outer_folds) # outer folds (3-10)
harrell_cindex = msr('surv.cindex')
n_evals = 100 # number of hyperparameter configurations searched
eval_trm = trm("evals", n_evals = n_evals)

# Glmnet learner ----
glmnet_lrn = lrn('surv.glmnet', standardize = FALSE,
  lambda = to_tune(p_dbl(1e-03, 1, logscale = TRUE)),
  alpha = to_tune(0, 1)) # from Ridge to Lasso penalty

glmnet_graph_lrn = ppl('distrcompositor',
  glmnet_lrn,
  estimator = baseline_estimator,
  form = distr_form,
  overwrite = TRUE,
  graph_learner = FALSE
) %>% GraphLearner$new(id = 'CoxNet')

glmnet_at = AutoTuner$new(
  learner = glmnet_graph_lrn,
  resampling = inner_resampling,
  measure = harrell_cindex,
  terminator = eval_trm,
  tuner = rand_search_tuner)

# Survival Tree learner ----
rpart_lrn = lrn('surv.rpart',
  minsplit = to_tune(1, 20),
  cp = to_tune(1e-04, 0.1, logscale = TRUE))

rpart_graph_lrn = ppl('distrcompositor',
  rpart_lrn,
  estimator = baseline_estimator,
  form = distr_form,
  overwrite = TRUE,
  graph_learner = TRUE,
  id = 'SurvivalTree'
)

rpart_at = AutoTuner$new(
  learner = rpart_graph_lrn,
  resampling = inner_resampling,
  measure = harrell_cindex,
  terminator = eval_trm,
  tuner = rand_search_tuner)

# Survival Forest learner ----
ranger_lrn = lrn('surv.ranger', verbose = FALSE, id = 'SurvivalForest',
  num.trees = to_tune(10, 2000),
  mtry.ratio = to_tune(0, 1),
  min.node.size = to_tune(1, 20),
  splitrule = to_tune(c('logrank', 'extratrees', 'C', 'maxstat')))
## no tuning for dependent parameters for `extratrees` (num.random.splits = [1,100])
## and 'maxstat' (alpha, minprop)

ranger_at = AutoTuner$new(
  learner = ranger_lrn,
  resampling = inner_resampling,
  measure = harrell_cindex,
  terminator = eval_trm,
  tuner = rand_search_tuner)

# XGBoost learner ----
# xgboost_lrn = lrn('surv.xgboost')
# id = 'XGBoostSurvival'

## consult:
#mlr3tuningspaces::lts('classif.xgboost.default')
#mlr3tuningspaces::lts('classif.xgboost.rbv2')

# Benchmark (Nested-CV) ----
learners = list(glmnet_at, rpart_at, ranger_at)
design = benchmark_grid(tasks, learners, outer_resampling)

# Less logging
lgr::get_logger("bbotk")$set_threshold("warn")
lgr::get_logger("mlr3")$set_threshold("warn")

# Runs the outer loop in parallel and the inner loop sequentially
future::plan(list("multisession", "sequential"))
set.seed(42)

tic()
# `store_backends = TRUE` for calculating 'surv.graf' and 'surv.intlogloss' scores
bm_res = benchmark(design, store_models = FALSE, store_backends = TRUE)
toc()

# save performance results
graf_score = msr('surv.graf')
intlogloss = msr('surv.intlogloss')
measures = list(harrell_cindex, graf_score, intlogloss)

perf_res = bm_res$score(measures) %>%
  as_tibble() %>%
  select(task_id, learner_id, surv.cindex, surv.graf, surv.intlogloss)
saveRDS(perf_res, file = 'results/perf_res_nestedCV_v2.rds')
