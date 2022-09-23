#'##################################################
#' Benchmark
#' - CoxPH using clinical features (baseline)
#' Data resamplings (x1000):
#' 1. Simple train/test split + bootstrap C-index
#' on test set (based on mRNA data split benchmarks)
#' 2. Bootstrap whole dataset
#' 3. Random permutation of time + status
#'   3.1 Shuffle all targets before resampling
#'   3.2 Shuffle only each train set's targets
#'##################################################
library(mlr3verse)
library(mlr3proba)
library(progressr)
library(tictoc)
library(tidyverse)
source('scripts/boot_mlr3.R')
source('scripts/pipeops.R')

# Reproducibility
set.seed(42)

# Progress bars
options(progressr.enable = TRUE)
handlers(global = TRUE)
handlers('progress')

# Logging (less)
lgr::get_logger('bbotk')$set_threshold('warn')
lgr::get_logger('mlr3')$set_threshold('warn')

# Clinical Task ----
# 6 features
task_clinical = readRDS(file = 'data/task_clinical.rds')

# CoxPH Learner ----
## Bootstrap test set ----
cox = lrn('surv.coxph')

## ~70%/30% train/test split (get indexes from previous benchmarks)
res = readRDS(file = 'results/xgboost/xgboost_cox.rds')
train_indx = res$train_indx
test_indx  = res$test_indx

## get model on train set
cox$train(task_clinical, train_indx)
summary(cox$model)

print('Boostrap test set (x1000)')
boot_test = boot_ci(data = task_clinical$data()[test_indx], nrsmps = 1000,
  learner = cox, num_threads = 4, include_boot_res = TRUE)

## Bootstrap resampling ----
future::plan('multisession')
print('Boostrap all dataset (x1000)')
tic()
rs = resample(task = task_clinical, learner = lrn('surv.coxph'),
  resampling = rsmp('bootstrap', repeats = 1000, ratio = 1))
toc()

boot_rs_harrel_cindex = rs$score(measures = msr('surv.cindex'))$surv.cindex
boot_rs_uno_cindex = rs$score(measures = msr('surv.cindex', weight_meth = 'G2'))$surv.cindex

## Random task benchmark ----
## Permute time + status in clinical_task (x1000)

#' @param x a {data.table} with columns `time` and `status`, representing
#' the survival target of a task
perm_trg = function(x) {
  x$time = sample(x$time)
  x$status = sample(x$status)
  x
}

features_dt = task_clinical$data(cols = task_clinical$feature_names)
targets_dt  = task_clinical$data(cols = task_clinical$target_names)

task_list = list()
num_tasks = 1000
for (i in 1:num_tasks) {
  task_list[[i]] = mlr3proba::as_task_surv(
    x = cbind(perm_trg(targets_dt), features_dt),
    time = 'time', event = 'status', id = 'Clinical'
  )
}

# check that some random tasks' targets are different
for (i in sample(num_tasks, 5)) {
  print(task_list[[i]]$head(4)[,c(1,2)])
}

# both train and test set targets will be shuffled now (which might not make sense)
bm_grid = benchmark_grid(task_list, lrn('surv.coxph'), rsmp('holdout', ratio = 0.7))

future::plan('multisession')
print('Start: Cox random task benchmark')
tic()
cox_bm = benchmark(bm_grid, store_models = TRUE, store_backends = TRUE)
toc()
print('Finish: Cox random task benchmark')

## save C-index results
rnd_har_cindex = cox_bm$score(msr('surv.cindex'))[['surv.cindex']]
rnd_uno_cindex = cox_bm$score(msr('surv.cindex', weight_meth = 'G2'))[['surv.cindex']]

## Random task benchmark (2) ----
## Permute time + status in clinical_task (x1000)
## Now, test set targets should not be shuffled, see
## 'scripts/pipeops.R' where `poss` PipeOp is implemented
learner = as_learner(poss %>>% po('learner', lrn('surv.coxph')))
cop = po('copy', 1000)
task_list = cop$train(list(task_clinical)) # list of 1000 identical `clinical_task`s
bm_grid2 = benchmark_grid(task_list, learner, rsmp('holdout', ratio = 0.7))

future::plan('multisession')
print('Start: Cox random task benchmark 2')
tic()
cox_bm2 = benchmark(bm_grid2, store_models = TRUE, store_backends = TRUE)
toc()
print('Finish: Cox random task benchmark 2')

rnd2_har_cindex = cox_bm2$score(msr('surv.cindex'))[['surv.cindex']]
rnd2_uno_cindex = cox_bm2$score(msr('surv.cindex', weight_meth = 'G2'))[['surv.cindex']]

# store all CoxPH results ----
cox_res = tibble(boot_rs_harrel_cindex = boot_rs_harrel_cindex,
  boot_rs_uno_cindex = boot_rs_uno_cindex,
  rnd_har_cindex = rnd_har_cindex, rnd_uno_cindex = rnd_uno_cindex,
  rnd2_har_cindex = rnd2_har_cindex, rnd2_uno_cindex = rnd2_uno_cindex)

saveRDS(object = list(boot_test = boot_test, cox_res = cox_res),
  file = 'results/coxph.rds')

# Boxplot C-index Comparison ----
res = readRDS(file = 'results/coxph.rds')

res$cox_res %>%
  add_column(boot_test_har_cindex = res$boot_test$boot_res$t[,1], .before = 1) %>%
  tidyr::pivot_longer(cols = everything(), names_to = 'method', values_to = 'cindex') %>%
  mutate(method = case_when(
    method == 'boot_test_har_cindex' ~ 'BootTest-HarC',
    method == 'boot_rs_harrel_cindex' ~ 'BootRS-HarC',
    method == 'boot_rs_uno_cindex' ~ 'BootRS-UnoC',
    method == 'rnd_har_cindex' ~ 'ShuffAll-HarC',
    method == 'rnd_uno_cindex' ~ 'ShuffAll-UnoC',
    method == 'rnd2_har_cindex' ~ 'ShuffTrain-HarC',
    method == 'rnd2_uno_cindex' ~ 'ShuffTrain-UnoC'
  )) %>%
  ggplot(aes(x = method, y = cindex, color = method)) +
  geom_boxplot(show.legend = FALSE) +
  geom_hline(yintercept = 0.5, linetype = 'dotted', color = 'red') +
  labs(x = 'Resampling Method - Metric', y = 'C-index') +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(filename = 'img/coxph.png', width = 4, height = 5, dpi = 350)
