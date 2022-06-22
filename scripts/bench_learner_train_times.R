# Getting an estimate of the train time for
# some (untuned) learners on the mRNA dataset
library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)

# Task ----
# first run `scripts/prepare_tasks.R`
tasks = readRDS(file = 'data/tasks.rds')

mRNA_task  = tasks$mRNA

# Learners ----
glmnet_lrn  = lrn('surv.glmnet', id = 'CoxLasso', standardize = FALSE,
  lambda = 0.001, alpha = 1, maxit = 10^3)
rpart_lrn   = lrn('surv.rpart', id = 'Survival Tree')
ranger_lrn  = lrn('surv.ranger', id = 'Survival Forest',
  verbose = TRUE, num.trees = 1000)
ranger_lrn_parallel = lrn('surv.ranger', id = 'Survival Forest',
  verbose = FALSE, num.trees = 1000, num.threads = 4)
xgboost_lrn = lrn('surv.xgboost', id = 'XGBoost Survival Learner') # nrounds = 1
xgboost_lrn_500 = lrn('surv.xgboost', id = 'XGBoost Survival Learner',
  nrounds = 500)

# Helper functions ----
get_train_times = function(learner_list) {
  sapply(learner_list, function(learner) {learner$timings['train']}) %>% unname()
}

print_times = function(times_vec, prefix_msg = '') {
  print(paste0(prefix_msg, paste0(round(times_vec, 3), collapse = ', ')))
}

print_mean_time = function(times_vec, prefix_msg = '') {
  print(paste0(prefix_msg, paste0(round(mean(times_vec), 3))))
}

# CV resampling
set.seed(42)
cv_rsmp = rsmp('cv', folds = 5)
cv_rsmp$instantiate(mRNA_task)

# Test train times ----
#if (FALSE) {
rr_glmnet = resample(task = mRNA_task, learner = glmnet_lrn, resampling = cv_rsmp)
glmnet_times = get_train_times(rr_glmnet$learners)
print_times(glmnet_times, 'Glmnet train times: ')
print_mean_time(glmnet_times, 'Glmnet mean train time: ') # 1-2 secs

rr_rpart = resample(task = mRNA_task, learner = rpart_lrn, resampling = cv_rsmp)
rpart_times = get_train_times(rr_rpart$learners)
print_times(rpart_times, 'Rpart train times: ')
print_mean_time(rpart_times, 'Rpart mean train time: ') # 8-12 secs

rr_ranger = resample(task = mRNA_task, learner = ranger_lrn, resampling = cv_rsmp)
ranger_times = get_train_times(rr_ranger$learners)
print_times(ranger_times, 'Ranger train times: ')
print_mean_time(ranger_times, 'Ranger mean train time: ') # 45-60 secs (1000 trees)

rr_ranger_parallel = resample(mRNA_task, ranger_lrn_parallel, cv_rsmp)
ranger_parallel_times = get_train_times(rr_ranger_parallel$learners)
print_times(ranger_parallel_times, 'Ranger train times (multi-threaded): ')
print_mean_time(ranger_parallel_times, 'Ranger mean train time (multi-threaded): ') # 17 secs

rr_xgboost = resample(task = mRNA_task, learner = xgboost_lrn, resampling = cv_rsmp)
xgboost_times = get_train_times(rr_xgboost$learners)
print_times(xgboost_times, 'XGBoost train times: ')
print_mean_time(xgboost_times, 'XGBoost mean train time: ') # 0.265 secs (nrounds = 1)

rr_xgboost_500 = resample(task = mRNA_task, learner = xgboost_lrn_500, resampling = cv_rsmp)
xgboost500_times = get_train_times(rr_xgboost_500$learners)
print_times(xgboost500_times, 'XGBoost train times: ')
print_mean_time(xgboost500_times, 'XGBoost mean train time: ') # 48 secs (nrounds = 500)

xgboost_lrn_500$param_set$values$nthread = 4 # parallelize
rr_xgboost_500 = resample(task = mRNA_task, learner = xgboost_lrn_500, resampling = cv_rsmp)
xgboost500_times = get_train_times(rr_xgboost_500$learners)
print_times(xgboost500_times, 'XGBoost train times (multi-threaded): ')
print_mean_time(xgboost500_times, 'XGBoost mean train time (multi-threaded): ') # 27 secs (nrounds = 500, nthreads = 4)
#}

# RSFs splitrule test ----
# RFs with different split rules
ranger_cindex = lrn('surv.ranger', id = 'Survival Forest', verbose = TRUE,
  num.trees = 1000, num.threads = 10, splitrule = 'C') # C-index
ranger_logrank = lrn('surv.ranger', id = 'Survival Forest', verbose = TRUE,
  num.trees = 1000, num.threads = 10, splitrule = 'logrank')
ranger_maxstat = lrn('surv.ranger', id = 'Survival Forest', verbose = TRUE,
  num.trees = 1000, num.threads = 10, splitrule = 'maxstat')
ranger_extratrees = lrn('surv.ranger', id = 'Survival Forest', verbose = TRUE,
  num.trees = 1000, num.threads = 10, splitrule = 'extratrees')

rranger_cindex  = resample(mRNA_task, ranger_cindex, cv_rsmp)
rranger_logrank = resample(mRNA_task, ranger_logrank, cv_rsmp)
rranger_maxstat = resample(mRNA_task, ranger_maxstat, cv_rsmp)
rranger_extratrees = resample(mRNA_task, ranger_extratrees, cv_rsmp)

rr_ranger = list(logrank = rranger_logrank, cindex = rranger_cindex,
  maxstat = rranger_maxstat, extratrees = rranger_extratrees)

saveRDS(rr_ranger, file = 'results/rr_ranger.rds')
#rr_ranger = readRDS(file = 'results/rr_ranger.rds')

rranger_cindex_times  = get_train_times(rr_ranger$cindex$learners)
rranger_logrank_times = get_train_times(rr_ranger$logrank$learners)
rranger_maxstat_times = get_train_times(rr_ranger$maxstat$learners)
rranger_extratrees_times = get_train_times(rr_ranger$extratrees$learners)

print_times(rranger_cindex_times, 'Ranger train times (C-index): ')
print_mean_time(rranger_cindex_times, 'Ranger mean train time (C-index): ')
print_times(rranger_logrank_times, 'Ranger train times (logrank): ')
print_mean_time(rranger_logrank_times, 'Ranger mean train time (logrank): ')
print_times(rranger_maxstat_times, 'Ranger train times (maxstat): ')
print_mean_time(rranger_maxstat_times, 'Ranger mean train time (maxstat): ')
print_times(rranger_extratrees_times, 'Ranger train times (extratrees): ')
print_mean_time(rranger_extratrees_times, 'Ranger mean train time (extratrees): ')
