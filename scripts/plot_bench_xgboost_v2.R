# plot results from `bench_xgboost_v2.R` script,
# so run that script first to generate the output files
library(mlr3verse)
library(tidyverse)

if (!dir.exists('img/xgboost/')) {
  dir.create(res_path)
}

# XGBoost Cox (survival-cox, 7 HPs) ----
res = readRDS(file = 'results/xgboost/xgboost_cox.rds')
#res = readRDS(file = 'results/xgboost/xgboost_cox_random_search.rds')

## Tuning results ----
as.data.table(res$xgboost_at_cox$tuning_instance$result_learner_param_vals)
as.data.table(res$xgboost_cox$param_set$values) # learner with best hpc
as.data.table(res$xgboost_at_cox$archive$best()$x_domain[[1L]])

## Performance plot ----
# that includes test, train + train CV C-index
set1_colors = RColorBrewer::brewer.pal(n = 3, name = 'Set1')
res$hpc_res %>% ggplot(aes(x = index)) +
  geom_line(aes(y = train_cindex, color = 'Train'), size = 0.3) +
  geom_line(aes(y = traincv_cindex, color = 'Train (CV)'), size = 0.3) +
  geom_line(aes(y = test_cindex, color = 'Test'), size = 0.7) +
  scale_color_manual(name = '',
    values = c(
      'Train' = set1_colors[1],
      'Train (CV)' = set1_colors[2],
      'Test' = set1_colors[3])
    ) +
  geom_hline(yintercept = 0.5, linetype = 'dotted', color = 'red') +
  labs(x = 'Number of evaluations/hpcs (BO)', y = 'C-index') +
  theme_bw(base_size = 14) + theme(legend.position = 'top')
ggsave(filename = 'img/xgboost/perf_cox.png', width = 7, height = 5, dpi = 450)

## Confidence Intervals ----
plot(res$boot_cox$boot_res)
res$boot_cox$bootci_res$bca

# verify C-index on test and train sets
tasks = readRDS(file = 'data/tasks.rds')
task_mRNA = tasks$mRNA
res$xgboost_cox$predict(task_mRNA, row_ids = res$test_indx)$score()
res$boot_cox$boot_res$t0 # same
res$xgboost_cox$predict(task_mRNA, row_ids = res$train_indx)$score()
# learner uses the hpc that resulted in the best average CV C-index
res$hpc_res %>% arrange(desc(traincv_cindex)) %>% slice(1)

# XGBoost AFT (survival-aft, 7 HPs) ----
res = readRDS(file = 'results/xgboost/xgboost_aft.rds')

## Performance plot ----
# that includes test, train + train CV C-index
res$hpc_res %>% ggplot(aes(x = index)) +
  geom_line(aes(y = train_cindex, color = 'Train'), size = 0.3) +
  geom_line(aes(y = traincv_cindex, color = 'Train (CV)'), size = 0.3) +
  geom_line(aes(y = test_cindex, color = 'Test'), size = 0.7) +
  scale_color_manual(name = '',
    values = c(
      'Train' = set1_colors[1],
      'Train (CV)' = set1_colors[2],
      'Test' = set1_colors[3])
  ) +
  geom_hline(yintercept = 0.5, linetype = 'dotted', color = 'red') +
  labs(x = 'Number of evaluations/hpcs (BO)', y = 'C-index') +
  theme_bw(base_size = 14) + theme(legend.position = 'top')
ggsave(filename = 'img/xgboost/perf_aft.png', width = 7, height = 5, dpi = 450)

## Confidence Intervals ----
plot(res$boot_aft$boot_res)
res$boot_aft$bootci_res$bca
