# plot results from `bench_xgboost2.R` script,
# so run that script first to generate the output files
library(mlr3verse)
library(tidyverse)

# XGBoost (survival-cox, 7 HPs) ----
res = readRDS(file = 'results/xgboost/xgboost_cox_veteran.rds')

## Tuning results
as.data.table(res$xgboost_at_cox$tuning_instance$result_learner_param_vals)
as.data.table(res$xgboost_cox$param_set$values) # learner with best hpc
as.data.table(res$xgboost_at_cox$archive$best()$x_domain[[1L]])

## Plots
### Performance plot
autoplot(res$xgboost_at_cox$tuning_instance, type = 'performance')

# include test, train + train CV C-index
set1_colors = RColorBrewer::brewer.pal(n = 3, name = 'Set1')
res$hpc_res %>% ggplot(aes(x = index)) +
  geom_line(aes(y = train_cindex, color = 'Train')) +
  geom_line(aes(y = traincv_cindex, color = 'Train (CV)')) +
  geom_line(aes(y = test_cindex, color = 'Test')) +
  scale_color_manual(name = '',
    values = c('Train' = set1_colors[1], 'Train (CV)' = set1_colors[2],
      'Test' = set1_colors[3])) +
  geom_hline(yintercept = 0.5, linetype = 'dotted', color = 'red') +
  labs(x = 'Number of evaluations/hpcs (BO)', y = 'C-index') +
  theme_bw(base_size = 14) + theme(legend.position = 'top')

res$hpc_res %>% tidyr::pivot_longer(cols = c(everything(),-index),
  values_to = c('cindex')) %>%
  mutate(name = case_when(
    name == 'train_cindex' ~ 'Train',
    name == 'test_cindex' ~ 'Test',
    name == 'traincv_cindex' ~ 'Train (CV)',
  )) %>%
  ggplot(aes(x = index, y = cindex, color = name)) +
  scale_color_brewer(palette = 'Set1') +
  geom_line(aes(linetype = name)) +
  geom_hline(yintercept = 0.5, linetype = 'dotted', color = 'red') +
  labs(x = 'Number of evaluations/hpcs (BO)', y = 'C-index') +
  theme_bw(base_size = 14)

### Parameter plot
autoplot(res$xgboost_at_cox$tuning_instance, type = 'parameter', trafo = TRUE)

### Example Surface plot
autoplot(res$xgboost_at_cox$tuning_instance, type = 'surface',
  cols_x = c('nrounds', 'x_domain_eta'))

# Confidence Interval plot ----
## QQ-plots seems quite normal to me :)
plot(xgboost1_res$boot_ci1$boot_res)

xgboost1_ci = xgboost1_res$boot_ci1$bootci_res$bca

# Calculate C-index on train set
test_indx = xgboost1_res$test_indx
tasks = readRDS(file = 'data/tasks.rds')
task_mRNA = tasks$mRNA
train_indx = setdiff(seq_len(task_mRNA$nrow), test_indx)

xgboost_ci = tibble::tibble(name =
  c('XGBoost-Cox-3HPs', 'XGBoost-Cox-6HPs', 'XGBoost-AFT-5HPs', 'XGBoost-AFT-8HPs'),
  # Average CV C-index for the hp config that BO choose (mean surrogate, not best C-index)
  cindex_cv_train = c(
    xgboost1_res$xgboost_at1$tuning_result$surv.cindex,
    xgboost2_res$xgboost_at2$tuning_result$surv.cindex,
    xgboost3_res$xgboost_at3$tuning_result$surv.cindex,
    xgboost4_res$xgboost_at4$tuning_result$surv.cindex),
  # C-index on train set (~1, overfit 100%!)
  cindex_train = c(
    xgboost1_res$xgboost_at1$predict(task_mRNA, train_indx)$score(),
    xgboost2_res$xgboost_at2$predict(task_mRNA, train_indx)$score(),
    xgboost3_res$xgboost_at3$predict(task_mRNA, train_indx)$score(),
    xgboost4_res$xgboost_at4$predict(task_mRNA, train_indx)$score()
  ),
  # C-index on test set
  cindex_test = c(xgboost1_res$boot_ci1$boot_res$t0, xgboost2_res$boot_ci2$boot_res$t0,
    xgboost3_res$boot_ci3$boot_res$t0, xgboost4_res$boot_ci4$boot_res$t0),
  # C-index intervals on test set
  low_ci = c(unname(xgboost1_ci)[,4], unname(xgboost2_ci)[,4],
    unname(xgboost3_ci)[,4], unname(xgboost4_ci)[,4]),
  high_ci = c(unname(xgboost1_ci)[,5], unname(xgboost2_ci)[,5],
    unname(xgboost3_ci)[,5], unname(xgboost4_ci)[,5]))

xgboost_ci %>%
  mutate(name = factor(name, levels = c('XGBoost-Cox-3HPs', 'XGBoost-Cox-6HPs',
    'XGBoost-AFT-5HPs', 'XGBoost-AFT-8HPs'))) %>%
  ggplot(aes(x = name, y = cindex_test, color = name)) +
  geom_pointrange(aes(ymin = low_ci, ymax = high_ci)) +
  geom_errorbar(aes(ymin = low_ci, ymax = high_ci), width = 0.1) +
  geom_point(aes(x = name, y = cindex_cv_train, size = 1)) +
  geom_point(aes(x = name, y = cindex_train, size = 1)) +
  geom_hline(yintercept = 0.5, linetype = 'dotted', color = 'red') +
  theme_bw(base_size = 16) + ylim(0.2, 1) +
  labs(x = '', y = 'C-index', title = 'Bootstrap 95% CIs on test set + C-indexes on train set') +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme(legend.position = 'none', title = element_text(size = 12))
ggsave(filename = 'img/xgboost_bootci.png', width = 7, height = 5, dpi = 450)

# Train times vs eta ----
tbl = tibble(name = c('XGBoost-Cox-3HPs', 'XGBoost-Cox-6HPs',
  'XGBoost-AFT-5HPs', 'XGBoost-AFT-8HPs'),
  train_times = c(xgboost1_res$xgboost_at1$timings['train'],
    xgboost2_res$xgboost_at2$timings['train'],
    xgboost3_res$xgboost_at3$timings['train'],
    xgboost4_res$xgboost_at4$timings['train']),
  sum_eta = c(sum(exp(xgboost1_res$xgboost_at1$archive$data$eta)),
    sum(exp(xgboost2_res$xgboost_at2$archive$data$eta)),
    sum(exp(xgboost3_res$xgboost_at3$archive$data$eta)),
    sum(exp(xgboost4_res$xgboost_at4$archive$data$eta))),
  sum_nrounds = c(sum(xgboost1_res$xgboost_at1$archive$data$nrounds),
    sum(xgboost2_res$xgboost_at2$archive$data$nrounds),
    sum(xgboost3_res$xgboost_at3$archive$data$nrounds),
    sum(xgboost4_res$xgboost_at4$archive$data$nrounds))
)
tbl %>% arrange(train_times) # the larger the train_times, the smaller the eta

# Overfit investigation plot ----
## Plot figure with all hpcs (hyperparameter configurations) on x-axis
## and train C-index + test C-index (2 lines) on y-axis

## XGBoost-Cox (3 hps) ----
hpc_list = xgboost1_res$xgboost_at1$archive$data$x_domain
stopifnot(length(hpc_list) == 250)

base_lrn = xgboost1_res$xgboost_at1$learner$clone(deep = TRUE)
base_lrn$reset()
base_lrn$param_set$values$nthread = 4

res = list()
i = 1
for (hpc in hp_list) {
  print(i)
  base_lrn$param_set$values =
    mlr3misc::insert_named(base_lrn$param_set$values, hpc)
  base_lrn$train(task_mRNA, train_indx)
  res[[i]] = tibble(
    train_cindex = base_lrn$predict(task_mRNA, train_indx)$score(),
    test_cindex  = base_lrn$predict(task_mRNA, test_indx) $score()
  )
  i = i + 1
  base_lrn$reset()
}

