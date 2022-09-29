# plot results from `coxboost_bench.R` script,
# so run that script first to generate the output files
library(mlr3verse)
library(tidyverse)

if (!dir.exists('img/coxboost/')) {
  dir.create('img/coxboost/')
}

# Tasks ----
## mRNA (~10000 features) ----
tasks = readRDS(file = 'data/tasks.rds')
task_mRNA = tasks$mRNA

## Clinical (8 features) ----
task_clinical = readRDS(file = 'data/task_clinical.rds')

## Clinical + mRNA ----
# make combined task with clinical first and then mRNA features
cmb_task = task_clinical$clone(deep = TRUE)
cmb_task$cbind(task_mRNA$data(cols = task_mRNA$feature_names))

# Read result data ----
## CoxPH (clinical)
res_coxph = readRDS(file = 'results/coxph.rds')
## CoxBoost (mRNA data only)
res_mrna = readRDS(file = 'results/coxboost/coxboost_mRNA_test.rds')
## CoxBoost (clinical + mRNA)
res_cmb = readRDS(file = 'results/coxboost/coxboost_clinical_mRNA_test.rds')

## verify train and test indexes are the same
all(res_coxph$train_indx == res_mrna$train_indx)
all(res_mrna$train_indx == res_cmb$train_indx)
all(res_coxph$test_indx == res_mrna$test_indx)
all(res_mrna$test_indx == res_cmb$test_indx)
train_indx = res_cmb$train_indx
test_indx  = res_cmb$test_indx

# CoxBoost - perf plot (mRNA only) ----
## test, train + train CV C-index for all hps
set1_colors = RColorBrewer::brewer.pal(n = 3, name = 'Set1')
res_mrna$hpc_res %>% ggplot(aes(x = index)) +
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
ggsave(filename = 'img/coxboost/perf_mRNA.png', width = 7, height = 5, dpi = 450)

# CoxBoost - perf plot (clinical + mRNA) ----
## test, train + train CV C-index for all hpcs
res_cmb$hpc_res %>% ggplot(aes(x = index)) +
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
ggsave(filename = 'img/coxboost/perf_cmb.png', width = 7, height = 5, dpi = 450)

# CI Plot ----
## Compare CoxPH vs CoxBoost (mRNA only) vs CoxBoost (clinical + mRNA)
## Boostrap test CIs
cox_ci = unname(res_coxph$boot_test$bootci_res$bca)
coxboost_mrna_ci = unname(res_mrna$boot_cox$bootci_res$bca)
coxboost_cmb_ci  = unname(res_cmb$boot_cox$bootci_res$bca)

res_ci = tibble::tibble(name = c('CoxPH-Clinical', 'CoxBoost-mRNA', 'CoxBoost-Clinical-mRNA'),
  # C-index on train set
  cindex_train = c(
    res_coxph$cox_lrn$predict(task_clinical, train_indx)$score(),
    res_mrna$coxboost_lrn$predict(task_mRNA, train_indx)$score(),
    res_cmb$coxboost_lrn$predict(cmb_task, train_indx)$score()
  ),
  # C-index on test set
  cindex_test = c(
    res_coxph$boot_test$boot_res$t0, # res_coxph$cox_lrn$predict(task_clinical, test_indx)$score() # same no need to recalculate
    res_mrna$boot_cox$boot_res$t0, # res_mrna$coxboost_lrn$predict(task_mRNA, train_indx)$score() # same no need to recalculate
    res_cmb$boot_cox$boot_res$t0),
  # C-index intervals on test set
  low_ci = c(cox_ci[,4], coxboost_mrna_ci[,4], coxboost_cmb_ci[,4]),
  high_ci = c(cox_ci[,5], coxboost_mrna_ci[,5], coxboost_cmb_ci[,5]))

res_ci %>%
  mutate(name = factor(name, levels = c('CoxPH-Clinical', 'CoxBoost-mRNA', 'CoxBoost-Clinical-mRNA'))) %>%
  ggplot(aes(x = name, y = cindex_test, color = name)) +
  geom_pointrange(aes(ymin = low_ci, ymax = high_ci)) +
  geom_errorbar(aes(ymin = low_ci, ymax = high_ci), width = 0.1) +
  geom_point(aes(x = name, y = cindex_train), size = 3) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'red') +
  theme_bw(base_size = 16) + ylim(0.2, 1) +
  labs(x = '', y = 'C-index', title = 'Bootstrap 95% CIs on test set +', subtitle = 'C-indexes on train set') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(legend.position = 'none', title = element_text(size = 11))
ggsave(filename = 'img/coxboost/ci_plot.png', width = 4, height = 5, dpi = 450)
