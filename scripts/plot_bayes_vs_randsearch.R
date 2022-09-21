# plot results from `bayesian_vs_randomsearch.R` script,
# so run that script first to generate the output files
library(mlr3verse)
library(dplyr)
library(ggplot2)

# Load results ----
res_path = 'results/bayes_vs_randsearch/'

## CoxNet ----
coxnet_tuning_res = readRDS(file = paste0(res_path, 'coxnet_tuning_res.rds'))
coxnet_bootci = readRDS(file = paste0(res_path, 'coxnet_bootci.rds'))

## Survival Forests (RSFs) ----
ranger_tuning_res = readRDS(file = paste0(res_path, 'ranger_tuning_res.rds'))
ranger_bootci = readRDS(file = paste0(res_path, 'ranger_bootci.rds'))

# Surrogates ----

## coxnet used KM (GP)
coxnet_tuning_res$coxnet_bayes_at$tuner$surrogate$model
## RSFs used ranger surrogate (due to `splitrule` being categorical)
ranger_tuning_res$ranger_bayes_at$tuner$surrogate$model

# Tuning results ----
## CoxNet alpha ~ 1, same results
coxnet_tuning_res$coxnet_bayes_at$tuning_result
coxnet_tuning_res$coxnet_rand_at$tuning_result

## RSF - C-splitrule prevailed, same results
ranger_tuning_res$ranger_bayes_at$tuning_result
ranger_tuning_res$ranger_rand_at$tuning_result

# Train times ----
# Random search 3 times faster with no parallelization for CoxNet
coxnet_tuning_res$coxnet_bayes_at$timings['train']/60 # 30 min
coxnet_tuning_res$coxnet_rand_at$timings['train']/60 # 12 min

# BO is faster with RSFs (with implicit parallelization on ranger)
ranger_tuning_res$ranger_bayes_at$timings['train']/3600 # 1 hours
ranger_tuning_res$ranger_rand_at$timings['train']/3600 # 2 hours

# Plots ----
## CoxNet ----
### Surface plots
autoplot(coxnet_tuning_res$coxnet_bayes_at$tuning_instance,
  type = 'surface', trafo = TRUE) +
  labs(title = 'CoxNet - Bayesian Optimization', x = 'lambda', y = 'alpha')
ggsave(filename = 'img/bayes_vs_randsearch/coxnet_bo_surface.png', width = 5, height = 5, dpi = 450)
autoplot(coxnet_tuning_res$coxnet_rand_at$tuning_instance,
  type = 'surface', trafo = TRUE) +
  labs(title = 'CoxNet - Random Search', x = 'lambda', y = 'alpha')
ggsave(filename = 'img/bayes_vs_randsearch/coxnet_rs_surface.png', width = 5, height = 5, dpi = 450)

### performance plots
# no need more than 50 evaluations/batches
autoplot(coxnet_tuning_res$coxnet_bayes_at$tuning_instance, type = 'performance') +
  labs(title = 'CoxNet - Bayesian Optimization')
ggsave(filename = 'img/bayes_vs_randsearch/coxnet_bo_perf.png', width = 5, height = 5, dpi = 450)

autoplot(coxnet_tuning_res$coxnet_rand_at$tuning_instance, type = 'performance') +
  labs(title = 'CoxNet - Random Search')
ggsave(filename = 'img/bayes_vs_randsearch/coxnet_rs_perf.png', width = 5, height = 5, dpi = 450)

### CIs
## QQ-plots seems quite normal to me :)
plot(coxnet_bootci$bootci_coxnet_rand$boot_res)
plot(coxnet_bootci$bootci_coxnet_bayes$boot_res)

coxnet_rs_ci = coxnet_bootci$bootci_coxnet_rand$bootci_res$bca
coxnet_bo_ci = coxnet_bootci$bootci_coxnet_bayes$bootci_res$bca

coxnet_ci = tibble::tibble(name = c('CoxNet-RS', 'CoxNet-BO'),
  cindex = c(coxnet_bootci$bootci_coxnet_rand$boot_res$t0,
             coxnet_bootci$bootci_coxnet_bayes$boot_res$t0),
  low_ci = c(unname(coxnet_rs_ci[,4]), unname(coxnet_bo_ci[,4])),
  high_ci = c(unname(coxnet_rs_ci[,5]), unname(coxnet_bo_ci[,5])))

coxnet_ci %>% mutate(name = factor(name, levels = c('CoxNet-RS', 'CoxNet-BO'))) %>%
  ggplot(aes(x = name, y = cindex, color = name)) +
  geom_pointrange(aes(ymin = low_ci, ymax = high_ci)) +
  geom_errorbar(aes(ymin = low_ci, ymax = high_ci), width = 0.1) +
  theme_bw(base_size = 16) + ylim(0.2, 0.8) +
  labs(x = '', y = 'C-index', title = 'Bootstrap 95% CIs on test set') +
  theme(legend.position = 'none')
ggsave(filename = 'img/bayes_vs_randsearch/coxnet_bootci.png', width = 5, height = 5, dpi = 450)

## RSFs ----
### Hyperparameter plots
# marginal => color == batches (not informative?)
autoplot(ranger_tuning_res$ranger_bayes_at$tuning_instance, type = 'marginal')
autoplot(ranger_tuning_res$ranger_rand_at$tuning_instance, type = 'marginal')

# parameter => color == cindex (better)
autoplot(ranger_tuning_res$ranger_bayes_at$tuning_instance, type = 'parameter')
ggsave(filename = 'img/bayes_vs_randsearch/ranger_bo_parameter.png', width = 9, height = 6, dpi = 450)
autoplot(ranger_tuning_res$ranger_rand_at$tuning_instance, type = 'parameter')
ggsave(filename = 'img/bayes_vs_randsearch/ranger_rs_parameter.png', width = 9, height = 6, dpi = 450)

autoplot(ranger_tuning_res$ranger_bayes_at$tuning_instance, type = 'parallel')
ggsave(filename = 'img/bayes_vs_randsearch/ranger_bo_parallel.png', width = 8, height = 5, dpi = 450)
autoplot(ranger_tuning_res$ranger_rand_at$tuning_instance, type = 'parallel')
ggsave(filename = 'img/bayes_vs_randsearch/ranger_rs_parallel.png', width = 8, height = 5, dpi = 450)

# Choose pair of variables to make surface plots
autoplot(ranger_tuning_res$ranger_bayes_at$tuning_instance, type = 'surface',
  cols_x = c('num.trees', 'mtry.ratio'))
autoplot(ranger_tuning_res$ranger_rand_at$tuning_instance, type = 'surface',
  cols_x = c('num.trees', 'mtry.ratio'))
autoplot(ranger_tuning_res$ranger_bayes_at$tuning_instance, type = 'points',
  cols_x = c('num.trees', 'splitrule'))

### Performance plots
autoplot(ranger_tuning_res$ranger_bayes_at$tuning_instance, type = 'performance') +
  labs(title = 'RSFs - Bayesian Optimization') + ylim(c(0.58, 0.72))
ggsave(filename = 'img/bayes_vs_randsearch/ranger_bo_perf.png', width = 5, height = 5, dpi = 450)
autoplot(ranger_tuning_res$ranger_rand_at$tuning_instance, type = 'performance') +
  labs(title = 'RSFs - Random Search') + ylim(c(0.58, 0.72))
ggsave(filename = 'img/bayes_vs_randsearch/ranger_rs_perf.png', width = 5, height = 5, dpi = 450)

### CIs
plot(ranger_bootci$bootci_ranger_rand$boot_res)
plot(ranger_bootci$bootci_ranger_bayes$boot_res)

ranger_rs_ci = ranger_bootci$bootci_ranger_rand$bootci_res$bca
ranger_bo_ci = ranger_bootci$bootci_ranger_bayes$bootci_res$bca

ranger_ci = tibble::tibble(name = c('RSF-RS', 'RSF-BO'),
  cindex = c(ranger_bootci$bootci_ranger_rand$boot_res$t0,
             ranger_bootci$bootci_ranger_bayes$boot_res$t0),
  low_ci = c(unname(ranger_rs_ci[,4]), unname(ranger_bo_ci[,4])),
  high_ci = c(unname(ranger_rs_ci[,5]), unname(ranger_bo_ci[,5])))

ranger_ci %>% mutate(name = factor(name, levels = c('RSF-RS', 'RSF-BO'))) %>%
  ggplot(aes(x = name %>% factor(), y = cindex, color = name)) +
  geom_pointrange(aes(ymin = low_ci, ymax = high_ci)) +
  geom_errorbar(aes(ymin = low_ci, ymax = high_ci), width = 0.1) +
  theme_bw(base_size = 16) + ylim(0.2,0.8) +
  labs(x = '', y = 'C-index', title = 'Bootstrap 95% CIs on test set') +
  theme(legend.position = 'none')
ggsave(filename = 'img/bayes_vs_randsearch/ranger_bootci.png', width = 5, height = 5, dpi = 450)

# check if C-index on test set is correct
# check C-indexes on train set
if (FALSE) {
  library(testthat)
  res = readRDS(file = 'results/xgboost/xgboost_cox.rds')
  test_indx = res$test_indx

  tasks = readRDS(file = 'data/tasks.rds')
  task_mRNA = tasks$mRNA

  expect_equal(coxnet_tuning_res$coxnet_rand_at$
      predict(task_mRNA, row_ids = test_indx)$score(),
    coxnet_bootci$bootci_coxnet_rand$boot_res$t0)

  expect_equal(ranger_tuning_res$ranger_bayes_at$
      predict(task_mRNA, row_ids = test_indx)$score(),
    ranger_bootci$bootci_ranger_bayes$boot_res$t0)

  # overfitting on train set
  train_indx = setdiff(seq_len(task_mRNA$nrow), test_indx)
  intersect(test_indx, train_indx)

  coxnet_tuning_res$coxnet_rand_at$predict(task_mRNA, train_indx)$score()
  coxnet_tuning_res$coxnet_bayes_at$predict(task_mRNA, train_indx)$score()

  ranger_tuning_res$ranger_rand_at$predict(task_mRNA, train_indx)$score()
  ranger_tuning_res$ranger_bayes_at$predict(task_mRNA, train_indx)$score()
}
