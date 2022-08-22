# plot results from `bench_xgboost.R` script,
# so run that script first to generate the output files
library(mlr3verse)
library(dplyr)
library(ggplot2)

res_path = 'results/xgboost/'

# XGBoost (survival-cox, 3 HPs) ----
xgboost1_res = readRDS(file = paste0(res_path, 'xgboost1_res.rds'))

## Tuning results
xgboost1_res$xgboost_at1$tuning_result$surv.cindex
data.table::rbindlist(xgboost1_res$xgboost_at1$tuning_result$x_domain)

xgboost1_res$xgboost_at1$archive$data %>%
  as_tibble() %>%
  select(nrounds, eta, runtime_learners) %>%
  mutate(eta = exp(eta)) # more eta + larger nrounds => more execution time

## Train time
xgboost1_res$xgboost_at1$timings['train']/60 # 48 min
sum(xgboost1_res$xgboost_at1$archive$data$runtime_learners)/60 # 21 min
# some time spend on BO?

## Plots
### Performance plot
autoplot(xgboost1_res$xgboost_at1$tuning_instance, type = 'performance')

### Parameter plot
autoplot(xgboost1_res$xgboost_at1$tuning_instance, type = 'parameter', trafo = TRUE)

### Surface plot
autoplot(xgboost1_res$xgboost_at1$tuning_instance, type = 'surface',
  cols_x = c('nrounds', 'x_domain_eta'))

# XGBoost (survival-cox, 6 HPs) ----
xgboost2_res = readRDS(file = paste0(res_path, 'xgboost2_res.rds'))

## Tuning results
xgboost2_res$xgboost_at2$tuning_result$surv.cindex
data.table::rbindlist(xgboost2_res$xgboost_at2$tuning_result$x_domain)

## Train time
xgboost2_res$xgboost_at2$timings['train']/3600 # 2.8h
sum(xgboost2_res$xgboost_at2$archive$data$runtime_learners)/3600 # 2.19h
# so some time spend on BO?

## Plots
### Performance plot
autoplot(xgboost2_res$xgboost_at2$tuning_instance, type = 'performance')

### Parameter plot
# regularization parameters changed the parameter plot landscape, allowing more
# configurations to be with okay performance?
# (alpha - L1, lambda = L2) - min_child_weight was less sensitive
autoplot(xgboost2_res$xgboost_at2$tuning_instance, type = 'parameter', trafo = TRUE)
# `min_child_weight` is close to default (don't tune?)
autoplot(xgboost2_res$xgboost_at2$tuning_instance, type = 'parameter', trafo = TRUE, cols_x = 'min_child_weight')

### Surface plot
autoplot(xgboost2_res$xgboost_at2$tuning_instance, type = 'surface',
  cols_x = c('nrounds', 'x_domain_eta'))

# XGBoost (survival-aft, 5 HPs) ----
xgboost3_res = readRDS(file = paste0(res_path, 'xgboost3_res.rds'))

## Tuning results
xgboost3_res$xgboost_at3$tuning_result$surv.cindex
data.table::rbindlist(xgboost3_res$xgboost_at3$tuning_result$x_domain)

## Train time
xgboost3_res$xgboost_at3$timings['train']/3600 # 2.6h
sum(xgboost3_res$xgboost_at3$archive$data$runtime_learners)/3600 # 2.53h
# so less time spend on BO compared to survival:cox?

## Plots
### Performance plot
# more stable
autoplot(xgboost3_res$xgboost_at3$tuning_instance, type = 'performance')

### Parameter plot
# logistic aft loss function seems a better choice for this dataset
autoplot(xgboost3_res$xgboost_at3$tuning_instance, type = 'parameter', trafo = TRUE)
autoplot(xgboost3_res$xgboost_at3$tuning_instance, type = 'parameter', trafo = TRUE, cols_x = 'aft_loss_distribution')

### Surface plot
# more concentrated nrounds values
autoplot(xgboost3_res$xgboost_at3$tuning_instance, type = 'surface',
  cols_x = c('nrounds', 'x_domain_eta'))

# XGBoost (survival-aft, 8 HPs) ----
xgboost4_res = readRDS(file = paste0(res_path, 'xgboost4_res.rds'))

## Tuning results
xgboost4_res$xgboost_at4$tuning_result$surv.cindex
data.table::rbindlist(xgboost4_res$xgboost_at4$tuning_result$x_domain)

## Train time
xgboost4_res$xgboost_at4$timings['train']/3600 # 1.84h
sum(xgboost4_res$xgboost_at4$archive$data$runtime_learners)/3600 # 1.76h
# so less time spend on BO compared to survival:cox?

## Plots
### Performance plot
# c-index = 0.5 sometimes...
autoplot(xgboost4_res$xgboost_at4$tuning_instance, type = 'performance')

### Parameter plot
# now more normal or extreme aft loss function seems best
autoplot(xgboost4_res$xgboost_at4$tuning_instance, type = 'parameter', trafo = TRUE)
autoplot(xgboost4_res$xgboost_at4$tuning_instance, type = 'parameter', trafo = TRUE, cols_x = 'aft_loss_distribution')
# `min_child_weight` is close to default (don't tune?)
autoplot(xgboost4_res$xgboost_at4$tuning_instance, type = 'parameter', trafo = TRUE, cols_x = 'min_child_weight')

### Surface plot
# all seem close to good prediction - alpha and lamba role?
autoplot(xgboost4_res$xgboost_at4$tuning_instance, type = 'surface',
  cols_x = c('nrounds', 'x_domain_eta'))
autoplot(xgboost4_res$xgboost_at4$tuning_instance, type = 'surface',
  cols_x = c('alpha', 'lambda'))
autoplot(xgboost4_res$xgboost_at4$tuning_instance, type = 'points',
  cols_x = c('x_domain_alpha', 'x_domain_lambda')) +
  scale_y_continuous(trans = scales::log10_trans()) +
  scale_x_continuous(trans = scales::log10_trans())

# Confidence Interval plot ----
## QQ-plots seems quite normal to me :)
plot(xgboost1_res$boot_ci1$boot_res)
plot(xgboost2_res$boot_ci2$boot_res)
plot(xgboost3_res$boot_ci3$boot_res)
plot(xgboost4_res$boot_ci4$boot_res)

xgboost1_ci = xgboost1_res$boot_ci1$bootci_res$bca
xgboost2_ci = xgboost2_res$boot_ci2$bootci_res$bca
xgboost3_ci = xgboost3_res$boot_ci3$bootci_res$bca
xgboost4_ci = xgboost4_res$boot_ci4$bootci_res$bca

xgboost_ci = tibble::tibble(name =
  c('XGBoost-Cox-3HPs', 'XGBoost-Cox-6HPs', 'XGBoost-AFT-5HPs', 'XGBoost-AFT-8HPs'),
  # Best C-index on train set after tuning
  cindex_train = c(xgboost1_res$xgboost_at1$tuning_result$surv.cindex,
    xgboost2_res$xgboost_at2$tuning_result$surv.cindex,
    xgboost3_res$xgboost_at3$tuning_result$surv.cindex,
    xgboost4_res$xgboost_at4$tuning_result$surv.cindex),
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
  geom_point(aes(x = name, y = cindex_train, size = 1)) +
  geom_hline(yintercept = 0.5, linetype = 'dotted', color = 'red') +
  theme_bw(base_size = 16) + ylim(0.2, 0.8) +
  labs(x = '', y = 'C-index', title = 'Bootstrap 95% CIs on test set + best C-index on trained set') +
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
