library(mlr3verse)
library(mlr3proba)
library(mlr3viz)
suppressMessages(library(dplyr))
library(readr)
library(ggplot2)
library(mlr3benchmark)
suppressMessages(library(rstatix))
library(ggpubr)
library(metR)

# bench_repeated_cv.R ----
# 4 untuned learners repeated-CV
perf_res = readRDS(file = 'results/perf_res_repeated_cv.rds')

## Performance Boxplot ----
perf_res %>%
  mutate(learner_id = factor(learner_id,
    levels = c('CoxLasso', 'Survival Tree', 'Survival Forest', 'XGBoost Survival Learner'))) %>%
  ggplot(aes(x = learner_id, y = surv.cindex, fill = learner_id)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(. ~ task_id) +
  mlr3viz::theme_mlr3(x.text.angle = 45) + xlab('') + ylab('C-index')
ggsave(filename = 'img/4learners_miRNA_mRNA_repeatedCV.png', width = 7, height = 5, dpi = 300)

# Is there are a significant difference in the rankings of the learners over all the tasks?
aggr_res = perf_res %>%
  group_by(task_id, learner_id) %>%
  summarise(c_index = mean(surv.cindex), .groups = 'drop') %>%
  mutate(task_id = factor(task_id), learner_id = factor(learner_id))

ba = BenchmarkAggr$new(aggr_res)
ba$friedman_test() # No!

# bench_glmnet_nestedCV.R ----
# glmnet Nested-CV + tuning
perf_res = readRDS(file = 'results/perf_res_glmnet_nestedCV.rds')

## Performance Boxplot ----
task_comp = list(c('CNA','mRNA'), c('Methyl','mRNA'), c('miRNA', 'mRNA'))

stat_test = perf_res %>%
  rstatix::wilcox_test(formula = surv.cindex ~ task_id, comparisons = task_comp) %>%
  rstatix::add_significance('p') %>%
  rstatix::add_y_position()

perf_res %>%
  ggplot(aes(x = task_id, y = surv.cindex)) +
  geom_boxplot(aes(fill = task_id), show.legend = FALSE) +
  scale_fill_brewer(palette = 'Set1') +
  mlr3viz::theme_mlr3(x.text.angle = 45) +
  xlab('') + ylab('C-index') + ggtitle('Tuned CoxNet (Nested-CV results)') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggpubr::stat_compare_means(hjust = -1.5) +
  ggpubr::stat_pvalue_manual(stat_test, label = '{p} ({p.signif})')
ggsave(filename = 'img/glmnet_4tasks_nestedCV.png', width = 7, height = 5, dpi = 300)

## Hyperparameter plot (a, lambda) ----
hpo_res = readRDS(file = 'results/hpo_res_glmnet_nestedCV.rds')
hp_data = hpo_res %>%
  filter(task_id == 'mRNA') %>%
  select(alpha, x_domain_lambda, surv.cindex) %>%
  rename(lambda = x_domain_lambda)
print(summary(hp_data)) # lambda in log scale

# Fill in a parameter grid
grid = interp::interp(x = hp_data$alpha, y = hp_data$lambda,
  nx = 100, ny = 100, z = hp_data$surv.cindex)
griddf = subset(data.frame(x = rep(grid$x, nrow(grid$z)),
  y = rep(grid$y, each = ncol(grid$z)),
  z = as.numeric(grid$z)),
  !is.na(z))

scale_fill_brewer_discretised = metR::as.discretised_scale(scale_fill_distiller)

griddf %>% ggplot(aes(x, y, z = z)) +
  metR::geom_contour_fill(aes(fill = stat(level))) +
  scale_fill_brewer_discretised(name = 'C-index', palette = 'RdYlBu', direction = 1) +
  labs(x = 'alpha', y = 'lambda') +
  geom_point(data = hp_data, aes(x = alpha, y = lambda), inherit.aes = F, size = 0.01) +
  theme_classic()
ggsave(filename = 'img/glmnet_hp_plot.png', width = 6, height = 5, dpi = 300)

# bench_nestedCV_v2.R ----
# glmnet,rpart,ranger Nested-CV + tuning
perf_res = readRDS(file = 'results/perf_res_nestedCV_v2.rds')

# better learner names
perf_res$learner_id = gsub(pattern = '\\.tuned', replacement = '', perf_res$learner_id)

perf_res = perf_res %>% mutate(learner_id = case_when(
  learner_id == 'SurvivalTree' ~ 'Survival Tree',
  learner_id == 'SurvivalForest' ~ 'Survival Forest',
  TRUE ~ learner_id
))

## Performance Boxplot (C-index) ----
perf_res %>%
  mutate(learner_id = factor(learner_id,
    levels = c('Survival Tree', 'CoxNet', 'Survival Forest'))) %>%
  ggplot(aes(x = learner_id, y = surv.cindex, fill = learner_id)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(. ~ task_id) +
  mlr3viz::theme_mlr3(x.text.angle = 45) + xlab('') + ylab('C-index')
ggsave(filename = 'img/bench_nestedCV_v2/cindex_boxplot.png', width = 7, height = 5, dpi = 300)

# get the boxplot colors for align color in later plots
hue_colors = scale_color_hue()$palette(n = 3)

## Stat. significance (C-index) ----
aggr_res = perf_res %>%
  mutate(learner_id = case_when(
    learner_id == 'Survival Tree' ~ 'Tree',
    learner_id == 'Survival Forest' ~ 'Forest',
    TRUE ~ learner_id
  )) %>%
  group_by(task_id, learner_id) %>%
  summarise(c_index = mean(surv.cindex), .groups = 'drop') %>%
  mutate(task_id = factor(task_id), learner_id = factor(learner_id))

ba = BenchmarkAggr$new(aggr_res)
ba$friedman_test() # Barely!

ba$friedman_posthoc(meas = 'c_index')
autoplot(ba, type = 'fn', meas = 'c_index') + theme(aspect.ratio = 0.5)
ggsave(filename = 'img/bench_nestedCV_v2/cindex_friedman_nemenyi.png', width = 3, height = 3, dpi = 300)

p = autoplot(ba, type = 'cd', meas = 'c_index', minimize = FALSE, style = 2) +
  theme(panel.background = element_rect(fill = 'white'),
        plot.background = element_rect(fill = 'white')) +
  scale_color_manual(values = c('Tree' = hue_colors[1],
    'Forest' = hue_colors[3], 'CoxNet' = hue_colors[2]))
p$layers[[6]] = NULL # manually remove the line segment below "Critical Difference"
p
ggsave(filename = 'img/bench_nestedCV_v2/cindex_cdplot.png', width = 4, height = 2, dpi = 300)

## Performance Boxplot (Integrated Brier score) ----
perf_res %>%
  mutate(learner_id = factor(learner_id,
    levels = c('Survival Tree', 'CoxNet', 'Survival Forest'))) %>%
  ggplot(aes(x = learner_id, y = surv.graf, fill = learner_id)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(. ~ task_id) +
  mlr3viz::theme_mlr3(x.text.angle = 45) + xlab('') + ylab('Integrated Brier Score')
ggsave(filename = 'img/bench_nestedCV_v2/brier_boxplot.png', width = 7, height = 5, dpi = 300)

## Stat. significance (Integrated Brier score) ----
aggr_res_brier = perf_res %>%
  mutate(learner_id = case_when(
    learner_id == 'Survival Tree' ~ 'Tree',
    learner_id == 'Survival Forest' ~ 'Forest',
    TRUE ~ learner_id
  )) %>%
  group_by(task_id, learner_id) %>%
  summarise(int_brier_score = mean(surv.graf), .groups = 'drop') %>%
  mutate(task_id = factor(task_id), learner_id = factor(learner_id))

ba_brier = BenchmarkAggr$new(aggr_res_brier)
ba_brier$friedman_test() # Barely!
ba_brier$friedman_posthoc(meas = 'int_brier_score') # same as C-index results

p = autoplot(ba_brier, type = 'cd', meas = 'int_brier_score', minimize = TRUE, style = 2) +
  theme(panel.background = element_rect(fill = 'white'),
    plot.background = element_rect(fill = 'white')) +
  scale_color_manual(values = c('Tree' = hue_colors[1],
    'Forest' = hue_colors[3], 'CoxNet' = hue_colors[2]))
p$layers[[6]] = NULL # remove the line segment below "Critical Difference"
p
ggsave(filename = 'img/bench_nestedCV_v2/brier_cdplot.png', width = 4, height = 2, dpi = 300) # identical as C-index

## Performance Boxplot (Integrated Log Loss) ----
perf_res %>%
  mutate(learner_id = factor(learner_id,
    levels = c('Survival Tree', 'CoxNet', 'Survival Forest'))) %>%
  ggplot(aes(x = learner_id, y = surv.intlogloss, fill = learner_id)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(. ~ task_id) +
  mlr3viz::theme_mlr3(x.text.angle = 45) + xlab('') + ylab('Integrated Log Loss')
ggsave(filename = 'img/bench_nestedCV_v2/logloss_boxplot.png', width = 7, height = 5, dpi = 300) # qualitative same results as the integrated brier score

## Stat. significance (Integrated Log Loss) ----
aggr_res_logloss = perf_res %>%
  mutate(learner_id = case_when(
    learner_id == 'Survival Tree' ~ 'Tree',
    learner_id == 'Survival Forest' ~ 'Forest',
    TRUE ~ learner_id
  )) %>%
  group_by(task_id, learner_id) %>%
  summarise(int_logloss = mean(surv.intlogloss), .groups = 'drop') %>%
  mutate(task_id = factor(task_id), learner_id = factor(learner_id))

ba_logloss = BenchmarkAggr$new(aggr_res_logloss)
ba_logloss$friedman_test() # Barely!
ba_logloss$friedman_posthoc() # same as C-index and int. brier score results

p = autoplot(ba_logloss, type = 'cd', meas = 'int_logloss', minimize = TRUE, style = 2) +
  theme(panel.background = element_rect(fill = 'white'),
    plot.background = element_rect(fill = 'white')) +
  scale_color_manual(values = c('Tree' = hue_colors[1],
    'Forest' = hue_colors[3], 'CoxNet' = hue_colors[2]))
p$layers[[6]] = NULL # remove the line segment below "Critical Difference"
p
ggsave(filename = 'img/bench_nestedCV_v2/logloss_cdplot.png', width = 4, height = 2, dpi = 300) # identical as C-index and int. brier score results


