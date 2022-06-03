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

# 4 untuned learners repeated-CV ----
# see `scripts/bench_repeated_cv.R`
perf_res = readRDS(file = 'results/perf_res_repeated_cv.rds')

## Performance Boxplot
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
  summarise(c_index = mean(surv.cindex), .groups = "drop") %>%
  mutate(task_id = factor(task_id), learner_id = factor(learner_id))

ba = BenchmarkAggr$new(aggr_res)
ba$friedman_test() # No!

# glmnet Nested-CV + tuning ----
# see `bench_glmnet_nestedCV.R`
perf_res = readRDS(file = 'results/perf_res_glmnet_nestedCV.rds')

## Performance Boxplot
task_comp = list(c('CNA','mRNA'), c('Methyl','mRNA'), c('miRNA', 'mRNA'))

stat_test = perf_res %>%
  rstatix::wilcox_test(formula = surv.cindex ~ task_id, comparisons = task_comp) %>%
  rstatix::add_significance("p") %>%
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

## plot (alpha, lambda) hyperparameter plot
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
