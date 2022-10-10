# plot results from `glmboost_bench.R` script,
# so run that script first to generate the output files
library(mlr3verse)
library(tidyverse)

if (!dir.exists('img/glmboost/')) {
  dir.create('img/glmboost/')
}

# Task mRNA ----
task_mRNA = readRDS(file = 'data/tasks.rds')$mRNA

# Read result data ----
## GlmBoost Cox
res_cox = readRDS(file = 'results/glmboost/glmboost_CoxPH_mRNA.rds')
## GlmBoost AFT
res_aft = readRDS(file = 'results/glmboost/glmboost_AFT_mRNA.rds')
## GlmBoost C-Uno
res_cindex = readRDS(file = 'results/glmboost/glmboost_Cindex_mRNA.rds')

## verify train and test indexes are the same
all(res_cox$train_indx == res_aft$train_indx)
all(res_aft$train_indx == res_cindex$train_indx)
all(res_cox$test_indx == res_aft$test_indx)
all(res_aft$test_indx == res_cindex$test_indx)
train_indx = res_cox$train_indx
test_indx  = res_cox$test_indx

# GlmBoost - perf plot (Cox, mRNA only) ----
## test, train + train CV C-index for all hps
set1_colors = RColorBrewer::brewer.pal(n = 3, name = 'Set1')
res_cox$hpc_res %>% ggplot(aes(x = index)) +
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
ggsave(filename = 'img/glmboost/perf_mRNA_cox.png', width = 7, height = 5, dpi = 450)

res_aft$hpc_res %>% ggplot(aes(x = index)) +
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
ggsave(filename = 'img/glmboost/perf_mRNA_aft.png', width = 7, height = 5, dpi = 450)

res_cindex$hpc_res %>% ggplot(aes(x = index)) +
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
ggsave(filename = 'img/glmboost/perf_mRNA_cindex.png', width = 7, height = 5, dpi = 450)
