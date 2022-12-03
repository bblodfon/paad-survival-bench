# Investigate model extrapolation in the mRNA dataset
# i.e. if test data is "too different" from train data
# Result: didn't found any evidence of this
library(mlr3verse)
library(tidyverse)

if (!dir.exists('img/model_extr/')) {
  dir.create('img/model_extr/')
}

# mRNA Task ----
# ~10000 features
tasks = readRDS(file = 'data/tasks.rds')
task_mRNA = tasks$mRNA

# get a benchmark xgboost result with the train and test indexes
res = readRDS(file = 'results/xgboost/xgboost_cox.rds')
intersect(res$train_indx, res$test_indx)

grp1 = sapply(1:task_mRNA$nrow,
  function(index) if (index %in% res$train_indx) "train" else "test")
all(sort(res$train_indx) == which(grp1 == 'train')) # check
all(sort(res$test_indx)  == which(grp1 == 'test')) # check

# PCA ----
pca = po('pca', center = FALSE, scale. = FALSE, rank. = 50) # already standardized
pca$param_set
pca$is_trained

task_pca = pca$train(list(task_mRNA))[[1L]]

# `$predict` using all data/task has the same output as `$train` above
pca$predict(list(task_mRNA$clone(deep = TRUE)$filter(1)))[[1L]]$data()
pca$predict(list(task_mRNA))[[1L]]$data()

p1 = task_pca$data() %>%
  add_column(group = grp1) %>%
  ggplot(aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  labs(title = 'PCA (mRNA)', subtitle = 'PCA performed on all data') +
  theme_bw(base_size = 14)
p1

# PCA (separately train/test) ----
pca = po('pca', center = FALSE, scale. = FALSE, rank. = 2)
pca$is_trained

train_mRNA_pca = pca$train(list(task_mRNA$clone(deep = TRUE)$filter(rows = res$train_indx)))[[1L]]
test_mRNA_pca  = pca$predict(list(task_mRNA$clone(deep = TRUE)$filter(rows = res$test_indx)))[[1L]]

grp = c(rep('train', length(res$train_indx)), rep('test', length(res$test_indx)))
tbl = bind_rows(
  train_mRNA_pca$data(cols = c('PC1', 'PC2')),
  test_mRNA_pca $data(cols = c('PC1', 'PC2'))) %>%
  add_column(group = grp)

p2 = tbl %>%
  ggplot(aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  labs(title = 'PCA (mRNA)', subtitle = 'PCA performed on train set') +
  theme_bw(base_size = 14)
p2

ggpubr::ggexport(plotlist = list(p1, p2), filename = 'img/model_extr/mRNA_pca.png',
  ncol = 2, nrow = 1, res = 450, width = 3700, height = 2000)

# Model extrapolation check ----
## Investigate the clinical data => are patients in the test set "different"
## from the ones in the train set? (NO)
clinical_var_mat = readRDS(file = 'data/clinical_var_mat.rds')

clinical_data_tbl = clinical_var_mat %>% as_tibble(rownames = 'patient_id') %>%
  rename(status = vital_status) %>%
  readr::type_convert() %>%
  mutate(time = ifelse(is.na(days_to_death), days_to_last_followup, days_to_death)) %>%
  rename(age = years_to_birth) %>%
  select(-c('patient_id', 'days_to_death', 'days_to_last_followup',
    'tumor_tissue_site', 'histological_type', # same for all
    'date_of_initial_pathologic_diagnosis',
    'number_pack_years_smoked', 'year_of_tobacco_smoking_onset', # many NAs, and
    # it may be not a risk factor for pancreatic cancer
    'race', 'ethnicity' # most are white, non-hispanic/latinos
    )) %>%
  filter(time > 0) %>% # remove one patient as was done in `prepare_tasks.R`
  add_column(group = grp1)

## convert to task for ease of use
paad_task = mlr3proba::as_task_surv(
  x = clinical_data_tbl,
  time = 'time', event = 'status', id = 'Clinical'
)
paad_task
paad_task$missings()

## Survival (KM)
autoplot(paad_task, rhs = 'group')
ggsave(filename = 'img/model_extr/survival.png', width = 6, height = 5, dpi = 350)

## KM plot of censored patients in train and test groups
paad_task_cens = mlr3proba::as_task_surv(
  x = clinical_data_tbl %>% mutate(status = 1 - status),
  time = 'time', event = 'status', id = 'Clinical-Cens'
)
autoplot(paad_task_cens, rhs = 'group')
ggsave(filename = 'img/model_extr/cens_survival.png', width = 6, height = 5, dpi = 350)

# check the time distributions of dead vs censored patients on all the dataset
paad_task$data(cols = c('time', 'status')) %>%
  mutate(status = as.factor(status)) %>%
  ggplot(aes(x = time, color = status)) +
  geom_density() +
  labs(title = 'Uncensored vs Censoring distribution') +
  theme_bw(base_size = 14)
ggsave(filename = 'img/model_extr/distr_uncens_vs_cens.png', width = 6, height = 5, dpi = 350)

paad_task$data(cols = c('time', 'status'))[status == 0, time]
plot(density(values))

## Time
clinical_data_tbl %>%
  ggplot(aes(x = group, y = time, color = group)) +
  geom_boxplot() +
  theme_bw(base_size = 14)
ggsave(filename = 'img/model_extr/time.png', width = 6, height = 5, dpi = 350)

## Status
clinical_data_tbl %>%
  group_by(group, status) %>%
  summarise(n=n(), .groups = "drop_last") %>%
  mutate(rel.freq = 100*(n/sum(n))) %>%
  mutate(status = factor(status)) %>%
  ggplot(aes(x = group, fill = status, y = rel.freq)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  ylim(0,80) +
  labs(y = 'Relative Frequency (%)') +
  theme_bw(base_size = 14)
ggsave(filename = 'img/model_extr/status.png', width = 6, height = 5, dpi = 350)

## age
clinical_data_tbl %>%
  ggplot(aes(x = group, y = age, color = group)) +
  geom_boxplot() +
  theme_bw(base_size = 14)
ggsave(filename = 'img/model_extr/age.png', width = 6, height = 5, dpi = 350)

## number of lymph nodes
clinical_data_tbl %>%
  ggplot(aes(x = group, y = number_of_lymph_nodes, color = group)) +
  geom_boxplot() +
  theme_bw(base_size = 14)
ggsave(filename = 'img/model_extr/lymph_nodes.png', width = 6, height = 5, dpi = 350)

## gender
clinical_data_tbl %>%
  group_by(group, gender) %>%
  summarise(n=n(), .groups = "drop_last") %>%
  mutate(rel.freq = 100*(n/sum(n))) %>%
  ggplot(aes(x = group, fill = gender, y = rel.freq)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  ylim(0,80) +
  labs(y = 'Relative Frequency (%)') +
  theme_bw(base_size = 14)
ggsave(filename = 'img/model_extr/gender.png', width = 6, height = 5, dpi = 350)

## pathologic_stage
clinical_data_tbl %>%
  group_by(group, pathologic_stage) %>%
  summarise(n=n(), .groups = "drop_last") %>%
  mutate(rel.freq = 100*(n/sum(n))) %>%
  ggplot(aes(x = group, fill = pathologic_stage, y = rel.freq)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  ylim(0,80) +
  labs(y = 'Relative Frequency (%)') +
  theme_bw(base_size = 14)
ggsave(filename = 'img/model_extr/path_stage.png', width = 6, height = 5, dpi = 350)

## pathology_T_stage
clinical_data_tbl %>%
  group_by(group, pathology_T_stage) %>%
  summarise(n=n(), .groups = "drop_last") %>%
  mutate(rel.freq = 100*(n/sum(n))) %>%
  ggplot(aes(x = group, fill = pathology_T_stage, y = rel.freq)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  #ylim(0,80) +
  labs(y = 'Relative Frequency (%)') +
  theme_bw(base_size = 14)
ggsave(filename = 'img/model_extr/stage_t.png', width = 6, height = 5, dpi = 350)

## pathology_N_stage
clinical_data_tbl %>%
  group_by(group, pathology_N_stage) %>%
  summarise(n=n(), .groups = "drop_last") %>%
  mutate(rel.freq = 100*(n/sum(n))) %>%
  ggplot(aes(x = group, fill = pathology_N_stage, y = rel.freq)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  ylim(0,80) +
  labs(y = 'Relative Frequency (%)') +
  theme_bw(base_size = 14)
ggsave(filename = 'img/model_extr/stage_n.png', width = 6, height = 5, dpi = 350)

## pathology_M_stage
clinical_data_tbl %>%
  group_by(group, pathology_M_stage) %>%
  summarise(n=n(), .groups = "drop_last") %>%
  mutate(rel.freq = 100*(n/sum(n))) %>%
  ggplot(aes(x = group, fill = pathology_M_stage, y = rel.freq)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  ylim(0,80) +
  labs(y = 'Relative Frequency (%)') +
  theme_bw(base_size = 14)
ggsave(filename = 'img/model_extr/stage_m.png', width = 6, height = 5, dpi = 350)

## radiation_therapy
clinical_data_tbl %>%
  group_by(group, radiation_therapy) %>%
  summarise(n=n(), .groups = "drop_last") %>%
  mutate(rel.freq = 100*(n/sum(n))) %>%
  ggplot(aes(x = group, fill = radiation_therapy, y = rel.freq)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  ylim(0,80) +
  labs(y = 'Relative Frequency (%)') +
  theme_bw(base_size = 14)
ggsave(filename = 'img/model_extr/radiotherapy.png', width = 6, height = 5, dpi = 350)

## residual tumor
clinical_data_tbl %>%
  group_by(group, residual_tumor) %>%
  summarise(n=n(), .groups = "drop_last") %>%
  mutate(rel.freq = 100*(n/sum(n))) %>%
  ggplot(aes(x = group, fill = residual_tumor, y = rel.freq)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  ylim(0,80) +
  labs(y = 'Relative Frequency (%)') +
  theme_bw(base_size = 14)
ggsave(filename = 'img/model_extr/residual_tumor.png', width = 6, height = 5, dpi = 350)
