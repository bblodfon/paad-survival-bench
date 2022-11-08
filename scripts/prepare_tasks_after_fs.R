library(mlr3verse)
library(mlr3proba)
library(tidyverse)
source('scripts/helpers.R')

# Load Data ----
#' Results from `fs_analysis.R`
fs_consensus_res = readRDS(file = 'results/fs/fs_consensus_res.rds')
#' Results from `fs_{Modality}.R`
fs_res_mRNA   = readRDS(file = 'results/fs/mRNA/fs_results.rds')
fs_res_miRNA  = readRDS(file = 'results/fs/miRNA/fs_results.rds')
fs_res_CNA    = readRDS(file = 'results/fs/CNA/fs_results.rds')
fs_res_Methyl = readRDS(file = 'results/fs/Methyl/fs_results.rds')
fs_res = list(mRNA = fs_res_mRNA$res, miRNA = fs_res_miRNA$res,
              CNA = fs_res_CNA$res, Methyl = fs_res_Methyl$res)
#' Tasks from `prepare_tasks.R`
tasks = readRDS(file = 'data/tasks.rds')
nfeats_mRNA = length(tasks$mRNA$feature_names)
nfeats_miRNA = length(tasks$miRNA$feature_names)
nfeats_CNA = length(tasks$CNA$feature_names)
nfeats_Methyl = length(tasks$Methyl$feature_names)

# Select #features per modality ----
#' Decide how many (top-selected) features to have per modality

#' 1st way: Using the consensus features results
data_stats1 = lapply(fs_consensus_res, function(l) {
  nfeats_total = nrow(l)

  # top 1%
  nfeats_top1per = ceiling(0.01 * nfeats_total)
  # top 3%
  nfeats_top3per = ceiling(0.03 * nfeats_total)
  # top 5%
  nfeats_top5per = ceiling(0.05 * nfeats_total)

  # freq_cutoff => 0.5
  nfeats_freq_cutoff0.5 = l %>% filter(freq > 0.5) %>% nrow()
  # freq_cutoff => 0.6
  nfeats_freq_cutoff0.6 = l %>% filter(freq > 0.6) %>% nrow()
  # freq_cutoff => 0.7
  nfeats_freq_cutoff0.7 = l %>% filter(freq > 0.7) %>% nrow()

  tibble::tibble(
    top10 = 10,
    top15 = 15,
    top1per = nfeats_top1per,
    top3per = nfeats_top3per,
    top5per = nfeats_top5per,
    cutoff0.5 = nfeats_freq_cutoff0.5,
    cutoff0.6 = nfeats_freq_cutoff0.6,
    cutoff0.7 = nfeats_freq_cutoff0.7
  )
}) %>% dplyr::bind_rows(.id = 'data_type')

row_sums = data_stats1 %>%
  select(-data_type) %>%
  colSums() %>%
  as_tibble_row() %>%
  add_column(data_type = 'all', .before = 1)

data_stats1 = dplyr::bind_rows(data_stats1, row_sums)
data_stats1 # could choose cutoff for example - no golden rule exists here!

# e.g. get features with cutoff = 0.6
fs_consensus_res$mRNA %>%
  filter(freq > 0.6) %>%
  pull(feat_name)

#' 2nd way: choose based on average performance score (C-index) and #selected features
data_stats2 = lapply(fs_res, function(res_tbl) {
  res_tbl %>%
    filter(method == 'rfe') %>%
    mutate(n_features = sapply(selected_features, length)) %>%
    summarise(avg_score = mean(score_rsmp), avg_feat_set_len = ceiling(mean(n_features)))
}) %>% dplyr::bind_rows(.id = 'data_type')
data_stats2

#' add more info to `data_stats2`
nfeats_selected = lapply(fs_consensus_res, nrow) %>% unlist()
nfeats_selected #features selected per modality during FS
stopifnot(names(nfeats_selected) == data_stats2$data_type)

total_nfeats = c(nfeats_mRNA, nfeats_miRNA, nfeats_CNA, nfeats_Methyl)
names(total_nfeats) = data_stats2$data_type
total_nfeats # total number of features per data type

#' play with some derived scores for choosing final #features to select
data_stats3 = data_stats2 %>%
  add_column(nfeats_selected) %>%
  add_column(total_nfeats) %>%
  mutate(avg_feat_set_per = avg_feat_set_len/total_nfeats) %>%
  #mutate(nfeats_selected_per = nfeats_selected/total_nfeats) %>%
  #mutate(s1 = avg_score/avg_feat_set_len, s1_norm = s1/max(s1)) %>%
  #mutate(s2 = avg_score/nfeats_selected, s2_norm = s2/max(s2)) %>%
  #mutate(s3 = avg_score*nfeats_selected, s3_norm = s3/max(s3)) %>%
  mutate(score = avg_score*avg_feat_set_len, score_norm = score/max(score)) # seems better

data_stats3 %>%
  arrange(desc(avg_score)) %>%
  select(-c(total_nfeats, avg_feat_set_per, avg_feat_set_per))

#' choose #features based on `max_nfeats` value and `score_norm`
max_nfeats = 13 # the modality which has the highest score will have this #features
nfeat_tbl = data_stats3 %>%
  select(data_type, avg_score, avg_feat_set_len, score_norm) %>%
  mutate(nfeats = ceiling(max_nfeats * score_norm)) %>%
  arrange(desc(score_norm))
nfeat_tbl
nfeat_tbl %>% summarise(sum(nfeats))

# Subset tasks to selected features ----
task_mRNA = tasks$mRNA

feat_list = list()
for (data_modality in names(fs_consensus_res)) {
  # selected feature frequency in descending order
  freq_tbl = fs_consensus_res[[data_modality]]
  # get #nfeats to extract from the frequency table
  nfeats = nfeat_tbl %>% filter(data_type == data_modality) %>% pull(nfeats)
  # get the feature names
  feat_list[[data_modality]] = freq_tbl %>% slice(1:nfeats) %>% pull(feat_name)
}

feat_list

#' We can't just subset a task, e.g. `tasks$mRNA` keeps the backend and takes
#' too much space, so we re-create them
task_list = list()
for (data_modality in names(feat_list)) {
  print(data_modality)

  task = tasks[[data_modality]]
  features = feat_list[[data_modality]]

  task_list[[data_modality]] = mlr3proba::as_task_surv(
    x = task$select(features)$data(),
    time = 'time', event = 'status', id = data_modality
  )
}

#' clinical data task (see `prepare_clinical_data_task.R`)
task_clinical = readRDS(file = 'data/task_clinical.rds')

# scale to unit variance as the other omics datasets
pos = po('scale')
task_clinical = pos$train(list(task_clinical))[[1L]]
task_clinical_scaled$data()
var_flt$calculate(task_clinical)
var_flt$scores # ok

# add clinical dataset to the list of tasks (as first element)
task_list = c(Clinical = task_clinical, task_list)
length(task_list) # 5 data types

# Create multimodal tasks ----
#' Combine the tasks (datasets) in any possible way, creating all possible
#' subsets, e.g. clinical + mRNA + miRNA
tasks = get_task_powerset(task_list)
length(tasks) # 31 tasks

# Save tasks ----
saveRDS(tasks, file = 'data/tasks_fs.rds')
