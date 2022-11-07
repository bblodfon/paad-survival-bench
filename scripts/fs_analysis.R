#'######################################################
#' Perform feature selection (FS) analysis regarding
#' stability, performance, #features selected using the
#' results of the scripts `fs_{datatype}.R`
#'
#' Main takeout is to not use GA (Genetic Algorithms)
#' since they are hard to configure and get stable
#' features with good performance.
#' RFE with RSFs worked fine.
#'######################################################
library(tidyverse)
library(mlr3verse)
library(stabm)
source('scripts/helpers.R')

# mRNA ----
#' see script `fs_mRNA.R`

## path to save plots
res_path = 'img/fs/mRNA'
if (!dir.exists(res_path)) {
  dir.create(res_path, recursive = TRUE)
}

## mRNA task
task_mRNA = readRDS(file = 'data/tasks.rds')$mRNA
nfeats_mRNA = length(task_mRNA$feature_names)

## mRNA eFS results
fs_results_mRNA = readRDS(file = 'results/fs/mRNA/fs_results.rds')
config_mRNA = fs_results_mRNA$config
res_mRNA    = fs_results_mRNA$res %>%
  mutate(n_features = sapply(selected_features, length))

## Performance plot ----
res_mRNA %>%
  ggplot(aes(x = method, y = score_rsmp, fill = learner_id)) +
  geom_boxplot() +
  labs(y = 'C-index (3 x 5-fold CV)', x = 'Wrapper-based FS method', title = 'mRNA') +
  #ylim(c(0.5,0.9)) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(res_path, '/perf.png'), width = 4.5, height = 5, dpi = 350)

## Number of features plot ----
res_mRNA %>%
  ggplot(aes(x = method, y = n_features, fill = learner_id)) +
  geom_boxplot() +
  labs(y = 'Number of selected features', x = 'Wrapper-based FS method', title = 'mRNA') +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(res_path, '/nfeat.png'), width = 4.5, height = 5, dpi = 350)

## Selection frequency plots ----
learners = unique(res_mRNA$learner_id)
data_list = list()
index = 1
for (wrapper_method in config_mRNA$wrapper_methods) {
  for (learner in learners) {
    selfeats_list = res_mRNA %>%
      filter(method == wrapper_method, learner_id == learner) %>%
      pull(selected_features)
    if (length(selfeats_list) > 0) {
      freq_tbl = get_consensus_features(selfeats_list)

      id = paste0(wrapper_method, '-', learner)
      p = feat_freq_barplot(freq_tbl, top_n = 10, title = paste0('mRNA (', id, ')'))
      ggsave(plot = p, filename = paste0(res_path, '/selfreq_', id, '.png'),
        width = 4.5, height = 5, dpi = 350)

      # calculate stability measures for each list of selected features
      data_list[[id]] = tibble(method = wrapper_method, learner_id = learner,
        jaccard = stabm::stabilityJaccard(features = selfeats_list, correction.for.chance = 'none'),
        nogueira = stabm::stabilityNogueira(features = selfeats_list, p = nfeats_mRNA))
    }
  }
}

stab_tbl = bind_rows(data_list)

## Stability assessment ----
stab_tbl %>%
  ggplot(aes(x = method, y = jaccard, fill = learner_id)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x = 'Wrapper-based FS method', y = 'Jaccard Similarity', title = 'mRNA') +
  ylim(c(0, 0.68)) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(res_path, '/jaccard_stab.png'), width = 4.5, height = 5, dpi = 350)

stab_tbl %>%
  ggplot(aes(x = method, y = nogueira, fill = learner_id)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x = 'Wrapper-based FS method', y = 'Nogueira Similarity', title = 'mRNA') +
  ylim(c(0, 0.68)) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(res_path, '/nogueira_stab.png'), width = 4.5, height = 5, dpi = 350)

## Consensus features and stability (RFE-RSFs) ----
rsf_lrns = stringr::str_subset(learners, 'rsf')
lrns_list = lapply(1:3, combinat::combn, x = rsf_lrns,
  simplify = FALSE) %>% unlist(recursive = FALSE)
lrns_list

#' `feat_ll` => list of list of features
feat_ll = lapply(lrns_list, function(lrns) {
  res_mRNA %>%
    filter(method == 'rfe', learner_id %in% lrns) %>%
    pull(selected_features)
})

data_list = list()
index = 1
for (index in 1:length(lrns_list)) {
  # learners combined
  print(lrns_list[[index]])

  feat_list = feat_ll[[index]]

  # Stability
  jaccard = stabm::stabilityJaccard(features = feat_list, correction.for.chance = 'none')
  nogueira = stabm::stabilityNogueira(features = feat_list, p = nfeats_mRNA)

  # Consensus features
  freq_tbl = get_consensus_features(feat_list)
  data_list[[index]] = list(
    id = paste0(str_remove(lrns_list[[index]], 'rsf_'), collapse = '-'),
    jaccard = jaccard,
    nogueira = nogueira,
    freq_tbl = freq_tbl)
}

rsfstab_tbl = lapply(data_list, function(l) { l$freq_tbl = NULL; l }) %>%
  dplyr::bind_rows()

rsfstab_tbl %>%
  mutate(id = forcats::fct_reorder(id, jaccard, .desc = TRUE)) %>%
  tidyr::pivot_longer(cols = c('jaccard', 'nogueira'), names_to = 'stab_measure') %>%
  ggplot(aes(x = id, y = value, fill = stab_measure)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x = 'RSF feature set combinations', y = 'Similarity score', title = 'mRNA') +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(res_path, '/rsfcmb_stab.png'), width = 6, height = 5, dpi = 350)

mRNA_cons_freq_tbl = data_list[[which(
  sapply(data_list, `[[`, 1) == 'cindex-logrank-maxstat'
)]]$freq_tbl
p = feat_freq_barplot(mRNA_cons_freq_tbl, top_n = 10,
  title = paste0('mRNA (RSF consensus)'))
ggsave(plot = p, filename = paste0(res_path, '/cons_selfreq.png'),
  width = 4.5, height = 5, dpi = 350)

# CNA ----
#' see script `fs_CNA.R`

## path to save plots
res_path = 'img/fs/CNA'
if (!dir.exists(res_path)) {
  dir.create(res_path, recursive = TRUE)
}

## CNA task
task_CNA = readRDS(file = 'data/tasks.rds')$CNA
nfeats_CNA = length(task_CNA$feature_names)

## CNA eFS results
fs_results_CNA = readRDS(file = 'results/fs/CNA/fs_results.rds')
config_CNA = fs_results_CNA$config
res_CNA    = fs_results_CNA$res %>%
  mutate(n_features = sapply(selected_features, length))

## Performance plot ----
res_CNA %>%
  ggplot(aes(x = method, y = score_rsmp, fill = learner_id)) +
  geom_boxplot() +
  labs(y = 'C-index (3 x 5-fold CV)', x = 'Wrapper-based FS method', title = 'CNA') +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(res_path, '/perf.png'), width = 4.5, height = 5, dpi = 350)

## Number of features plot ----
res_CNA %>%
  ggplot(aes(x = method, y = n_features, fill = learner_id)) +
  geom_boxplot() +
  labs(y = 'Number of selected features', x = 'Wrapper-based FS method', title = 'CNA') +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(res_path, '/nfeat.png'), width = 4.5, height = 5, dpi = 350)

## Selection frequency plots ----
learners = unique(res_CNA$learner_id)
data_list = list()
index = 1
for (wrapper_method in config_CNA$wrapper_methods) {
  for (learner in learners) {
    selfeats_list = res_CNA %>%
      filter(method == wrapper_method, learner_id == learner) %>%
      pull(selected_features)
    if (length(selfeats_list) > 0) {
      freq_tbl = get_consensus_features(selfeats_list)

      id = paste0(wrapper_method, '-', learner)
      p = feat_freq_barplot(freq_tbl, top_n = 10, title = paste0('CNA (', id, ')'))
      ggsave(plot = p, filename = paste0(res_path, '/selfreq_', id, '.png'),
        width = 4.5, height = 5, dpi = 350)

      # calculate stability measures for each list of selected features
      data_list[[id]] = tibble(method = wrapper_method, learner_id = learner,
        jaccard = stabm::stabilityJaccard(features = selfeats_list, correction.for.chance = 'none'),
        nogueira = stabm::stabilityNogueira(features = selfeats_list, p = nfeats_CNA))
    }
  }
}

stab_tbl = bind_rows(data_list)

## Stability assessment ----
stab_tbl %>%
  ggplot(aes(x = method, y = jaccard, fill = learner_id)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x = 'Wrapper-based FS method', y = 'Jaccard Similarity', title = 'CNA') +
  ylim(c(0, 0.72)) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(res_path, '/jaccard_stab.png'), width = 4.5, height = 5, dpi = 350)

stab_tbl %>%
  ggplot(aes(x = method, y = nogueira, fill = learner_id)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x = 'Wrapper-based FS method', y = 'Nogueira Similarity', title = 'CNA') +
  ylim(c(0, 0.72)) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(res_path, '/nogueira_stab.png'), width = 4.5, height = 5, dpi = 350)

## Consensus features and stability (RFE-RSFs) ----

#' `feat_ll` => list of list of features
feat_ll = lapply(lrns_list, function(lrns) {
  res_CNA %>%
    filter(method == 'rfe', learner_id %in% lrns) %>%
    pull(selected_features)
})

data_list = list()
index = 1
for (index in 1:length(lrns_list)) {
  # learners combined
  print(lrns_list[[index]])

  feat_list = feat_ll[[index]]

  # Stability
  jaccard = stabm::stabilityJaccard(features = feat_list, correction.for.chance = 'none')
  nogueira = stabm::stabilityNogueira(features = feat_list, p = nfeats_CNA)

  # Consensus features
  freq_tbl = get_consensus_features(feat_list)
  data_list[[index]] = list(
    id = paste0(str_remove(lrns_list[[index]], 'rsf_'), collapse = '-'),
    jaccard = jaccard,
    nogueira = nogueira,
    freq_tbl = freq_tbl)
}

rsfstab_tbl = lapply(data_list, function(l) { l$freq_tbl = NULL; l }) %>%
  dplyr::bind_rows()

rsfstab_tbl %>%
  mutate(id = forcats::fct_reorder(id, jaccard, .desc = TRUE)) %>%
  tidyr::pivot_longer(cols = c('jaccard', 'nogueira'), names_to = 'stab_measure') %>%
  ggplot(aes(x = id, y = value, fill = stab_measure)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x = 'RSF feature set combinations', y = 'Similarity score', title = 'CNA') +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(res_path, '/rsfcmb_stab.png'), width = 6, height = 5, dpi = 350)

CNA_cons_freq_tbl = data_list[[which(
  sapply(data_list, `[[`, 1) == 'cindex-logrank-maxstat'
)]]$freq_tbl
p = feat_freq_barplot(CNA_cons_freq_tbl, top_n = 10,
  title = paste0('CNA (RSF consensus)'))
ggsave(plot = p, filename = paste0(res_path, '/cons_selfreq.png'),
  width = 4.5, height = 5, dpi = 350)

# miRNA ----
#' see script `fs_miRNA.R`

## path to save plots
res_path = 'img/fs/miRNA'
if (!dir.exists(res_path)) {
  dir.create(res_path, recursive = TRUE)
}

## miRNA task
task_miRNA = readRDS(file = 'data/tasks.rds')$miRNA
nfeats_miRNA = length(task_miRNA$feature_names)

## miRNA eFS results
fs_results_miRNA = readRDS(file = 'results/fs/miRNA/fs_results.rds')
config_miRNA = fs_results_miRNA$config
res_miRNA    = fs_results_miRNA$res %>%
  mutate(n_features = sapply(selected_features, length))

## Performance plot ----
res_miRNA %>%
  ggplot(aes(x = method, y = score_rsmp, fill = learner_id)) +
  geom_boxplot() +
  labs(y = 'C-index (3 x 5-fold CV)', x = 'Wrapper-based FS method', title = 'miRNA') +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(res_path, '/perf.png'), width = 4.5, height = 5, dpi = 350)

## Number of features plot ----
res_miRNA %>%
  ggplot(aes(x = method, y = n_features, fill = learner_id)) +
  geom_boxplot() +
  labs(y = 'Number of selected features', x = 'Wrapper-based FS method', title = 'miRNA') +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(res_path, '/nfeat.png'), width = 4.5, height = 5, dpi = 350)

## Selection frequency plots ----
learners = unique(res_miRNA$learner_id)
data_list = list()
index = 1
for (wrapper_method in config_miRNA$wrapper_methods) {
  for (learner in learners) {
    selfeats_list = res_miRNA %>%
      filter(method == wrapper_method, learner_id == learner) %>%
      pull(selected_features)
    if (length(selfeats_list) > 0) {
      freq_tbl = get_consensus_features(selfeats_list)

      id = paste0(wrapper_method, '-', learner)
      p = feat_freq_barplot(freq_tbl, top_n = 10, title = paste0('miRNA (', id, ')'))
      ggsave(plot = p, filename = paste0(res_path, '/selfreq_', id, '.png'),
        width = 4.5, height = 5, dpi = 350)

      # calculate stability measures for each list of selected features
      data_list[[id]] = tibble(method = wrapper_method, learner_id = learner,
        jaccard = stabm::stabilityJaccard(features = selfeats_list, correction.for.chance = 'none'),
        nogueira = stabm::stabilityNogueira(features = selfeats_list, p = nfeats_miRNA))
    }
  }
}

stab_tbl = bind_rows(data_list)

## Stability assessment ----
stab_tbl %>%
  ggplot(aes(x = method, y = jaccard, fill = learner_id)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x = 'Wrapper-based FS method', y = 'Jaccard Similarity', title = 'miRNA') +
  ylim(c(0, 0.72)) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(res_path, '/jaccard_stab.png'), width = 4.5, height = 5, dpi = 350)

stab_tbl %>%
  ggplot(aes(x = method, y = nogueira, fill = learner_id)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x = 'Wrapper-based FS method', y = 'Nogueira Similarity', title = 'miRNA') +
  ylim(c(0, 0.72)) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(res_path, '/nogueira_stab.png'), width = 4.5, height = 5, dpi = 350)

## Consensus features and stability (RFE-RSFs) ----

#' `feat_ll` => list of list of features
feat_ll = lapply(lrns_list, function(lrns) {
  res_miRNA %>%
    filter(method == 'rfe', learner_id %in% lrns) %>%
    pull(selected_features)
})

data_list = list()
index = 1
for (index in 1:length(lrns_list)) {
  # learners combined
  print(lrns_list[[index]])

  feat_list = feat_ll[[index]]

  # Stability
  jaccard = stabm::stabilityJaccard(features = feat_list, correction.for.chance = 'none')
  nogueira = stabm::stabilityNogueira(features = feat_list, p = nfeats_miRNA)

  # Consensus features
  freq_tbl = get_consensus_features(feat_list)
  data_list[[index]] = list(
    id = paste0(str_remove(lrns_list[[index]], 'rsf_'), collapse = '-'),
    jaccard = jaccard,
    nogueira = nogueira,
    freq_tbl = freq_tbl)
}

rsfstab_tbl = lapply(data_list, function(l) { l$freq_tbl = NULL; l }) %>%
  dplyr::bind_rows()

rsfstab_tbl %>%
  mutate(id = forcats::fct_reorder(id, jaccard, .desc = TRUE)) %>%
  tidyr::pivot_longer(cols = c('jaccard', 'nogueira'), names_to = 'stab_measure') %>%
  ggplot(aes(x = id, y = value, fill = stab_measure)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x = 'RSF feature set combinations', y = 'Similarity score', title = 'miRNA') +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(res_path, '/rsfcmb_stab.png'), width = 6, height = 5, dpi = 350)

miRNA_cons_freq_tbl = data_list[[which(
  sapply(data_list, `[[`, 1) == 'cindex-logrank-maxstat'
)]]$freq_tbl
p = feat_freq_barplot(miRNA_cons_freq_tbl, top_n = 10,
  title = paste0('miRNA (RSF consensus)'))
ggsave(plot = p, filename = paste0(res_path, '/cons_selfreq.png'),
  width = 5, height = 5, dpi = 350)

# Methyl ----
#' see script `fs_meth.R`

## path to save plots
res_path = 'img/fs/Methyl'
if (!dir.exists(res_path)) {
  dir.create(res_path, recursive = TRUE)
}

## Methyl task
task_Methyl = readRDS(file = 'data/tasks.rds')$Methyl
nfeats_Methyl = length(task_Methyl$feature_names)

## Methyl eFS results
fs_results_Methyl = readRDS(file = 'results/fs/Methyl/fs_results.rds')
config_Methyl = fs_results_Methyl$config
res_Methyl    = fs_results_Methyl$res %>%
  mutate(n_features = sapply(selected_features, length))

## Performance plot ----
res_Methyl %>%
  ggplot(aes(x = method, y = score_rsmp, fill = learner_id)) +
  geom_boxplot() +
  labs(y = 'C-index (3 x 5-fold CV)', x = 'Wrapper-based FS method', title = 'Methyl') +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(res_path, '/perf.png'), width = 4.5, height = 5, dpi = 350)

## Number of features plot ----
res_Methyl %>%
  ggplot(aes(x = method, y = n_features, fill = learner_id)) +
  geom_boxplot() +
  labs(y = 'Number of selected features', x = 'Wrapper-based FS method', title = 'Methyl') +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(res_path, '/nfeat.png'), width = 4.5, height = 5, dpi = 350)

## Selection frequency plots ----
learners = unique(res_Methyl$learner_id)
data_list = list()
index = 1
for (wrapper_method in config_Methyl$wrapper_methods) {
  for (learner in learners) {
    selfeats_list = res_Methyl %>%
      filter(method == wrapper_method, learner_id == learner) %>%
      pull(selected_features)
    if (length(selfeats_list) > 0) {
      freq_tbl = get_consensus_features(selfeats_list)

      id = paste0(wrapper_method, '-', learner)
      p = feat_freq_barplot(freq_tbl, top_n = 10, title = paste0('Methyl (', id, ')'))
      ggsave(plot = p, filename = paste0(res_path, '/selfreq_', id, '.png'),
        width = 4.5, height = 5, dpi = 350)

      # calculate stability measures for each list of selected features
      data_list[[id]] = tibble(method = wrapper_method, learner_id = learner,
        jaccard = stabm::stabilityJaccard(features = selfeats_list, correction.for.chance = 'none'),
        nogueira = stabm::stabilityNogueira(features = selfeats_list, p = nfeats_Methyl))
    }
  }
}

stab_tbl = bind_rows(data_list)

## Stability assessment ----
stab_tbl %>%
  ggplot(aes(x = method, y = jaccard, fill = learner_id)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x = 'Wrapper-based FS method', y = 'Jaccard Similarity', title = 'Methyl') +
  ylim(c(0, 0.72)) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(res_path, '/jaccard_stab.png'), width = 4.5, height = 5, dpi = 350)

stab_tbl %>%
  ggplot(aes(x = method, y = nogueira, fill = learner_id)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x = 'Wrapper-based FS method', y = 'Nogueira Similarity', title = 'Methyl') +
  ylim(c(0, 0.72)) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(res_path, '/nogueira_stab.png'), width = 4.5, height = 5, dpi = 350)

## Consensus features and stability (RFE-RSFs) ----

#' `feat_ll` => list of list of features
feat_ll = lapply(lrns_list, function(lrns) {
  res_Methyl %>%
    filter(method == 'rfe', learner_id %in% lrns) %>%
    pull(selected_features)
})

data_list = list()
index = 1
for (index in 1:length(lrns_list)) {
  # learners combined
  print(lrns_list[[index]])

  feat_list = feat_ll[[index]]

  # Stability
  jaccard = stabm::stabilityJaccard(features = feat_list, correction.for.chance = 'none')
  nogueira = stabm::stabilityNogueira(features = feat_list, p = nfeats_Methyl)

  # Consensus features
  freq_tbl = get_consensus_features(feat_list)
  data_list[[index]] = list(
    id = paste0(str_remove(lrns_list[[index]], 'rsf_'), collapse = '-'),
    jaccard = jaccard,
    nogueira = nogueira,
    freq_tbl = freq_tbl)
}

rsfstab_tbl = lapply(data_list, function(l) { l$freq_tbl = NULL; l }) %>%
  dplyr::bind_rows()

rsfstab_tbl %>%
  mutate(id = forcats::fct_reorder(id, jaccard, .desc = TRUE)) %>%
  tidyr::pivot_longer(cols = c('jaccard', 'nogueira'), names_to = 'stab_measure') %>%
  ggplot(aes(x = id, y = value, fill = stab_measure)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x = 'RSF feature set combinations', y = 'Similarity score', title = 'Methyl') +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(res_path, '/rsfcmb_stab.png'), width = 6, height = 5, dpi = 350)

Methyl_cons_freq_tbl = data_list[[which(
  sapply(data_list, `[[`, 1) == 'cindex-logrank-maxstat'
)]]$freq_tbl
p = feat_freq_barplot(Methyl_cons_freq_tbl, top_n = 10,
  title = paste0('Methyl (RSF consensus)'))
ggsave(plot = p, filename = paste0(res_path, '/cons_selfreq.png'),
  width = 4.5, height = 5, dpi = 350)

# Save FS frequency tables ----
fs_consensus_res = list(mRNA = mRNA_cons_freq_tbl, miRNA = miRNA_cons_freq_tbl,
  CNA = CNA_cons_freq_tbl, Methyl = Methyl_cons_freq_tbl)
saveRDS(fs_freq, file = 'results/fs/fs_consensus_res.rds')

#' Examples how to get consensus features per data modality

#' All selected mRNA features listed in descending frequency order
mRNA_cons_feats = mRNA_cons_freq_tbl$feat_name
head(mRNA_cons_feats)
#' most frequently selected features - top 5%
mRNA_cons_feats[1:(0.05*length(mRNA_cons_feats))]
#' most frequently selected features - top 20
mRNA_cons_feats[1:20]
#' most frequently selected features - above `freq_cutoff`
freq_cutoff = 0.5

mRNA_cons_freq_tbl %>% filter(freq > freq_cutoff) %>% pull(feat_name)
miRNA_cons_freq_tbl %>% filter(freq > freq_cutoff) %>% pull(feat_name)
Methyl_cons_freq_tbl %>% filter(freq > freq_cutoff) %>% pull(feat_name)
CNA_cons_freq_tbl %>% filter(freq > freq_cutoff) %>% pull(feat_name)
