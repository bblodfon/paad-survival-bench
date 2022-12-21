#' Plot benchmarking results from `bench_powerset.R`
library(mlr3verse)
library(mlr3benchmark)
library(tidyverse)
source('scripts/helpers.R')
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

# Load and tidy data ----
#' Get results file
res_path = 'results/powerset_bench/'
#' using `max_nfeats` = 13 in `prepare_tasks_after_fs.R`
res_file = paste0(res_path, 'bench_res_22112022_234158.rds')
#' using `max_nfeats` = 26 in `prepare_tasks_after_fs.R`
#' too large file to upload in GitHub
# res_file = paste0(res_path, 'bench_res_more_features_26112022_173300.rds')
res = readRDS(file = res_file)
config = res$config
bm_res = res$bench_res

#' Get benchmark name
bm_name = sub(pattern = '.rds', x = basename(res_file), replacement = '')
img_path = paste0('img/powerset_bench/', bm_name)
if (!dir.exists(img_path)) {
  dir.create(img_path)
}

if (!dir.exists(paste0(img_path, '/tuning_plots'))) {
  dir.create(paste0(img_path, '/tuning_plots'))
}

#' Get all bootstrap test results in a nice tidy format
all_boot_scores = list()
for (row_id in 1:nrow(bm_res)) {
  results = bm_res$results[[row_id]]
  task_id = bm_res$task_id
  lrn_name = nice_lrn_label(results$trained_learner$id)

  tbl = tibble(
    task_id = bm_res$task_id[row_id],
    lrn_id  = bm_res$lrn_id[row_id]
  )

  # for each test measure, get all test bootstrap scores + median value
  for (msr in results$test_boot) {
    msr_label = msr$msr_label
    msr_label_median = paste0(msr_label, '_median')
    tbl = dplyr::bind_cols(tbl, tibble(!!(msr$msr_label) := msr$t,
      !!(msr_label_median) := msr$t_median))
  }

  all_boot_scores[[row_id]] = tbl
}
all_boot_res = dplyr::bind_rows(all_boot_scores)
all_boot_res

# Tuning-validation plots ----
if (FALSE) { # too many plots, produce some manually
  for (row_id in 1:nrow(bm_res)) {
    results = bm_res$results[[row_id]]
    task_id = bm_res$task_id[row_id]
    lrn_name = nice_lrn_label(results$trained_learner$id)

    get_hpc_perf_plot(hpc_res = results$hpc_res,
      measure_name = nice_msr_label(results$train_measure),
      title = paste0(lrn_name, ' - ', task_id))

    ggsave(filename = paste0(img_path, '/tuning_plots/', results$trained_learner$id, '_',
      task_id, '.png'), width = 6, height = 5, dpi = 300)
  }
}

# Plot Learner performance across tasks ----
for (learner_id in unique(all_boot_res$lrn_id)) {
  message('Learner: ', learner_id)

  # Harrell's C-index
  all_boot_res %>%
    filter(lrn_id == learner_id) %>%
    mutate(task_id = forcats::fct_reorder(task_id, HarrellC_median, .desc = TRUE)) %>%
    ggplot(aes(x = task_id, y = HarrellC, fill = task_id)) +
    geom_boxplot(show.legend = FALSE) +
    geom_hline(yintercept = 0.5, linetype = 'dotted', color = 'red') +
    labs(x = 'Tasks', y = nice_msr_label('Harrell'),
      title = paste0('Bootstrap Test Set (', nice_lrn_label(learner_id), ')')) +
    ylim(c(0,1)) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  ggsave(filename = paste0(img_path, '/', learner_id, '_harc.png'),
    width = 8, height = 5, dpi = 300)

  # Uno's C-index
  all_boot_res %>%
    filter(lrn_id == learner_id) %>%
    mutate(task_id = forcats::fct_reorder(task_id, UnoC_median, .desc = TRUE)) %>%
    ggplot(aes(x = task_id, y = UnoC, fill = task_id)) +
    geom_boxplot(show.legend = FALSE) +
    geom_hline(yintercept = 0.5, linetype = 'dotted', color = 'red') +
    labs(x = 'Tasks', y = nice_msr_label('Uno'),
      title = paste0('Bootstrap Test Set (', nice_lrn_label(learner_id), ')')) +
    ylim(c(0,1)) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  ggsave(filename = paste0(img_path, '/', learner_id, '_unoc.png'),
    width = 8, height = 5, dpi = 300)
}

# Learner ranking ----
#' get CoxPH results (see `coxph_powerset.R`)
coxph_res = readRDS(file = paste0(res_path, 'coxph_res.rds'))

#' homogenize CoxPH results to `all_boot_res` tibble
coxph_res = coxph_res %>%
  select(task_id, learner_id, harc, harc_median, unoc, unoc_median) %>%
  rename(lrn_id = learner_id, HarrellC = harc, HarrellC_median = harc_median,
    UnoC = unoc, UnoC_median = unoc_median)

#' and merge them
all_boot_res = dplyr::bind_rows(all_boot_res, coxph_res)

#' aggregate by taking only the median scores
aggr_res = all_boot_res %>%
  select(task_id, lrn_id, HarrellC_median, UnoC_median) %>%
  distinct() %>%
  rename(learner_id = lrn_id) %>%
  mutate(task_id = factor(task_id), learner_id = factor(learner_id))

#' check performance
aggr_res %>% arrange(desc(HarrellC_median))
aggr_res %>% arrange(desc(UnoC_median))

ba = BenchmarkAggr$new(aggr_res)

## Global Friedman's test ----
ba$friedman_test()
# Learners significantly different using Harrell's C, not for Uno's!

#' Learners aggregated performance across tasks (boxplots)
for (measure in ba$measures) {
  aggr_res %>%
    ggplot(aes(x = learner_id, y = .data[[measure]], fill = learner_id)) +
    geom_boxplot(show.legend = FALSE) +
    theme_bw(base_size = 14) +
    labs(title = 'Aggregated Performance across tasks',
      y = paste0(nice_msr_label(measure), ' (median)')) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(filename = paste0(img_path, '/aggr_perf_', measure, '.png'),
    width = 7, height = 5, dpi = 300)
}

## Nemenyi and Critical Difference plots ----
meas = 'HarrellC_median' # only Harrell's C-index for the rest of the analysis
ranks = get_ranks(ba, measure = meas, minimize = FALSE)
rmat_lrn = ranks$rmat_lrn
sort(rowMeans(rmat_lrn)) # from best to worst

cd = get_cd(ba, rmat_lrn, measure = meas)
cd$test # post-hoc Friedman-Nemenyi pairwise tests

fn_plot = plot_fn(cd)
ggsave(plot = fn_plot, filename = paste0(img_path, '/fn_', meas, '.png'),
  width = 5, height = 5, dpi = 300)

cd_plot = plot_cd(cd)
ggsave(plot = cd_plot, filename = paste0(img_path, '/cd_', meas, '.png'),
  width = 7, height = 5, dpi = 300)

# Clinical + CoxPH vs best learner (and task)
best_lrn = names(sort(rowMeans(rmat_lrn)))[1]
best_lrn

#' get best and worst performing task for best learner
task_ids = aggr_res %>%
  filter(learner_id == best_lrn) %>%
  arrange(desc(!!sym(meas))) %>% # from best to worse (C-index)
  pull(task_id) %>%
  as.character()

best_task = task_ids[1]
best_task
worst_task = task_ids[length(task_ids)]
worst_task

harc_best = all_boot_res %>%
  filter(lrn_id == best_lrn, task_id == best_task) %>%
  pull(HarrellC)
harc_worse = all_boot_res %>%
  filter(lrn_id == best_lrn, task_id == worst_task) %>%
  pull(HarrellC)

# significant difference between the two "end" tasks!
wilcox.test(x = harc_best, y = harc_worse)

# Performance comparison with baseline model (CoxPH-Clinical) ----
harc_coxph_clinical = all_boot_res %>%
  filter(lrn_id == 'CoxPH', task_id == 'Clinical') %>%
  pull(HarrellC)
wilcox.test(x = harc_coxph_clinical, y = harc_best) # significant, but does it matter?

# Make a boxplot to show
harc_tbl = tibble::tibble(
  harc_coxph_clinical = harc_coxph_clinical,
  #harc_worse = harc_worse,
  harc_best = harc_best
)
tbl = harc_tbl %>%
  tidyr::pivot_longer(cols = everything(), names_to = 'name', values_to = 'vals') %>%
  mutate(name = case_when(
    #name == 'harc_worse' ~ paste0(best_lrn, '-', worst_task),
    name == 'harc_best' ~ paste0(best_lrn, '-', best_task),
    name == 'harc_coxph_clinical' ~ 'CoxPH-Clinical'
  )
)

comp = list(
  #c('RSF-worst-task', paste0(best_lrn, '-', worst_task)),
  c('CoxPH-Clinical', paste0(best_lrn, '-', best_task)))
stat_test = tbl %>%
  rstatix::wilcox_test(formula = vals ~ name, comparisons = comp) %>%
  rstatix::add_significance('p') %>%
  rstatix::add_y_position()

# boxplot
tbl %>%
  ggplot(aes(x = name, y = vals, color = name)) +
  geom_boxplot(show.legend = FALSE) +
  ggpubr::stat_pvalue_manual(stat_test, label = "p = {p} ({p.signif})") +
  labs(x = 'Learner+Task combo', y = 'Harrell\'s C-index') +
  theme_bw(base_size = 14) #+
  #theme(axis.text.x = element_text(angle = 15, hjust = 1))
ggsave(filename = paste0(img_path, '/best_learner_vs_coxph_', meas, '.png'),
  width = 5, height = 5, dpi = 300)

# Correlation between Harrell's and Uno's C-index results per learner ----
l = aggr_res %>%
  group_by(learner_id) %>%
  group_split()

cor_list = lapply(l, function(tbl) {
  harc = tbl$HarrellC_median
  unoc = tbl$UnoC_median

  tibble(
    learner_id = unique(tbl$learner_id),
    kend_cor = cor(harc, unoc, method = 'kendall'),
    pval = psych::corr.test(harc, unoc, method = 'kendall')$p
  )
})

cor_tbl = dplyr::bind_rows(cor_list)
cor_tbl = cor_tbl %>% # add mean_rank
  mutate(mean_rank = cd$data[levels(cor_tbl$learner_id),]$mean_rank)
cor_tbl %>% arrange(kend_cor)

# Task ranking ----
res = task_ranking(ranks)
#' normalized rank sum score and ks statistic correlate
cor(c(res$nrs_mat), c(res$ks_mat))
cor(c(res$nrs_pval), c(res$ks_pval))

# define colors for C-index
col_fun = colorRamp2(c(0, 0.5, 1), c("red", "white", "green"))
png(filename = paste0(img_path, '/task_ranking_rs.png'), width = 7, height = 5,
  units = 'in', res = 450)
Heatmap(matrix = res$nrs_mat, name = 'RS score',
  row_title = 'Learners', column_title = 'Single omics', col = col_fun,
  cell_fun = function(j, i, x, y, width, height, fill) {
    pval = res$nrs_pval[i, j]
    if(pval < 0.0001) {
      grid.text('****', x, y)
    } else if(pval < 0.001) {
      grid.text('***', x, y)
    } else if(pval < 0.01) {
      grid.text('**', x, y)
    } else if(pval < 0.05) {
      grid.text('*', x, y)
    }
  }
)
dev.off()

png(filename = paste0(img_path, '/task_ranking_ks.png'), width = 7, height = 5,
  units = 'in', res = 450)
Heatmap(matrix = res$ks_mat, name = 'KS score',
  row_title = 'Learners', column_title = 'Single omics', col = col_fun,
  cell_fun = function(j, i, x, y, width, height, fill) {
    pval = res$ks_pval[i, j]
    if(pval < 0.0001) {
      grid.text('****', x, y)
    } else if(pval < 0.001) {
      grid.text('***', x, y)
    } else if(pval < 0.01) {
      grid.text('**', x, y)
    } else if(pval < 0.05) {
      grid.text('*', x, y)
    }
  }
)
dev.off()
