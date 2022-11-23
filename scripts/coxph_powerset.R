#'##################################################
#' A simple benchmark where the CoxPH model is
#' trained and tested on all possible combinations of
#' PAAD tasks (task powerset), after FS has been
#' applied (so a lot less features in total per modality)
#'
#' - Simple train/test split + bootstrap C-index
#' on test set (based on mRNA data split benchmarks)
#'##################################################
library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
library(corrplot)
library(psych)
source('scripts/helpers.R')

res_path = 'results/powerset_bench/'
if (!dir.exists(res_path)) {
  dir.create(res_path)
}
img_path = 'img/powerset_bench/'
if (!dir.exists(img_path)) {
  dir.create(img_path)
}

# Reproducibility
set.seed(42)

# Logging (less)
lgr::get_logger('bbotk')$set_threshold('warn')
lgr::get_logger('mlr3')$set_threshold('warn')

# Measures of performance
harrell_c = msr('surv.cindex')
harrell_c$label = 'HarrellC'
uno_c = msr('surv.cindex', weight_meth = 'G2')
uno_c$label = 'UnoC'
ibrier = msr('surv.graf')
ibrier$label = 'IBrier'

# Tasks (see `prepare_tasks_after_fs.R`)
tasks = readRDS(file = 'data/tasks_fs.rds')

# ~70%/30% train/test split
data_split = readRDS(file = 'data/data_split.rds')
train_indx = data_split$train_indx
test_indx  = data_split$test_indx

# CoxPH benchmark
learner = lrn('surv.coxph', id = 'CoxPH', fallback = lrn('surv.kaplan'))

res_list = list()
index = 1
for (task in tasks) {
  message('Train start: (', learner$id, ', ', task$id, ')')
  trained_learner = learner$clone(deep = TRUE)$train(task, train_indx)

  message('Calculate bootstrap CIs (Harrell\'s C-index)')
  test_boot_harc = get_boot_ci(task = task, train_indx = train_indx,
    test_indx = test_indx, learner = trained_learner, measure = harrell_c)
  message('Calculate bootstrap CIs (Uno\'s C-index)')
  test_boot_unoc = get_boot_ci(task = task, train_indx = train_indx,
    test_indx = test_indx, learner = trained_learner, measure = uno_c)
  message('Calculate bootstrap CIs (IBrier Score)')
  test_boot_ibrier = get_boot_ci(task = task, train_indx = train_indx,
    test_indx = test_indx, learner = trained_learner, measure = ibrier)

  res_list[[index]] = tibble::tibble(
    learner_id = learner$id,
    task_id = task$id,
    harc = test_boot_harc$t,
    unoc = test_boot_unoc$t,
    ibrier = test_boot_ibrier$t,
    harc_median = test_boot_harc$t_median,
    unoc_median = test_boot_unoc$t_median,
    ibrier_median = test_boot_ibrier$t_median
  )
  index = index + 1
}

coxph_res = dplyr::bind_rows(res_list)

# save results
saveRDS(coxph_res, file = paste0(res_path, 'coxph_res.rds'))

coxph_res = readRDS(file = paste0(res_path, 'coxph_res.rds'))

# Plot performance box-plots
coxph_res %>%
  mutate(task_id = forcats::fct_reorder(task_id, harc_median, .desc = TRUE)) %>%
  ggplot(aes(x = task_id, y = harc, fill = task_id)) +
  geom_boxplot(show.legend = FALSE) +
  geom_hline(yintercept = 0.5, linetype = 'dotted', color = 'red') +
  labs(x = 'Tasks', y = 'Harrell\'s C-index', title = 'Bootstrap Test Set (CoxPH)') +
  ylim(c(0,1)) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(filename = paste0(img_path, 'coxph_harc.png'), width = 8, height = 5, dpi = 300)

coxph_res %>%
  mutate(task_id = forcats::fct_reorder(task_id, unoc_median, .desc = TRUE)) %>%
  ggplot(aes(x = task_id, y = unoc, fill = task_id)) +
  geom_boxplot(show.legend = FALSE) +
  geom_hline(yintercept = 0.5, linetype = 'dotted', color = 'red') +
  labs(x = 'Tasks', y = 'Uno\'s C-index', title = 'Bootstrap Test Set (CoxPH)') +
  ylim(c(0,1)) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(filename = paste0(img_path, 'coxph_unoc.png'), width = 8, height = 5, dpi = 300)

coxph_res %>%
  mutate(task_id = forcats::fct_reorder(task_id, ibrier_median)) %>%
  ggplot(aes(x = task_id, y = ibrier, fill = task_id)) +
  geom_boxplot(show.legend = FALSE) +
  geom_hline(yintercept = 0.25, linetype = 'dotted', color = 'red') +
  labs(x = 'Tasks', y = 'Integrated Brier Score', title = 'Bootstrap Test Set (CoxPH)') +
  #ylim(c(0,1)) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(filename = paste0(img_path, 'coxph_ibrier.png'), width = 8, height = 5, dpi = 300)

# Correlation between different performance scores
harc_median_tbl = coxph_res %>% select(task_id, harc_median) %>% unique() #%>% arrange(desc(harc_median))
unoc_median_tbl = coxph_res %>% select(task_id, unoc_median) %>% unique()
ibrier_median_tbl = coxph_res %>% select(task_id, ibrier_median) %>% unique()
stopifnot(harc_median_tbl$task_id == unoc_median_tbl$task_id,
          unoc_median_tbl$task_id == ibrier_median_tbl$task_id)
med_tbl = tibble(harc = harc_median_tbl$harc_median, unoc = unoc_median_tbl$unoc_median,
  ibrier = ibrier_median_tbl$ibrier_median)

## get the p-values
corr_test = psych::corr.test(x = med_tbl, method = 'kendall')

if (FALSE) {
  png(filename = paste0(img_path, 'coxph_cor.png'), width = 5, height = 5, units = 'in', res = 300)
  corrplot::corrplot(cor(med_tbl, method = 'kendall'), type = 'upper',
    p.mat = corr_test$p, insig = 'pch') # label-sig
  dev.off()
}
