library(mlr3verse)
library(mlr3proba)
suppressMessages(library(dplyr))
library(readr)
library(combinat)

# Get survival data
clinical_var_mat = readRDS(file = 'data/clinical_var_mat.rds')
surv_tbl = clinical_var_mat %>%
  as_tibble(rownames = 'patient_id') %>%
  rename(status = vital_status) %>%
  readr::type_convert() %>%
  select(patient_id, status, starts_with('days_to')) %>%
  mutate(time = ifelse(is.na(days_to_death), days_to_last_followup, days_to_death))

# Remove patient with time = 0 (causes problems with some models)
surv_tbl %>% arrange(time)
surv_tbl_flt = surv_tbl %>% filter(time > 0) # 145 patients left

# Load preprocessed TCGA data
mrna_mat   = readRDS(file = 'data/mrna_mat.rds')
mirna_mat  = readRDS(file = 'data/mirna_mat.rds')
methyl_mat = readRDS(file = 'data/methyl_mat.rds')
cna_mat    = readRDS(file = 'data/cna_snp_mat.rds')

# filter to patients with sufficient survival information
mrna_mat   = mrna_mat[surv_tbl_flt %>% pull(patient_id),]
mirna_mat  = mirna_mat[surv_tbl_flt %>% pull(patient_id),]
methyl_mat = methyl_mat[surv_tbl_flt %>% pull(patient_id),]
cna_mat    = cna_mat[surv_tbl_flt %>% pull(patient_id),]

# check patient ID order
stopifnot((surv_tbl_flt %>% pull(patient_id)) == rownames(mrna_mat))
stopifnot((surv_tbl_flt %>% pull(patient_id)) == rownames(mirna_mat))
stopifnot((surv_tbl_flt %>% pull(patient_id)) == rownames(methyl_mat))
stopifnot((surv_tbl_flt %>% pull(patient_id)) == rownames(cna_mat))

# make syntactic proper names otherwise task creation fails (dashes => dots)
colnames(mrna_mat)  = make.names(colnames(mrna_mat))
colnames(mirna_mat) = make.names(colnames(mirna_mat))
colnames(cna_mat)   = make.names(colnames(cna_mat))

# add data type prefix to each matrix's colnames
# (making colnames unique across data types)
colnames(mrna_mat)   = paste0('mRNA.', colnames(mrna_mat))
colnames(mirna_mat)  = paste0('miRNA.', colnames(mirna_mat))
colnames(cna_mat)    = paste0('cna.', colnames(cna_mat))
colnames(methyl_mat) = paste0('met.', colnames(methyl_mat))

# Create tasks
mrna_task = as_task_surv(
  x = cbind(surv_tbl_flt %>% select(status, time), mrna_mat),
  time = 'time', event = 'status', id = 'mRNA'
)
mirna_task = as_task_surv(
  x = cbind(surv_tbl_flt %>% select(status, time), mirna_mat),
  time = 'time', event = 'status', id = 'miRNA'
)
cna_task = as_task_surv(
  x = cbind(surv_tbl_flt %>% select(status, time), cna_mat),
  time = 'time', event = 'status', id = 'CNA'
)
methyl_task = as_task_surv(
  x = cbind(surv_tbl_flt %>% select(status, time), methyl_mat),
  time = 'time', event = 'status', id = 'Methyl'
)

# Create combined tasks (all subsets)
po_featureunion = po("featureunion")
task_set = c(mrna_task, mirna_task, methyl_task, cna_task)

task_subsets = lapply(1:length(task_set), combinat::combn, x = task_set,
  simplify = FALSE) %>% unlist(recursive = FALSE)

tasks = list()
for (task_subset in task_subsets) {
  task = po_featureunion$train(task_subset)$output
  taskid = paste0(sapply(task_subset, function(l) l$id), collapse = '-')
  task$id = taskid
  tasks[[taskid]] = task
}

stopifnot(length(tasks) == 15) # all subset tasks with 4 data types (excluding empty set)

# Save list of tasks
saveRDS(tasks, file = 'data/tasks.rds')
