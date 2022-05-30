library(mlr3verse)
library(mlr3proba)
library(dplyr)
library(readr)

# Make Tasks ----
# get survival data
clinical_var_mat = readRDS(file = 'data/clinical_var_mat.rds')
surv_tbl = clinical_var_mat %>% as_tibble(rownames = 'patient_id') %>%
  rename(status = vital_status) %>%
  readr::type_convert() %>%
  select(patient_id, status, starts_with('days_to')) %>%
  mutate(time = ifelse(is.na(days_to_death), days_to_last_followup, days_to_death))

# Remove patient with time = 0 (causes problems with some models)
surv_tbl %>% arrange(time)
surv_tbl_flt = surv_tbl %>% filter(time > 0) # 145 patients left

# load preprocessed TCGA data
mrna_mat  = readRDS(file = 'data/mrna_mat.rds')
mirna_mat = readRDS(file = 'data/mirna_mat.rds')

# filter to patients with sufficient survival information
mrna_mat  = mrna_mat[surv_tbl_flt %>% pull(patient_id),]
mirna_mat = mirna_mat[surv_tbl_flt %>% pull(patient_id),]

# check patient ID order
stopifnot((surv_tbl_flt %>% pull(patient_id)) == rownames(mrna_mat))
stopifnot((surv_tbl_flt %>% pull(patient_id)) == rownames(mirna_mat))

# make syntactic proper names otherwise task creation fails (dashes => dots)
colnames(mrna_mat)  = make.names(colnames(mrna_mat))
colnames(mirna_mat) = make.names(colnames(mirna_mat))

# Create tasks
mrna_task = as_task_surv(
  x = cbind(surv_tbl_flt %>% select(status, time), mrna_mat),
  time = 'time', event = 'status', id = 'mRNA'
)
mrna_task$col_info$label = 'mRNA'

mirna_task = as_task_surv(
  x = cbind(surv_tbl_flt %>% select(status, time), mirna_mat),
  time = 'time', event = 'status', id = 'miRNA'
)
mirna_task$col_info$label = 'miRNA'

# Make list of tasks
tasks = list(mRNA = mrna_task, miRNA = mirna_task)

# Save list
saveRDS(tasks, file = 'data/tasks.rds')
