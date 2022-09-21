# PAAD clinical features as a Survival mlr3 task (mlr3proba::TaskSurv)
library(mlr3verse)
library(mlr3proba)
library(tidyverse)

# Read clinical data ----
clinical_var_mat = readRDS(file = 'data/clinical_var_mat.rds')

# Convert to tibble, filter and preprocess features ----
clinical_var_tbl = clinical_var_mat %>%
  as_tibble(rownames = 'patient_id') %>%
  rename(status = vital_status) %>%
  readr::type_convert(guess_integer = TRUE) %>%
  mutate(time = ifelse(is.na(days_to_death), days_to_last_followup, days_to_death)) %>%
  rename(age = years_to_birth) %>%
  select(-c('patient_id', 'days_to_death', 'days_to_last_followup',
    'tumor_tissue_site', 'histological_type', # same value for all patients
    'date_of_initial_pathologic_diagnosis',
    'number_pack_years_smoked', 'year_of_tobacco_smoking_onset', # many NAs, and
    # it may be not a risk factor for pancreatic cancer
    'pathology_M_stage', 'pathology_N_stage', 'pathology_T_stage', # redundant,
    # since we have the `pathologic_stage` which is a combination assessment of
    # the 3 features above (T,N,M stages)
    'race', 'ethnicity' # most are white, non-hispanic/latinos
  )) %>%
  filter(time > 0) %>% # remove one patient with time = 0 (see also `prepare_tasks.R`)
  mutate(across(where(is.character), forcats::as_factor)) %>% # convert chars to factors
  # fix reference level for `pathologic_stage`
  #mutate(pathologic_stage = relevel(x = pathologic_stage, ref = 'stage ia'))
  # `pathologic_stage` is an ordinal categorical variable so encode it with integers
  mutate(pathologic_stage = as.integer(case_when(
    pathologic_stage == 'stage ia' ~ 1,
    pathologic_stage == 'stage iia' ~ 2,
    pathologic_stage == 'stage iib' ~ 3,
    pathologic_stage == 'stage iii' ~ 4,
    pathologic_stage == 'stage ib' ~ 5,
    pathologic_stage == 'stage iv' ~ 6
  )))

# Convert to Survival task ----
paad_task = mlr3proba::as_task_surv(
  x = clinical_var_tbl,
  time = 'time', event = 'status', id = 'Clinical'
)
paad_task
paad_task$col_info # check the baseline levels if they are okay!
# `method=treatment` in `encode` leaves out the first (baseline) factor level

paad_task$missings()

# Impute + encode factors ----
impute_fct = po('imputelearner', lrn('classif.rpart'), id = 'impute_fct')
impute_num = po('imputelearner', lrn('regr.rpart'), id = 'impute_num')
fix_fct = po('fixfactors') # remove strange xyz.factor columns that stay after imputation
encode_fct = po('encode', method = 'treatment') # n-1 columns for n factor levels

pre = impute_fct %>>% fix_fct %>>% impute_num %>>% encode_fct

task_clinical = pre$train(paad_task)[[1L]]
task_clinical
all(task_clinical$missings() == 0)

## check that we have same order in the `Surv` target of the mRNA task
tasks = readRDS(file = 'data/tasks.rds')
task_mRNA = tasks$mRNA
all(task_mRNA$truth()[,1] == task_clinical$truth()[,1]) # time
all(task_mRNA$truth()[,2] == task_clinical$truth()[,2]) # status

# Save task ----
saveRDS(task_clinical, file = 'data/task_clinical.rds')
