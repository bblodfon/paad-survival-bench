#' Due to recent mlr3 update (14.0.1), a task's `row_roles` must have
#' a `test` slot instead of a `early_stopping`, otherwise you can't train
#' on them!
library(mlr3verse)

# Original tasks ----
tasks = readRDS(file = 'data/tasks.rds')

tasks = lapply(tasks, function(task) {
  task$row_roles = list(
    use = task$row_roles$use,
    holdout = task$row_roles$holdout,
    test = integer(0))
  task
})

# keep only the ginle modality tasks and a dual-modality that I had used
tasks$`miRNA-Methyl` = NULL
tasks$`mRNA-miRNA` = NULL
tasks$`mRNA-Methyl` = NULL
tasks$`Methyl-CNA` = NULL
tasks$`mRNA-miRNA-Methyl` = NULL
tasks$`mRNA-miRNA-CNA` = NULL
tasks$`mRNA-Methyl-CNA` = NULL
tasks$`miRNA-Methyl-CNA` = NULL
tasks$`mRNA-miRNA-Methyl-CNA` = NULL

saveRDS(tasks, file = 'data/tasks.rds')

# mRNA filtered task ----
task_mRNA = readRDS(file = 'data/task_mRNA_flt.rds')
task_mRNA$row_roles = list(
  use = task_mRNA$row_roles$use,
  holdout = task_mRNA$row_roles$holdout,
  test = integer(0))
saveRDS(task_mRNA, file = 'data/task_mRNA_flt.rds')

# Clinical task ----
task_clinical = readRDS(file = 'data/task_clinical.rds')
task_clinical$row_roles = list(
  use = task_clinical$row_roles$use,
  holdout = task_clinical$row_roles$holdout,
  test = integer(0))
saveRDS(task_clinical, file = 'data/task_clinical.rds')

# FS-tasks ----
tasks_fs = readRDS(file = 'data/tasks_fs.rds')

tasks_fs = lapply(tasks_fs, function(task) {
  task$row_roles = list(
    use = task$row_roles$use,
    holdout = task$row_roles$holdout,
    test = integer(0))
  task
})

saveRDS(tasks_fs, file = 'data/tasks_fs.rds')

