library(mlr3verse)

# mRNA task ----
tasks = readRDS(file = 'data/tasks.rds')
task_mRNA = tasks$mRNA
task_mRNA$ncol

## Keep 1000 most variant mRNA features
po_varflt = po('filter', filter = flt('variance'), filter.nfeat = 1000)

task_mRNA_flt = po_varflt$train(list(task_mRNA))[[1]]
task_mRNA_flt$ncol

task_mRNA = task_mRNA_flt