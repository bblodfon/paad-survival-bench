# Bootstrap functions to use with an mlr3 test task and get
# Confidence Intervals for a performance metric (e.g. C-index)

#' `data` is a data.table/data.frame object with the test data
#' and has to have the same structure (features and target columns)
#' as the task that was used to train the provided learner
boot_fun = function(data, index, learner, measure) {
  learner$predict_newdata(data[index])$score(measure)
}

boot_ci = function(data, learner, measure = mlr3::msr('surv.cindex'),
  num_threads = 1, nrsmps = 1000, parallel = 'multicore', include_boot_res = TRUE) {

  boot_res = boot::boot(data, statistic = boot_fun, R = nrsmps,
    parallel = parallel, ncpus = num_threads, learner = learner, measure = measure)

  bootci_res = boot::boot.ci(boot_res, type = c('basic', 'norm', 'perc', 'bca'))

  if (include_boot_res)
    return(list(boot_res = boot_res, bootci_res = bootci_res))
  else
    return(list(bootci_res = bootci_res))
}
