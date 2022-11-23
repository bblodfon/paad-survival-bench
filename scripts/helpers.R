library(mlr3verse)
library(purrr)
library(ggplot2)
suppressMessages(library(dplyr))

# remove features that have more than cutoff% NAs
remove_NAs = function(mat, cutoff = 0.2) {
  percent_NAs = colSums(is.na(mat))/nrow(mat)
  features_num = sum(percent_NAs > cutoff)
  message(paste0('Remove ', features_num, ' features with more than ',
    cutoff * 100, '% NAs'))

  # remove features
  if (features_num > 0) {
    mat = mat[, -which(percent_NAs > cutoff)]
  }

  mat
}

# Remove features that have more than cutoff% 0's
# Use only for read count expression matrices (e.g. mRNA, miRNA)
remove_zeros = function(mat, cutoff = 0.2) {
  percent_zeros = colSums(mat == 0)/nrow(mat)
  features_num = sum(percent_zeros > cutoff)
  print(paste0('Remove ', features_num, ' features with more than ',
    cutoff * 100, '% zeros'))

  # remove features
  if (features_num > 0) {
    mat = mat[, -which(percent_zeros > cutoff)]
  }

  mat
}

# Log2 transform matrix
log2_trans = function(mat, offset_value = 1) {
  log(mat + offset_value, base = 2)
}

# sample your data
sample_df = function(df, n = 10) {
  df[sample(1:nrow(df), n), sample(1:ncol(df), n)]
}

# feature selection based on variance
# keep top %percentage features
# by default, var doesn't remove NAs
# better use with a matrix with no NAs
flt_var = function(mat, percentage = 0.5) {
  variances = apply(mat, 2, var)
  features_num = round(percentage * length(variances)) # how many features to keep?
  features_to_keep = names(sort(variances, decreasing = TRUE)[1:features_num])

  # subset matrix
  print(paste0('Keeping top ', percentage*100, '% of features with higher variance'))
  mat[,colnames(mat) %in% features_to_keep]
}

# filter a data.table from `benchmark_grid()` to specific learners (manual version)
# Better do: grid[mlr3misc::map_chr(grid$learner, `[[`, 'id') %in% learner$id]
flt_learners = function(grid, learner_ids) {
  learner_list = grid$learner
  rows_to_keep = list()
  j = 1
  for (i in 1:length(learner_list)) {
    learner = learner_list[[i]]
    if (learner$id %in% learner_ids) {
      rows_to_keep[[j]] = i
      j = j + 1
    }
  }
  grid[purrr::flatten_dbl(rows_to_keep)]
}

#' Does the measure require the extra parameters `task` and `train_set`?
#'
#' @param measure an mlr3 `Measure` object. Currently C-index and Integrated
#' Brier Score are supported.
extra_params_required = function(measure) {
  use_extra_params = FALSE
  # Uno's C-index
  if (measure$id == 'surv.cindex' && measure$param_set$values$weight_meth == 'G2')
    use_extra_params = TRUE

  # Integrated Brier Score
  if (measure$id == 'surv.graf')
    use_extra_params = TRUE

  use_extra_params
}

#' Calculate C-indexes for all hyperparameter configurations tested in an
#' AutoTuner for both the train and the test set of the Task used to train
#' the AutoTuner
#' @param at AutoTuner object, already trained, with non-NULL `at$archive`
#' @param task Task used to train and test AutoTuner's learner
#' @param train_indx row_ids from `task` used for training the AutoTuner `at`
get_cindex_all_hps = function(at, task, train_indx) {
  hpc_list = at$archive$data$x_domain

  # empty base learner (delete trained model)
  base_lrn = at$learner$clone(deep = TRUE)
  base_lrn$reset()

  # test set
  test_indx = setdiff(seq_len(task$nrow), train_indx)

  res = list()
  for (i in 1:length(hpc_list)) {
    hpc = hpc_list[[i]]
    base_lrn$param_set$values = mlr3misc::insert_named(base_lrn$param_set$values, hpc)

    base_lrn$train(task, train_indx)
    res[[i]] = dplyr::tibble(
      index = i,
      train_cindex = base_lrn$predict(task, train_indx)$score(),
      test_cindex  = base_lrn$predict(task, test_indx) $score()
    )

    base_lrn$reset()
  }

  hpc_res = dplyr::bind_rows(res)

  # add the stored average (CV) C-index on the train data
  hpc_res = hpc_res %>%
    mutate(traincv_cindex = at$archive$data$surv.cindex)
}

#' Calculate a performance score for all hyperparameter configurations used in
#' an `AutoTuner`. This function calculates the performance score on both the
#' train and the test set of the `Task` used to train the `AutoTuner`.
#' The aim is to check for overfitting in the train set and how the performance
#' on the test set changes during selection of different hyperparameter configurations.
#'
#' @param at AutoTuner object, already trained, with non-NULL `at$archive`
#' @param task Task used to train and test `AutoTuner`'s learner
#' @param train_indx row_ids from `task` used for training the AutoTuner `at`
#' @param test_indx row_ids from `task` not used for training the AutoTuner `at`
#' @param measure an mlr3 `Measure` object. Used to caulcate the performance of
#' a particular hpc configuration on the train and test set. Currently it can be
#' Harrell's C-index, Uno's C-index or the Integrated Brier Score.
#'
#' @return a tibble with performance scores
#'
#' @note You have to keep track of which type of score was used during the
#' training of the AutoTuner (`rsmp_score`) and the calculation of performance
#' on the train and test sets (`train_score`, `test_score`)!
get_scores_hps = function(at, task, train_indx, test_indx, measure) {
  # some checks
  mlr3::assert_task(task)
  mlr3::assert_measure(measure)

  # do we need {task, train_indx} in the `$score()`?
  use_extra_params = extra_params_required(measure)

  # get hyperparameter config list
  hpc_list = at$archive$data$x_domain

  # empty base learner (delete trained model)
  base_lrn = at$learner$clone(deep = TRUE)

  res = list()
  for (i in 1:length(hpc_list)) {
    hpc = hpc_list[[i]]

    base_lrn$reset()
    base_lrn$param_set$values = mlr3misc::insert_named(base_lrn$param_set$values, hpc)
    base_lrn$train(task, train_indx)

    if (!use_extra_params) {
      train_score = base_lrn$predict(task, train_indx)$score(measure)
      test_score  = base_lrn$predict(task, test_indx)$score(measure)
    } else {
      train_score = base_lrn$predict(task, train_indx)$
        score(measure, task = task, train_set = train_indx)
      test_score  = base_lrn$predict(task, test_indx)$
        score(measure, task = task, train_set = train_indx)
    }

    res[[i]] = dplyr::tibble(index = i, train_score = train_score,
      test_score = test_score)
  }

  hpc_res = dplyr::bind_rows(res)

  # add the stored resampling score on the train data
  # This assumes the `measure$id` used when training the AutoTuner starts with 'surv'
  measure_name = grep(x = colnames(at$archive$data), pattern = '^surv',
    value = TRUE)
  hpc_res = hpc_res %>%
    mutate(rsmp_score = at$archive$data[[measure_name]], .after = 1)

  hpc_res
}

#' Plot the results of `get_scores_hps`
get_hpc_perf_plot = function(hpc_res, rand_perf_score = 0.5,
  xlab = 'Number of evaluations/hpcs (BO)', measure_name = 'C-index') {
  hpc_res %>% ggplot2::ggplot(aes(x = index)) +
    geom_line(aes(y = train_score, color = 'Train'), linewidth = 0.3) +
    geom_line(aes(y = rsmp_score, color = 'Train (rsmp)'), linewidth = 0.3) +
    geom_line(aes(y = test_score, color = 'Test'), linewidth = 0.7) +
    scale_color_manual(name = '',
      values = c(
        # colors from RColorBrewer::brewer.pal(3, 'Set1')
        'Train' = '#E41A1C',
        'Train (rsmp)' = '#377EB8',
        'Test' = '#4DAF4A')
    ) +
    geom_hline(yintercept = rand_perf_score, linetype = 'dotted', color = 'red') +
    labs(x = xlab, y = measure_name) +
    theme_bw(base_size = 14) + theme(legend.position = 'top')
}

#' @param `learner` for the `AutoFSelector`
#' @param `task` for the training of the `AutoFSelector`
#' @param `resampling` for the `AutoFSelector`. Default: 3x 5-fold CV
#' @param `measure` for the `AutoFSelector`. Default: Harrell's C-index
#' @param `method` which FS wrapper method to use. Either 'rfe' (recursive feature
#' elimination) or 'ga' (genetic algorithm search). Default: 'rfe'
#' @param `repeats` how many times should we run the FS method
#' @param `rfe_n_features` number of features that signals the termination of the RFE. Default: 2
#' @param `rfe_feature_fraction` % of features to keep in each iteration of the RFE. Default: 0.85
#' @param `ga_iters` how many feature subsets to search in GA (equal to number of generations)
#' @param `ga_zeroToOneRatio` control how sparse are the feature subsets in GA
#' (i.e. less features selected). As a rule of thumb, `task$ncol`/`zeroToOneRatio`
#' will be approx. equal to the number of features in each subset chosen by GA
#' @param `ga_popSize` initial population size (number of feature subsets considered)
run_wrapper_fs = function(learner, task,
  resampling = rsmp('repeated_cv', folds = 5, repeats = 3),
  measure = msr('surv.cindex'), method = 'rfe', repeats = 100,
  rfe_n_features = 2, rfe_feature_fraction = 0.85,
  ga_iters = 100, ga_zeroToOneRatio = 200, ga_popSize = 1000
) {

  if (!method %in% c('rfe', 'ga')) {
    stop('Only RFE and GA wrapper methods allowed')
  }

  if (method == 'rfe') { # Recursive Feature Elimination
    if (!'importance' %in% learner$properties) {
      message('Learner ' , learner$id, ' doesn\'t provide importance scores')
      return()
    }

    at = mlr3fselect::AutoFSelector$new(
      learner = learner,
      resampling = resampling,
      measure = measure,
      terminator = trm('none'),
      fselector = fs('rfe', n_features = rfe_n_features,
        feature_fraction = rfe_feature_fraction),
      store_models = TRUE
    )
  } else if (method == 'ga') { # Genetic Algorithm Search
    at = mlr3fselect::AutoFSelector$new(
      learner = learner,
      resampling = resampling,
      measure = measure,
      terminator = trm('evals', n_evals = ga_iters),
      fselector = fs('genetic_search', zeroToOneRatio = ga_zeroToOneRatio)
    )
  }

  # Run Wrapper Algorithm
  res = list()
  for (i in 1:repeats) {
    message(i, '/', repeats, ': ', method, ', ', learner$id)

    at$train(task)

    sel_features = at$fselect_instance$result_feature_set
    score_rsmp   = at$archive$best()[[measure$id]]

    res[[i]] = tibble::tibble(method = method, learner_id = learner$id,
      selected_features = list(sel_features), score_rsmp = score_rsmp)
  }

  dplyr::bind_rows(res)
}

#' @param `sel_feat_list` a list containing vectors of features
#' @return a `tibble` with the features in descending selection frequency order
get_consensus_features = function(selfeats_list) {
  repeats = length(selfeats_list)
  res = sort(table(unlist(selfeats_list)), decreasing = TRUE)
  times = as.vector(unname(res))
  tibble(feat_name = names(res), times = times, freq = times/repeats)
}

#' plots the result (`freq_tbl`) of `get_consensus_features()`
feat_freq_barplot = function(freq_tbl, top_n = 10, title = '') {
  freq_tbl %>%
    slice(1:top_n) %>%
    mutate(feat_name = forcats::fct_reorder(feat_name, times, .desc = FALSE)) %>%
    ggplot(aes(x = feat_name, y = freq)) +
    geom_bar(stat = "identity", fill = '#377EB8', show.legend = FALSE) +
    ggpubr::theme_classic2(base_size = 14) +
    labs(x = 'Feature name', y = 'Selection Frequency', title = title) +
    scale_y_continuous(labels = scales::label_percent()) +
    coord_flip()
}

#' @param `task_list` named list of mlr3 `Task`s
#' @return a list of `Task`s with every possible combination of features from `task_list`
get_task_powerset = function(task_list) {
  po_featureunion = mlr3pipelines::po('featureunion')

  task_subsets = lapply(1:length(task_list), combinat::combn, x = task_list,
                        simplify = FALSE) %>% unlist(recursive = FALSE)

  task_powerset = list()
  for (task_subset in task_subsets) {
    task = po_featureunion$train(task_subset)[[1L]]
    task_id = paste0(sapply(task_subset, function(l) l$id), collapse = '-')
    task$id = task_id
    task_powerset[[task_id]] = task
  }

  # excluding the empty subset
  stopifnot(length(task_powerset) == (2^length(task_list)-1))

  task_powerset
}

#' Convenient function to check if the input is an mlr3 measure +
#' has an assigned label
#'
#' @param measure mlr3 `Measure`
check_measure = function(measure) {
  mlr3::assert_measure(measure)

  # stop if there is no label!
  if (is.na(measure$label)) {
    stop('Please add a label to measure ', measure$id)
  }
}

#' Bootstrap function to use with an mlr3 test task and get
#' Confidence Intervals for a performance metric (e.g. C-index)
#' @param task an mlr3 `Task`
#' @param train_indx row ids of the training set of the given `task`
#' @param test_indx row ids of the test set of the given `task`
#' @param learner a mlr3 `Learner`, trained on the given `task` using the `train_indx`
#' rows
#' @param measure an mlr3 `Measure` object. Currently C-index (Harrell's and
#' Uno's) and the Integrated Brier Score are supported.
#' @param nthreads how many threads to use? Passed on to the `boot` function's
#' `ncpus` argument. Default is to use all available cores via `parallelly::availableCores()`
#' @param nrsmps The number of bootstrap replicates/resamplings
#' @param parallel Type of parallel operation to be used in `boot` function
#' @param full_result (FALSE) By default we return a list with only what we deem
#' necessary, i.e. the confidence intervals of `boot::boot.ci()` function, the
#' observed score on the test set (`t0`), the scores on the bootstrap resamplings
#' of the test set (`t`) and if there were any warnings or unstable intervals
#' calculated.
#' If `TRUE`, a list is returned, having the outputs of both the `boot::boot()`
#' and the `boot::boot.ci()` function.
get_boot_ci = function(task, train_indx, test_indx, learner, measure = measure,
  nthreads = parallelly::availableCores(), nrsmps = 1000, full_result = FALSE) {

  # some checks
  mlr3::assert_task(task)
  mlr3::assert_learner(learner)
  check_measure(measure)

  # get the test dataset
  data = task$data(rows = test_indx)

  boot_res = boot::boot(
    data, statistic = boot_fun, R = nrsmps, parallel = 'multicore', ncpus = nthreads,
    learner = learner, measure = measure, task = task, train_indx = train_indx
  )

  bootci_res = suppressWarnings(
    boot::boot.ci(boot_res, type = c('basic', 'norm', 'perc', 'bca'))
  )

  if (full_result)
    return(list(boot_res = boot_res, bootci_res = bootci_res))
  else {
    # any warnings/errors? (hacked version, see `boot:::print.bootci`)
    output = utils::capture.output(bootci_res)
    warn = grep(pattern = 'Warning', x = output, value = TRUE)
    unstable = grep(pattern = 'unstable', x = output, value = TRUE)

    return(
      list(
        msr_label = measure$label, msr_id = measure$id,
        t0 = bootci_res$t0, t = boot_res$t[,1], t_mean = mean(boot_res$t, na.rm = T),
        t_median = median(boot_res$t, na.rm = T), R = bootci_res$R,
        normal = bootci_res$normal, basic = bootci_res$basic,
        percent = bootci_res$percent, bca = bootci_res$bca,
        warn = warn, unstable = unstable
      )
    )
  }
}

#' @param `data` data.table/data.frame object with the test data
#' that has the same structure (features and target columns)
#' as the `task` that was used to train the provided `learner`
#' on the `train_indx` rows
#' @param `learner` trained Learner object
#' @param `measure` an mlr3 `Measure`
#' @param task an mlr3 `Task`
#' @param train_indx row ids of the training set of the given `task`
boot_fun = function(data, index, learner, measure, task, train_indx) {
  # do we need {task, train_indx} in the `$score()`?
  use_extra_params = extra_params_required(measure)

  if (!use_extra_params)
    learner$predict_newdata(data[index])$score(measure)
  else
    learner$predict_newdata(data[index])$score(measure, task = task, train_set = train_indx)
}

#' Get a list of Survival Learners
#'
#' @param nthreads for implicit parallelization (RSFs, XGBoost)
get_surv_lrns = function(nthreads = parallelly::availableCores()) {
  # CoxPH
  coxph = lrn('surv.coxph', id = 'CoxPH', fallback = lrn('surv.kaplan'))

  # CoxNet
  coxnet = lrn('surv.glmnet', id = 'CoxNet', fallback = lrn('surv.kaplan'),
    standardize = FALSE, maxit = 10^3)
  coxlasso = lrn('surv.glmnet', id = 'CoxLasso', fallback = lrn('surv.kaplan'),
    standardize = FALSE, maxit = 10^3)

  # Survival Tree
  surv_tree = lrn('surv.rpart', id = 'SurvivalTree', fallback = lrn('surv.kaplan'))

  ## Random Survival Forests
  rsf_cindex  = lrn('surv.ranger', verbose = FALSE, id = 'SurvivalForestCIndex',
    fallback = lrn('surv.kaplan'),
    num.threads = nthreads,
    splitrule = 'C') # C-index
  rsf_logrank = lrn('surv.ranger', verbose = FALSE, id = 'SurvivalForestLogRank',
    fallback = lrn('surv.kaplan'),
    num.threads = nthreads,
    splitrule = 'logrank')
  rsf_maxstat = lrn('surv.ranger', verbose = FALSE, id = 'SurvivalForestMaxStat',
    fallback = lrn('surv.kaplan'),
    num.threads = nthreads,
    splitrule = 'maxstat')

  # CoxBoost
  coxboost = lrn('surv.coxboost', id = 'CoxBoost',
    fallback = lrn('surv.kaplan'),
    standardize = FALSE, # data already standardized
    return.score = FALSE) # don't need this in the output

  # XGBoost
  xgboost_cox = lrn('surv.xgboost', nthread = nthreads, booster = 'gbtree',
    fallback = lrn('surv.kaplan'), objective = 'survival:cox', id = 'XGBoostCox')
  xgboost_aft = lrn('surv.xgboost', nthread = nthreads, booster = 'gbtree',
    fallback = lrn('surv.kaplan'), objective = 'survival:aft', id = 'XGBoostAFT')

  list(coxph = coxph, coxnet = coxnet, coxlasso = coxlasso, surv_tree = surv_tree,
    rsf_cindex = rsf_cindex, rsf_logrank = rsf_logrank, rsf_maxstat = rsf_maxstat,
    coxboost = coxboost, xgboost_cox = xgboost_cox, xgboost_aft = xgboost_aft)
}

#' Add a `distr` prediction for survival learners that don't support it
#' or overwrite the one that some learners have.
#'
#' @param lrn_list a list of survival mlr3 learners (`LearnerSurv`)
#' @param lrn_ids a vector of learner ids, for which the `distr` will be forcibly
#' added or overwritten
#'
#' @note Deciding on the form of the composed survival prediction: if the
#' learner's id has 'Cox' somewhere, then we decide on a Proportional Hazards
#' (`ph`) form. If the learner's id has 'AFT' somewhere, then we decide on an
#' Accelerated Failure Time (`aft`) form. Otherwise, we go with a PH!
add_distr_pred = function(lrn_list, lrn_ids = c('SurvivalTree', 'CoxBoost')) {
  lapply(lrn_list, function(learner) {
    # if learner doesn't have a `distr` or is in the `lrn_list`, it gets one!
    if ((!'distr' %in% learner$predict_types) || (learner$id %in% lrn_ids)) {
      learner = mlr3pipelines::ppl('distrcompositor',
        learner = learner,
        estimator = 'kaplan',
        form = ifelse(grepl(pattern = 'AFT', x = learner$id), 'aft', 'ph'),
        overwrite = TRUE,
        graph_learner = FALSE
      ) %>% mlr3pipelines::GraphLearner$new(id = learner$id)
    }

    learner
  })
}

#' @return a list with suggested parameter spaces for some survival learners
#' (to be used for tuning)
get_param_spaces = function() {
  coxnet = paradox::ps(
    lambda = p_dbl(1e-03, 10, logscale = TRUE),
    alpha  = p_dbl(0, 1)) # from Ridge to Lasso penalty

  coxlasso = paradox::ps(
    lambda = p_dbl(1e-03, 10, logscale = TRUE)
  )

  surv_tree = paradox::ps(
    minsplit = p_int(1, 64, logscale = TRUE),
    cp = p_dbl(1e-04, 1, logscale = TRUE))

  rsf = paradox::ps(
    num.trees = p_int(100, 1500),
    mtry.ratio = p_dbl(0.1, 0.9),
    min.node.size = p_int(1, 20))

  coxboost = paradox::ps(
    stepno = p_int(50, 500),
    # leave at default => 9 * sum(status == 1)?
    penalty = p_int(10, 1000, logscale = TRUE),
    # leave at default => 1?
    stepsize.factor = p_dbl(1e-01, 10, logscale = TRUE),
    # use the penalized scores - `hpscore` is faster to calculate
    criterion = p_fct(c('pscore', 'hpscore')))

  xgboost_cox = paradox::ps(
    nrounds = p_int(150, 1500),
    eta = p_dbl(1e-04, 1, logscale = TRUE),
    max_depth = p_int(2, 8),
    min_child_weight = p_dbl(1, 128, logscale = TRUE))

   xgboost_cox_reg = paradox::ps(
     nrounds = p_int(150, 1500),
     eta = p_dbl(1e-04, 1, logscale = TRUE),
     max_depth = p_int(2, 8),
     min_child_weight = p_dbl(1, 128, logscale = TRUE),
     # more possibility for regularization
     alpha  = p_dbl(1e-03, 10, logscale = TRUE),
     lambda = p_dbl(1e-03, 10, logscale = TRUE),
     gamma  = p_dbl(0, 5),
     # early stopping => 10% of nrounds
     .extra_trafo = function(x, param_set) {
       x$early_stopping_rounds = as.integer(ceiling(0.1 * x$nrounds))
       x
     }
     # depending on the version of xgboost learner in mlr3extralearners,
     # you may have to set `early_stopping_set = 'train' or 'test'` for
     # this parameter space to work
   )

   xgboost_aft = xgboost_cox$
     clone(deep = TRUE)$
     add(paradox::ps(
       aft_loss_distribution = p_fct(c('normal', 'logistic', 'extreme')),
       aft_loss_distribution_scale = p_dbl(0.5, 2.0)
       )
     )

   xgboost_aft_reg = xgboost_cox_reg$
     clone(deep = TRUE)$
     add(paradox::ps(
       aft_loss_distribution = p_fct(c('normal', 'logistic', 'extreme')),
       aft_loss_distribution_scale = p_dbl(0.5, 2.0)
       )
     )

  list(coxnet = coxnet, coxlasso = coxlasso, surv_tree = surv_tree, rsf = rsf,
    coxboost = coxboost, xgboost_cox = xgboost_cox, xgboost_cox_reg = xgboost_cox_reg,
    xgboost_aft = xgboost_aft, xgboost_aft_reg = xgboost_aft_reg)
}

#' @param nthreads for implicit parallelization (RSFs, XGBoost)
#' @param ids learners ids to filter (these are manually set)
#'
#' @return a `data.table` object with learner ids, `Learner` objects and their
#' corresponding `ParamSet`s to be used for tuning
get_lrns_and_ps = function(nthreads, ids) {
  learners = suppressWarnings(get_surv_lrns(nthreads))
  pss = get_param_spaces()

  lrn_tbl = tibble::tribble(
    ~id, ~learner, ~param_set,
    'coxnet', learners$coxnet, pss$coxnet,
    'coxlasso', learners$coxlasso, pss$coxlasso,
    'surv_tree', learners$surv_tree, pss$surv_tree,
    'rsf_cindex', learners$rsf_cindex, pss$rsf,
    'rsf_logrank', learners$rsf_logrank, pss$rsf,
    'rsf_maxstat', learners$rsf_maxstat, pss$rsf,
    'coxboost', learners$coxboost, pss$coxboost,
    'xgboost_cox', learners$xgboost_cox, pss$xgboost_cox,
    # need to set `early_stopping_set` appropriately
    'xgboost_cox_reg', learners$xgboost_cox, pss$xgboost_cox_reg,
    'xgboost_aft', learners$xgboost_aft, pss$xgboost_aft,
    # need to set `early_stopping_set` appropriately
    'xgboost_aft_reg', learners$xgboost_aft, pss$xgboost_aft_reg,
  )

  # assert that ids are unique
  stopifnot(length(unique(lrn_tbl$id)) == length(lrn_tbl$id))

  if (!missing(ids)) {
    lrn_tbl = lrn_tbl %>% filter(id %in% ids)
  }

  data.table::as.data.table(lrn_tbl)
}

#' Run `AutoTuner` on train set of a given task and get bootstrap estimates of
#' performance on given test set.
#'
#' Tune a `learner` with a given hyperparameter `search_space` on the train set
#' (`train_indx`) of a task, using the given `resampling` scheme and performance
#' `measure`. By default the tuner uses Bayesian Optimization and stops after
#' `nevals` evaluations. The hyperparameter configuration with the best average
#' performance across the resamplings is chosen and a `trained_learner` using
#' that configuration is returned.
#'
#' With `all_hpcs_perf` enabled (default), you get a `tibble` with the
#' performance scores (`measure`) of all hyperparameter configurations
#' tested (`nevals`) on the train set, the test set and the resampling estimate
#' (e.g. average CV performance).
#'
#' With `boot_test` enabled (default), you get the performance estimates
#' on bootstrapped test datasets (`nrsmps` = 1000 by default) using `nthreads`
#' to speed up computations. If `nthreads` is missing, we use all available cores.
#' Parallelization is not currently working with `boot()` so we suggest 1 thread.
#' `test_measures` is a list of `mlr3` measures. Currently supported Harrell's
#' C-index, Uno's C-index and the Integrated Brier score.
#'
#' @note `test_measures` list can be different from `measure` - e.g. one could be
#' Harrell's C-index and the other the Integrated Brier Score. Note that not
#' all learners have predictions that can be used by any survival metric, so be
#' careful which learner + measures you use!
#' If `test_measures` is not specified at all, it becomes the same as `measure`.
#'
run_at = function(learner, task, train_indx, test_indx, resampling, measure,
  nevals, search_space, all_hpcs_perf = TRUE, boot_test = TRUE, test_measures,
  nrsmps = 1000, nthreads = nthreads) {

  if (missing(test_measures))
    test_measures = list(measure)

  if (missing(nthreads))
    nthreads = parallelly::availableCores()

  # some checks
  mlr3::assert_task(task)
  mlr3::assert_learner(learner)
  mlr3::assert_resampling(resampling)
  check_measure(measure)
  for (test_measure in test_measures) {
    check_measure(test_measure)
  }

  # Create AutoTuner
  at = mlr3tuning::AutoTuner$new(
    learner = learner, # initial object isn't changed during training
    resampling = resampling, # initial object isn't changed during training
    measure = measure,
    search_space = search_space,
    terminator = trm('evals', n_evals = nevals),
    tuner = tnr('mbo') # careful, need to define it here
  )

  # Train AutoTuner
  message('Train start: (', learner$id, ', ', task$id, ')')
  at$train(task, train_indx)
  message('Train finished: (', learner$id, ', ', task$id, ')')

  # Create trained learner using the best hyperparameter config (hpc) from BO
  # Surrogate doesn't assign the best hpc by default, see
  # https://github.com/mlr-org/mlr3mbo/issues/84#issuecomment-1230278238
  best_hpc = at$archive$best()$x_domain[[1L]]
  trained_learner = learner$clone(deep = TRUE)
  trained_learner$param_set$values = mlr3misc::insert_named(
    trained_learner$param_set$values, best_hpc)
  trained_learner$train(task, train_indx)

  if (all_hpcs_perf) {
    # Performance scores for all hpcs (rsmp, train, test)
    message('Calculate ', measure$id, ' (', measure$label, ')',
      ' for all hpcs (test + train set)')
    hpc_res = get_scores_hps(at = at, task = task, train_indx = train_indx,
      test_indx = test_indx, measure = measure)
  }

  if (boot_test) {
    # get bootstrap test set results of learner with best hpc
    # using all specified measures
    test_boot = list()
    for (test_measure in test_measures) {
      message('Calculate bootstrap CIs for ', test_measure$id, ' (',
        test_measure$label, ')')
      test_boot[[test_measure$label]] = get_boot_ci(task = task,
        train_indx = train_indx, test_indx = test_indx, learner = trained_learner,
        measure = test_measure, nrsmps = nrsmps, nthreads = nthreads)
    }
  }

  if (all_hpcs_perf && boot_test) {
    list(
      trained_learner = trained_learner,
      train_measure = measure$label,
      hpc_res  = hpc_res,
      test_boot = test_boot
    )
  } else if (all_hpcs_perf && !boot_test) {
    list(
      trained_learner = trained_learner,
      train_measure = measure$label,
      hpc_res  = hpc_res
    )
  } else if (!all_hpcs_perf && boot_test) {
    list(
      trained_learner = trained_learner,
      train_measure = measure$label,
      test_boot = test_boot
    )
  } else {
    list(
      trained_learner = trained_learner,
      train_measure = measure$label
    )
  }
}
