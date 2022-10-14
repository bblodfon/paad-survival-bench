library(purrr)
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
# Better do: grid[mlr3misc::map_chr(grid$learner, `[[`, 'id') %in% 'CoxNet.tuned']
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

#' @param `learner` for the `AutoFSelector`
#' @param `task` for the training of the `AutoFSelector`
#' @param `resampling` for the `AutoFSelector`. Default: 5x 5-fold CV
#' @param `measure` for the `AutoFSelector`. Default: Harrell's C-index
#' @param `method` which FS wrapper method to use. Either 'rfe' (recursive feature
#' elimination) or 'ga' (genetic algorithm search). Default: 'rfe'
#' @param `repeats` how many times should we run the FS method
#' @param `rfe_n_features` number of features that signals the termination of the RFE. Default: 2
#' @param `rfe_feature_fraction` % of features to keep in each iteration of the RFE. Default: 0.85
#' @param `ga_iters` how many feature subsets to search in GA
#' @param `ga_zeroToOneRatio` control how sparse are the feature subsets in GA
#' (i.e. less features selected). As a rule of thumb, nfeatures-of-task/zeroToOneRatio
#' will be approx. equal to the number of features in each subset chosen by GA
run_wrapper_fs = function(learner, task, resampling =
    rsmp('repeated_cv', folds = 5, repeats = 5),
  measure = msr('surv.cindex'), method = 'rfe', repeats = 100,
  rfe_n_features = 2, rfe_feature_fraction = 0.85,
  ga_iters = 100, ga_zeroToOneRatio = 125
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
