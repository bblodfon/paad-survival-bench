library(purrr)
library(dplyr)

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
      train_cindex = base_lrn$predict(task_mRNA, train_indx)$score(),
      test_cindex  = base_lrn$predict(task_mRNA, test_indx) $score()
    )

    base_lrn$reset()
  }

  hpc_res = dplyr::bind_rows(res)

  # add the stored average (CV) C-index on the train data
  hpc_res = hpc_res %>%
    mutate(traincv_cindex = at$archive$data$surv.cindex)
}
