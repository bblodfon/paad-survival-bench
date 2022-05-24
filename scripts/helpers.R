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
