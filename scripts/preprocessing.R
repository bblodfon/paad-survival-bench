# Preprocessing of mRNA, miRNA, CNA and Methylation data
# Make sure you have enough free memory when running this script!
library(curatedTCGAData)
library(TCGAutils)
library(dplyr)
library(readr)

source(file = 'scripts/helpers.R')

# read data
cancer_data = readRDS(file = 'data/cancer_data.rds')

set.seed(42)
#'###########################
# Preprocess mRNA matrix ----
#'###########################
mrna_mat = t(assay(cancer_data[,,"PAAD_RNASeq2GeneNorm-20160128"]))
rownames(mrna_mat) = cancer_data$patientID

print(paste0('#mRNA features: ', dim(mrna_mat)[2])) # 20501 features

# remove features with more than 20% of NAs
mrna_mat = remove_NAs(mrna_mat) # no NAs found

# remove features with more than 20% of zeros
mrna_mat = remove_zeros(mrna_mat)

# Keep top 65% most variant features (reduce dimensionality)
mrna_mat = flt_var(mrna_mat, percentage = 0.65)

# remove near-zero variance and highly correlated features
trans1 = caret::preProcess(mrna_mat, method = c("nzv", "corr"),
  cutoff = 0.9, freqCut = 95/5, uniqueCut = 20, verbose = T)
trans1
mrna_mat = predict(trans1, mrna_mat)

# log2(x+1) transform
mrna_mat = log2_trans(mrna_mat)

# center and scale
center_scale_trans = caret::preProcess(mrna_mat,
  method = c("center", "scale"), verbose = T)
mrna_mat = predict(center_scale_trans, mrna_mat)

print(paste0('#mRNA features: ', dim(mrna_mat)[2])) # 10465 features

# save to file
saveRDS(object = mrna_mat, file = 'data/mrna_mat.rds')

#'############################
# Preprocess miRNA matrix ----
#'############################
mirna_mat = t(assay(cancer_data[,,"PAAD_miRNASeqGene-20160128"]))
rownames(mirna_mat) = cancer_data$patientID

print(paste0('#miRNA features: ', dim(mirna_mat)[2])) # 1046 features

# remove features with more than 20% of NAs
mirna_mat = remove_NAs(mirna_mat) # no NAs found

# remove features with more than 20% of zeros
mirna_mat = remove_zeros(mirna_mat)

# remove near-zero variance and highly correlated features
trans2 = caret::preProcess(mirna_mat, method = c("nzv", "corr"),
  cutoff = 0.9, freqCut = 95/5, uniqueCut = 20, verbose = T)
trans2
mirna_mat = predict(trans2, mirna_mat)

# log2(x+1) transform
mirna_mat = log2_trans(mirna_mat)

# center and scale
center_scale_trans = caret::preProcess(mirna_mat,
  method = c("center", "scale"), verbose = T)
mirna_mat = predict(center_scale_trans, mirna_mat)

print(paste0('#miRNA features: ', dim(mirna_mat)[2])) # 339 features

# save to file
saveRDS(object = mirna_mat, file = 'data/mirna_mat.rds')

#'##########################
# Preprocess CNA matrix ----
#'##########################
# CNA data are converted from `RaggedExperiment` to `RangedSummarizedExperiment`
# Copy number regions become genes using a weighted mean function
cancer_data_simplified = TCGAutils::simplifyTCGA(cancer_data)

# https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/#copy-number-variation-analysis-pipeline
# The GDC further transforms these copy number values into segment mean values, which are equal to log2(copy-number/2)
# Diploid regions will have a segment mean of zero, amplified regions will have positive values, and deletions will have negative values.
cna_snp_mat = t(assay(cancer_data_simplified[,,"PAAD_CNASNP-20160128_simplified"]))
rownames(cna_snp_mat) = cancer_data$patientID
print(paste0('#CNA features: ', dim(cna_snp_mat)[2])) # 22911 features

# remove features with at least 1 NA value
cna_snp_mat = remove_NAs(cna_snp_mat, cutoff = 0) # 1250 features removed

# remove near-zero variance and highly correlated features
# choose cutoff = 0.99 to remove mostly identical features
# center and scale (no need to log-transform as the values are already in log scale)
trans3 = caret::preProcess(cna_snp_mat,
  method = c("nzv", "corr", "center", "scale"),
  cutoff = 0.99, freqCut = 95/5, uniqueCut = 20, verbose = T)
trans3
cna_snp_mat = predict(trans3, cna_snp_mat)

print(paste0('#CNA features: ', dim(cna_snp_mat)[2])) # 2831 features

# save to file
saveRDS(object = cna_snp_mat, file = 'data/cna_snp_mat.rds')

#'#######################################
# Preprocess Methylation data matrix ----
#'#######################################
# Beta values between 0 (low levels of methylation) and 1 (high levels)
methyl_data = cancer_data[,,"PAAD_Methylation-20160128"]
methyl_mat = t(assay(methyl_data))
rownames(methyl_mat) = cancer_data$patientID
print(paste0('#Methylation features: ', dim(methyl_mat)[2])) # 485577 features

# Get CpGs associated with sex chromosomes (X,Y)
methyl_metadata = rowData(experiments(methyl_data)[[1]]) %>%
  as_tibble(rownames = 'CpGs')
methyl_metadata %>% count(Chromosome) %>% print(n = 25)
somatic_cpgs = methyl_metadata %>%
  filter(!is.na(Gene_Symbol), !Chromosome %in% c('X','Y',NA)) %>%
  pull(CpGs)

# check
all(somatic_cpgs %in% colnames(methyl_mat))

# subset CpGs
methyl_mat = methyl_mat[, somatic_cpgs]
methyl_mat = methyl_mat %>% as.matrix()

# remove features with at least 1 NA
methyl_mat = remove_NAs(methyl_mat, cutoff = 0) # ~64K features removed

# Keep top 5% most variant features (reduce dimensionality)
methyl_mat = flt_var(methyl_mat, percentage = 0.05)

# remove near-zero variance and highly correlated features
trans4 = caret::preProcess(methyl_mat,
  method = c("nzv", "corr"),
  cutoff = 0.9, freqCut = 95/5, uniqueCut = 20, verbose = T)
trans4
methyl_mat = predict(trans4, methyl_mat)

# log2(x+1) transform
methyl_mat = log2_trans(methyl_mat)

# center and scale
center_scale_trans = caret::preProcess(methyl_mat,
  method = c("center", "scale"), verbose = T)
methyl_mat = predict(center_scale_trans, methyl_mat)

print(paste0('#Methylation features: ', dim(methyl_mat)[2])) # 11805 features

# save to file
saveRDS(object = methyl_mat, file = 'data/methyl_mat.rds')
