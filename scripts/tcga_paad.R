# With this script we download the dataset and find a suitable assay subset
# to use for our ML task
library(curatedTCGAData)
library(TCGAutils)
library(dplyr)
library(readr)
library(UpSetR)

# Download Data ----

# check available disease codes
data('diseaseCodes', package = 'TCGAutils')
head(diseaseCodes)

disease_code = 'PAAD'

# see all assays (no download => `dry.run = TRUE`)
curatedTCGAData(diseaseCode = disease_code, version = '2.0.1', dry.run = TRUE)

# Download selected assays
my_assays = c('*RPPA*', '*Methylation*', '*miRNA*', '*RNASeq2GeneNorm*', '*Mutation*', '*CNA*')
cancer_data = curatedTCGAData(diseaseCode = disease_code, assays = my_assays, version = '2.0.1', dry.run = FALSE)

# Some checks
nrow(colData(cancer_data)) # number of patients (185)
head(cancer_data$patientID) # patientIDs

# Data filtering ----

# check sample types
data('sampleTypes', package = 'TCGAutils')
sampleTypes # Primary Solid Tumor = '01', Solid Tissue Normal = '11'

# what sample types do we have in the dataset/study across different data types?
TCGAutils::sampleTables(cancer_data)

# keep only primary tumor samples
cancer_data = TCGAutils::TCGAprimaryTumors(cancer_data)
TCGAutils::sampleTables(cancer_data) # verify: (only '01')

# keep only patients with the major `histological_type` (PDAC)
tbl_hist = table(cancer_data$histological_type, useNA = 'ifany')
tbl_hist # 154
major_hist_type = names(tbl_hist[1])
major_hist_type

logical_vec = cancer_data$histological_type == major_hist_type
cancer_data = cancer_data[, logical_vec & !is.na(logical_vec), ]

# see common subsets of samples
MultiAssayExperiment::upsetSamples(cancer_data)

# Decide which assays to keep (criterion: as many common samples as possible)
# (`upsetSamples()` arg `nsets` doesn't work)
patient_membership_df = do.call(
  function(...) { data.frame(..., check.names = TRUE) },
  lapply(mapToList(sampleMap(cancer_data)),
    function(minimap) {
      rownames(colData(cancer_data)) %in% minimap[["primary"]] * 1L
    }
  )
)
rownames(patient_membership_df) = rownames(colData(cancer_data)) # patient IDs
assay_names = colnames(patient_membership_df)
# change '.' back to '-'
assay_names = stringr::str_replace(string = assay_names, pattern = '\\.', replacement = '-')
colnames(patient_membership_df) = stringr::str_match(colnames(patient_membership_df),
  pattern = 'PAAD_([^.]+)')[,2] # compress assay names

# `nsets` keeps the n top largest sets
UpSetR::upset(data = patient_membership_df, nsets = 3, order.by = 'freq') # 146
UpSetR::upset(data = patient_membership_df, nsets = 4, order.by = 'freq') # 146
UpSetR::upset(data = patient_membership_df, nsets = 5, order.by = 'freq') # 132
UpSetR::upset(data = patient_membership_df, nsets = 6, order.by = 'freq') # 89

# choose specific sets
UpSetR::upset(data = patient_membership_df, sets = c('RNASeq2GeneNorm', 'miRNASeqGene', 'RPPAArray')) # 101
UpSetR::upset(data = patient_membership_df, sets = c('RNASeq2GeneNorm', 'RPPAArray')) # 101
UpSetR::upset(data = patient_membership_df, sets = c('Methylation', 'RPPAArray')) # 107
# RPPA keeps the common elements low so don't use it!

# keep the 4 largest assays
colnames(patient_membership_df) = assay_names
assays_to_keep = names(sort(colSums(patient_membership_df), decreasing = TRUE))[1:4]
assays_to_keep # CNA, Methyl, RNA, miRNA

cancer_data = cancer_data[,,assays_to_keep]
cancer_data

# keep only the common patients across assays
cancer_data = MultiAssayExperiment::intersectColumns(cancer_data)
nrow(colData(cancer_data)) # 146 patients now

# check that patientIDs are in correct order across assays
for(i in 1:length(cancer_data)) {
  stopifnot(substr(colnames(cancer_data[[i]]), 1, 12) == cancer_data$patientID)
}

# Save data ----
saveRDS(object = cancer_data, file = 'data/cancer_data.rds')

#'###################################################
# Get key "level 4" clinical & pathological data ----
#'###################################################
clinical_vars = TCGAutils::getClinicalNames(diseaseCode = disease_code)
length(clinical_vars) # 19
clinical_var_mat = colData(cancer_data)[, clinical_vars] %>% as.matrix()
saveRDS(object = clinical_var_mat, file = 'data/clinical_var_mat.rds')
