# paad-survival-bench

The aim of this repo is to benchmark ML survival models (available via [mlr3proba](https://github.com/mlr-org/mlr3proba/)) on the TCGA [PAAD dataset](https://www.cbioportal.org/study/summary?id=paad_tcga_pan_can_atlas_2018) from the PanCancer Atlas project.

- TCGA data download and filter: [tcga_paad.R](https://github.com/bblodfon/paad-survival-bench/blob/main/scripts/tcga_paad.R)
- Data preprocessing: [preprocessing.R](https://github.com/bblodfon/paad-survival-bench/blob/main/scripts/preprocessing.R)
- The [scripts](https://github.com/bblodfon/paad-survival-bench/tree/main/scripts) directory has several benchmarks, with output [results](https://github.com/bblodfon/paad-survival-bench/tree/main/results) and produced [plots](https://github.com/bblodfon/paad-survival-bench/tree/main/img)
