# paad-survival-bench

The aim of this repo is to benchmark ML survival models (available via [mlr3proba](https://github.com/mlr-org/mlr3proba/)) on the TCGA [PAAD dataset](https://www.cbioportal.org/study/summary?id=paad_tcga_pan_can_atlas_2018) from the PanCancer Atlas project.

-   TCGA data download and filter: [tcga_paad.R](https://github.com/bblodfon/paad-survival-bench/blob/main/scripts/tcga_paad.R)
-   Data preprocessing: [preprocessing.R](https://github.com/bblodfon/paad-survival-bench/blob/main/scripts/preprocessing.R)
-   The [scripts](https://github.com/bblodfon/paad-survival-bench/tree/main/scripts) directory has several benchmarks, with some output [results](https://github.com/bblodfon/paad-survival-bench/tree/main/results) stored and the most important produced [plots](https://github.com/bblodfon/paad-survival-bench/tree/main/img).
The most important scripts/investigations are the following:
  - Benchmark CoxNet, Survival Trees and Survival Forests using nested-CV - [script](https://github.com/bblodfon/paad-survival-bench/blob/main/scripts/bench_nestedCV_v4.R)
  - Tuning strategy investigation (Random search vs Bayesian Optimization) using CoxNet and Survival Forests - [script](https://github.com/bblodfon/paad-survival-bench/blob/main/scripts/bayesian_vs_randomsearch.R)
  - XGBoost survival learner performance on mRNA dataset - [script](https://github.com/bblodfon/paad-survival-bench/blob/main/scripts/bench_xgboost.R)

