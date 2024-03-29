# paad-survival-bench

The aim of this repo is to benchmark ML survival models (available via [mlr3proba](https://github.com/mlr-org/mlr3proba/)) on the TCGA [PAAD dataset](https://www.cbioportal.org/study/summary?id=paad_tcga_pan_can_atlas_2018) from the PanCancer Atlas project.

-   **TCGA data** download and filter: [tcga_paad.R](https://github.com/bblodfon/paad-survival-bench/blob/main/scripts/tcga_paad.R)
-   **Data preprocessing**: [preprocessing.R](https://github.com/bblodfon/paad-survival-bench/blob/main/scripts/preprocessing.R)
-   The [scripts](https://github.com/bblodfon/paad-survival-bench/tree/main/scripts) directory has several benchmarks, with some output [results](https://github.com/bblodfon/paad-survival-bench/tree/main/results) stored and the most important produced [plots](https://github.com/bblodfon/paad-survival-bench/tree/main/img).
The most important scripts/investigations are the following:
    - Benchmark **CoxNet**, **Survival Trees** and **Survival Forests** using nested-CV - [script](https://github.com/bblodfon/paad-survival-bench/blob/main/scripts/bench_nestedCV_v4.R)
    - Tuning strategy investigation (Random search vs Bayesian Optimization) using CoxNet and Survival Forests - [script](https://github.com/bblodfon/paad-survival-bench/blob/main/scripts/bayesian_vs_randomsearch.R)
    - **XGBoost** survival learner performance on mRNA dataset - [script](https://github.com/bblodfon/paad-survival-bench/blob/main/scripts/bench_xgboost_v2.R)
    - **CoxPH** baseline performance using clinical features and several resampling strategies - [script](https://github.com/bblodfon/paad-survival-bench/blob/main/scripts/coxph_bench.R)
    - **CoxBoost** (mRNA only and mRNA + clinical) vs CoxPH (clinical) - [script](https://github.com/bblodfon/paad-survival-bench/blob/main/scripts/coxboost_bench.R)
    - **Glmboost** survival learner performance on mRNA dataset - [script](https://github.com/bblodfon/paad-survival-bench/blob/main/scripts/glmboost_bench.R)
    - Wrapper-based **Ensemble Feature Selection (eFS)** per data modality - see [script](https://github.com/bblodfon/paad-survival-bench/blob/main/scripts/fs_mRNA.R) for mRNA data
    - **Task powerset benchmark** after eFS is applied (using simple [CoxPH](https://github.com/bblodfon/paad-survival-bench/blob/main/scripts/coxph_powerset.R) or multiple [learners](https://github.com/bblodfon/paad-survival-bench/blob/main/scripts/bench_powerset.R))
