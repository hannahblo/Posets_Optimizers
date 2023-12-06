# Partially Ranking of Optimizers
## Introduction
This anonymous repository contains R-code and data sets corresponding to the "Robust Statistical Comparison of Random Variables with Locally Varying Scale of Measurement" article. We apply the introduced tests on three examples: BBOB suite, DeepOPS suite, and a benchmarking suite on multi-objective evolutionary algorithms.

The structure of the repository is as follows:
- File _setup_session.R installs all needed R-packages.
- Folder DeepOPS_Analysis contains the files to evaluate the DeepOPS benchmarking suite
- Folder BBOB_Analysis contains the files to evaluate the BBOB suite
- Folder Multi_Objective_Evolutionary_Algorithms_Analysis contains the files to evaluate the benchmarking suite corresponding to TODOOOO


## Setup
First, please install all necessary R-packages by sourcing the file _setup_session.R
Second, download the main.R file in the corresponding folder you liked to analyze the benchmarking suite.
Third, 
- for DeepOps analysis: download deepops_Schneideretal_data_test_loss.csv file and save it in the save folder as main.R
- for Multi-Objective Evolutionary Algoritms:download 3_stages_DMOP_Wuetal_data_no_std_err.csv file and save in the same folder as main.R
- for BBOB analysis: download bbob_ranking.Rdata file in the same folder as main.R
Finally, run main.R.
