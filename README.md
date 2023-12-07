# Partial Rankings of Optimizers - Corresponding R Code

## Introduction
This anonymous repository contains R-code and data sets corresponding to the "Partial Rankings of Optimizers" article. We apply the introduced tests on three examples: BBOB suite, DeepOPS suite, and a benchmarking suite on multi-objective evolutionary algorithms.

The structure of the repository is as follows:
- File _setup_session.R installs all needed R-packages.
- Folder DeepOPS_Analysis contains the files to evaluate the DeepOPS benchmarking suite, see Schneider et.al. (2019).
- Folder BBOB_Analysis contains the files to evaluate the BBOB suite, see Mersmann et.al. (2015) and Hansen et.al. (2010).
- Folder Multi_Objective_Evolutionary_Algorithms_Analysis contains the files to evaluate the benchmarking suite corresponding to Wu et.al. (2023).

The computation has been tested on Linux Ubuntu 20.04.5 with R version 2.2.4.

## Setup
First, install all necessary R packages by downloading and sourcing the file _setup_session.R

Second, download the main.R file into the appropriate folder where you want to analyze the benchmarking suite.

Third, (depending on what you chose in the second step)
- for DeepOps analysis: download the deepops_Schneideretal_data_test_loss.csv file and save it in the save folder as main.R (approximate runtime: 1 second)
- for Multi-Objective Evolutionary Algorithms: download the 3_stages_DMOP_Wuetal_data_no_std_err.csv file and save it in the same folder as main.R (approximate run time: 10 seconds)
- for BBOB analysis: download the bbob_ranking.Rdata file into the same folder as main.R (approximate runtime: 7 hours)
  
Finally, run main.R.

## Explanation of all the produced plots
Running the code results in some automatically generated pdf/tex/jpeg files. They all start with an initial that indicates which calculation they belong to:
- DeepOPS
- MOEA (Multi-Objective Evolutionary Algoritms).
- BBOB_dim_2 (analysis restricted to test functions with dimension 2).
- BBOB_dim_all (analysis on all dimensions of the test functions).


The resulting files are as follows (note that the first page of the file is often blank)
- ..._all_observed.pdf: plots all observed performance orders
- ..._heatmap.jpg: counts how often optimizer i is above optimizer j in all observations
- ..._boxplot_depth.pdf: Boxplot of the calculated ufg depth values
- ..._summary.txt: Prints summary statistics (like mean, standard deviation,...) of the observed ufg depth values.
- ... _matrix_nonedges_intersect_from_high_to_low.txt: Represents which non-edges the posets corresponding to the k, k in natural numbers, highest depth values have in common. TRUE means that they all agree to have a non-edge, false means that there is at least one poset that contains that edge (note that the row optimizer is above the column optimizer).
- ..._plots_from_highest_to_lowest.pfd: Plots all observed posets sorted by their corresponding ufg depth (in descending order)
- ..._plots_intersect_from_highest_to_lowest.pdf: Plots the intersection of the posets corresponding to the k,k natural numbers that have the highest depth values in common. Page 3 shows the intersection of the posets corresponding to the 2 highest depth values (because the first page is empty).
- ..._plots_intersect_from_lowest_to_highest.pfg: Plot the intersection of the posets corresponding to the k, k natural numbers that have the lowest depth values in comon. Page 3 shows the intersection of the posets corresponding to the 2 smallest depth values (because the first page is empty).
- ..._depth_premises.rds: stores the calculation of the ufg depth values.
- ..._total_time_depth_premises.rds: stores the calculation time of the ufg depth values.

## Literatur:
- O. Mersmann, M. Preuss, H. Trautmann, B. Bischl, and C. Weihs. Analyzing the bbob results by means of benchmarking concepts. Evolutionary Computation, 23(1):161–185, 2015.
-  Schneider, L Balles, and P Hennig. Deepobs: A deep learning optimizer benchmark suite. In International Conference on Learning Representations (ICLR 2019). Amerst, MA, 2019
-  Hansen et al 2010
-  Wu etal 2023
