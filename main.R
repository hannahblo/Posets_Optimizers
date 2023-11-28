# blablaba


# In what follows is ERT the mean expected running time and FT mean expected target value
################################################################################
# R Session
################################################################################

### Information about packages ddandrda. This package is under developement on
# git. Installation can be done by:
# remove.packages("ddandrda")
# install.packages("devtools")
# devtools::install_github("hannahblo/ddandrda")
library(ddandrda)

### Information about packages oofos. This package is under developement on git.
# Note that to install this package, the R-package gurobi is needed. This is an
# remove.packages("oofos")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("RBGL")
# install.packages("gurobi")
# install.packages("devtools")
# devtools::install_github("schollmeyer/oofos")
library(oofos)

### All the other R-packages are on CRAN packages (06.10.2023)
library(reshape)
library(dplyr)
library(data.table)
library(ggplot2) # bis hier
library(forcats)
library(gridExtra)
library(stargazer)
library(prefmod)
library(reshape2)
library(utils)
library(hasseDiagram) # therefore one needs Rgraphviz and graph packe --> to install this one needs Bioconductor

# setwd("RCode/")


################################################################################
# Functions needed later
################################################################################
# Function which converts to to partial orders
convert_to_matrix <- function(single_data_eval) {
  list_optimizer <- single_data_eval[["optimizer"]]
  number_optimizer <- length(list_optimizer)
  graph_mat <- matrix(rep(0, number_optimizer * number_optimizer),
                      nrow = number_optimizer)
  rownames(graph_mat) <- colnames(graph_mat) <- list_optimizer
  diag(graph_mat) <- 1

  single_data_eval <- subset(single_data_eval, select = -c(optimizer))
  for (opt_i in seq(1, number_optimizer)) {
    learner_base <- single_data_eval[opt_i, ]
    for (opt_j in seq(1, number_optimizer)[-opt_i]) {
      learner_comp <- single_data_eval[opt_j, ]
      difference <- as.numeric(learner_base) - as.numeric(learner_comp)
      difference <- difference[!is.na(difference)]
      if (length(difference) == 0) {
        print(paste0("Optimizer ", list_optimizer[opt_i], " and Optimizer ", list_optimizer[opt_j], " have both in all components Inf."))
      } else if (all(difference == 0)) {
        print(paste0("Optimizer ", list_optimizer[opt_i], " and Optimizer ", list_optimizer[opt_j], " have both in all components same values."))
      } else if (!any(difference < 0)) {
        graph_mat[opt_i, opt_j] <- 1
      }
    }
  }
  return(graph_mat)
}

#####
#
# Preparation of computing the ufg premises
# Function to compute the weights of the fc
get_weighted_representation <- function(x, y = rep(1, dim(x)[1])) {
  ## computes weighted representation of a data matrix x with duplicated rows,
  ##  returns unique(x) together with counts: how often appears the column,
  # mean_y: mean of y in the set of the duplicated columns
  xd <- data.frame(cbind(x, y))
  names(xd)[1] <- "v1"
  v1 <- "v1"
  p <- dim(x)[2]
  result <- as.matrix(plyr::ddply(xd, names(xd[(1:p)]), dplyr::summarise, count = length(v1), mean.y = mean(y), sum.y = sum(y)))
  x_weighted <- result[, (1:p)]
  colnames(x_weighted) <- colnames(x)
  return(list(x_weighted = x_weighted, y_weighted = result[, p + 3], mean_y = result[, p + 2], counts = result[, p + 1]))
}


prepare_ufg_premises <- function(list_mat_porders_ml,
                                 number_items) {

  fc_ml_porder <- ddandrda::compute_conceptual_scaling(input_porder = list_mat_porders_ml)
  # porder_all <- ddandrda::compute_all_partial_orders(number_items, list = FALSE, complemented = TRUE)

  data_context <- get_weighted_representation(fc_ml_porder) # duplication
  n_row_context <- nrow(data_context$x_weighted)
  count_dup <- data_context$counts
  number_obs <- sum(data_context$counts)
  list_porder_premises <- ddandrda::convert_context_to_list(data_context$x_weighted[ ,(1:25)],  complemented = FALSE)

  # whole_context <- rbind(data_context$x_weighted, porder_all) # context of all posets
  # index <- which(!duplicated(whole_context))
  # whole_context <- whole_context[index,]
  return(list(count_dup = count_dup,
              number_obs = number_obs,
              whole_context = data_context$x_weighted,
              list_porder_premises = list_porder_premises,
              n_row_context = n_row_context))
}



# Computing the ufg depth based on already computed premises
compute_ufg_exist_premises <- function(poset_interest, ufg_premises,
                                       prep_ufg_premises) {

  emp_prob <- prep_ufg_premises$count_dup / prep_ufg_premises$number_obs
  depth_ufg <- rep(0, length(poset_interest))
  constant_c <- 0

  for (i in 1:length(ufg_premises)) {
    # print(paste0("Iteration ", i,  " of ", dim(ufg_premises)[1]))
    index_premise <- ufg_premises[[i]]
    prod_emp_ufg <- prod(emp_prob[index_premise])
    concl_ufg <- ddandrda::test_porder_in_concl(prep_ufg_premises$list_porder_premises[index_premise], poset_interest) * 1

    depth_ufg <- depth_ufg + concl_ufg * prod_emp_ufg
    constant_c <- constant_c + prod_emp_ufg
  }

  depth_value <- depth_ufg / constant_c

  return(depth_value)

}
######

# this plot function returns summaries of the result
plot_result <- function(names_columns, item_number, depth_value,
                        list_mat_porders_ml,
                        max_plot_number = 100,
                        file_name_add) {

  max_plot_number <- min(max_plot_number, length(list_mat_porders_ml))

  sink(file = paste0(file_name_add, "_summary.txt"))
  print(paste0("The minimal value is ", min(depth_value)))
  print(paste0("The maximal value is ", max(depth_value)))
  print(paste0("The mean value is ", mean(depth_value)))
  print(paste0("The standard deviation is ", sd(depth_value)))
  print(paste0("The median is ", median(depth_value)))
  print(paste0("The number of depth value duplicates (reduced by duplicates given by the data) are ", length(depth_value) -
                 length(unique(depth_value))))
  sink(file = NULL)

  ### Distribution of Depth Values
  pdf(paste0( file_name_add,"_boxplot_depth.pdf"), onefile = TRUE)
  boxplot(depth_value, main = "Boxplot of the depth values")
  dev.off()



  ## partial orders
  max_depth_index <- sort(depth_value, index.return = TRUE, decreasing = TRUE)$ix
  pdf(paste0(file_name_add, "_plots_from_highest_to_lowest.pdf"), onefile = TRUE)
  for (i in max_depth_index[seq(1, max_plot_number)]) {
    mat <- matrix(as.logical(list_mat_porders_ml[[i]]), ncol = item_number)
    colnames(mat) <- rownames(mat) <- names_columns
    hasse(t(mat), parameters = list(arrow = "backward", shape = "roundrect"))
  }
  dev.off()



  ## Intersections (high to low)
  max_depth_index <- sort(depth_value, index.return = TRUE, decreasing = TRUE)$ix
  pdf(paste0(file_name_add, "_plots_intersect_from_highest_to_lowest.pdf"), onefile = TRUE)
  for (i in 1:max_plot_number) {
    intersect <- matrix(rep(TRUE, item_number * item_number), ncol = item_number)
    for (j in seq(1, i)) {
      intersect <- intersect & matrix(as.logical(list_mat_porders_ml[[max_depth_index[j]]]), ncol = item_number)
    }
    colnames(intersect) <- rownames(intersect) <- names_columns
    hasse(t(intersect), parameters = list(arrow = "backward", shape = "roundrect"))
  }
  dev.off()



  ## Intersections (low to high)
  min_depth_index <- sort(depth_value, index.return = TRUE, decreasing = FALSE)$ix
  pdf(paste0(file_name_add,"_plots_intersect_from_lowest_to_highest.pdf"), onefile = TRUE)
  for (i in 1:max_plot_number) {
    intersect <- matrix(rep(TRUE, item_number * item_number), ncol = item_number)
    for (j in seq(1, i)) {
      intersect <- intersect & matrix(as.logical(list_mat_porders_ml[[min_depth_index[j]]]), ncol = item_number)
    }
    colnames(intersect) <- rownames(intersect) <- names_columns
    hasse(t(intersect), parameters = list(arrow = "backward", shape = "roundrect"))
  }
  dev.off()

}

################################################################################
#
# PART 0: POSET VALUED DATA SET CONSTRUCTION
#
################################################################################
# Step 1: Aim to obtain a data frame which following columns
# - Function id
# - Optimizer
# - next columns are the values of the performance measures of interest


### Version 1
### We use the ERT and FV obtained from https://iohanalyzer.liacs.nl/ (accessed: 18.11.2023)
# version_computation <- "_1"
# name_plots <- "1"
#
# # fixed target results
# TR <- read.csv("ERT_Table_Multi_18_11_23.csv")
# # fixed budget results
# BR <- read.csv("FV_Table_Multi_18_11_23.csv")
#
# full_res <- merge(BR, TR, by = c("funcId", "ID"))
# full_res <- full_res[, c("funcId", "ID", "mean.x", "mean.y")]
# full_res[ ,"funcId"] <- as.factor(full_res[, "funcId"])
# colnames(full_res) <- c("funcId", "optimizer", "mean_tr", "mean_br")
#
# dim(full_res)
# # when looking at the dataframe we observe that for each function each optimizer
# # is evaluated by both performance measures.
# # Thus in what follow we build based on each function ID one single poset
#
# # Problem -> NAs exists. We set in the following NA to zero
# full_res[is.na(full_res)] <- Inf
#
# funcId <- as.factor(seq(1, 24))



### Version 2
### We use the data given by https://dl.p-value.net/2013-ecj_benchmarking/
### corresponding to Mersmann, 0. etal (2015): Analyzing the BBOB results by
### means of benchmarking concepts, Evolutionary Compuation
###
### The selected algorithms correspond to these selected in Section 4.3 of the upper
### paper
###
### The idea is that for each function we use the ERT values for different dimensions
### as multiple performance criteria.
# version_computation <- "_2"
# name_plots <- "2"
#
# load("bbob_ranking.Rdata")
# # View(data)
#
# # Step 1: filter those data needed
# unique(data$algorithm)
# optimizer_interest <- as.factor(c("BFGS", "(1+2_m^s) CMA-ES", "BIPOP-CMA-ES", "MOS", "PSO",
#                         "RANDOMSEARCH", "FULLNEWUOA", "Nelder-Doerr", "iAMALGAM",
#                         "IPOP-CMA-ES", "G3PCX"))
# data_filter <- data %>% filter(algorithm %in% optimizer_interest)
# full_res <- data.frame(funcId = sort(rep(seq(1,24), 11)),
#                        optimizer = rep(optimizer_interest, 24),
#                        ERT_2 = rep(FALSE, 11*24),
#                        ERT_3 = rep(FALSE, 11*24),
#                        ERT_5 = rep(FALSE, 11*24),
#                        ERT_10 = rep(FALSE, 11*24),
#                        ERT_20 = rep(FALSE, 11*24),
#                        ERT_40 = rep(FALSE, 11*24))
#
# for (func_index in unique(data_filter$funcId)) {
#   for (dim_index in unique(data_filter$dimension)) {
#     for (algo_index in optimizer_interest) {
#     data_inner <- data_filter %>%
#       filter(dimension == dim_index) %>%
#       filter(funcId == func_index) %>%
#       filter(algorithm %in% algo_index)
#     precision <- min(unique(data_inner$precision))
#     data_inner <- data_inner %>% filter(precision == precision)
#     row_full_res <- intersect(which(full_res$optimizer == algo_index),
#                               which(full_res$funcId == func_index))
#     full_res[row_full_res, paste0("ERT_", dim_index)] <- data_inner[1, "ert"]
#     }
#   }
# }
# funcId <- as.factor(seq(1, 24))
#
# # NA means that the optimizer did never reach the necessary precision, thus we set
# # them all to Inf (infinity)
# colSums(is.na(full_res))
# full_res[is.na(full_res)] <- Inf
# colSums(is.na(full_res))
# dim(full_res)
# unique(full_res[which(is.na(full_res$ERT_2)), ]$optimizer)
# unique(full_res[which(is.na(full_res$ERT_5)), ]$optimizer)





### Version 3
### We use the data given by https://dl.p-value.net/2013-ecj_benchmarking/
### corresponding to Mersmann, 0. etal (2015): Analyzing the BBOB results by
### means of benchmarking concepts, Evolutionary Compuation
###
### The selected algorithms correspond to these selected in Section 4.3 of the upper
### paper
###
### Looking at Version 2 we obtain that for two funcID there exist two optimizers where all ERT criteria
### equal Inf as they never reached the necessary
### Thus, setting these two to incomparable analogously as we define it, is not
### necessary meaningful. Thus, we drop these two functIds (that are 18 and 24)
### in the following and discuss the same analysis as above
# version_computation <- "_3"
# name_plots <- "3"
#
# load("bbob_ranking.Rdata")
# # View(data)
#
# # Step 1: filter those data needed
# unique(data$algorithm)
# optimizer_interest <- as.factor(c("BFGS", "(1+2_m^s) CMA-ES", "BIPOP-CMA-ES", "MOS", "PSO",
#                                   "RANDOMSEARCH", "FULLNEWUOA", "Nelder-Doerr", "iAMALGAM",
#                                   "IPOP-CMA-ES", "G3PCX"))
# data_filter <- data %>% filter(algorithm %in% optimizer_interest)
# full_res <- data.frame(funcId = sort(rep(seq(1,24)[-c(18,24)], 11)),
#                        optimizer = rep(optimizer_interest, 22),
#                        ERT_2 = rep(FALSE, 11*22),
#                        ERT_3 = rep(FALSE, 11*22),
#                        ERT_5 = rep(FALSE, 11*22),
#                        ERT_10 = rep(FALSE, 11*22),
#                        ERT_20 = rep(FALSE, 11*22),
#                        ERT_40 = rep(FALSE, 11*22))
#
# for (func_index in unique(data_filter$funcId)[-c(18,24)]) {
#   for (dim_index in unique(data_filter$dimension)) {
#     for (algo_index in optimizer_interest) {
#       data_inner <- data_filter %>%
#         filter(dimension == dim_index) %>%
#         filter(funcId == func_index) %>%
#         filter(algorithm %in% algo_index)
#       precision <- max(unique(data_inner$precision))
#       data_inner <- data_inner %>% filter(precision == precision)
#       row_full_res <- intersect(which(full_res$optimizer == algo_index),
#                                 which(full_res$funcId == func_index))
#       full_res[row_full_res, paste0("ERT_", dim_index)] <- data_inner[1, "ert"]
#     }
#   }
# }
# funcId <- as.factor(seq(1, 24)[-c(18,24)])
#
# # NA means that the optimizer did never reach the necessary precision, thus we set
# # them all to Inf (infinity)
# colSums(is.na(full_res))
# full_res[is.na(full_res)] <- Inf
# colSums(is.na(full_res))
# dim(full_res)
# unique(full_res[which(is.na(full_res$ERT_2)), ]$optimizer)
# unique(full_res[which(is.na(full_res$ERT_5)), ]$optimizer)



### Version 4
### We use the data given by https://dl.p-value.net/2013-ecj_benchmarking/
### corresponding to Mersmann, 0. etal (2015): Analyzing the BBOB results by
### means of benchmarking concepts, Evolutionary Compuation
###
### The selected algorithms correspond to these selected in Section 4.3 of the upper
### paper
###
### The idea is that for each function we use the ERT values of dimension 2 and 3
### as multiple performance criteria.
# version_computation <- "_4"
# name_plots <- "4"
#
# load("bbob_ranking.Rdata")
# # View(data)
#
# # Step 1: filter those data needed
# unique(data$algorithm)
# optimizer_interest <- as.factor(c("BFGS", "(1+2_m^s) CMA-ES", "BIPOP-CMA-ES", "MOS", "PSO",
#                         "RANDOMSEARCH", "FULLNEWUOA", "Nelder-Doerr", "iAMALGAM",
#                         "IPOP-CMA-ES", "G3PCX"))
# data_filter <- data %>% filter(algorithm %in% optimizer_interest)
# full_res <- data.frame(funcId = sort(rep(seq(1,24), 11)),
#                        optimizer = rep(optimizer_interest, 24),
#                        ERT_2 = rep(FALSE, 11*24),
#                        ERT_3 = rep(FALSE, 11*24))
#                        # ERT_5 = rep(FALSE, 11*24),
#                        # ERT_10 = rep(FALSE, 11*24),
#                        # ERT_20 = rep(FALSE, 11*24),
#                        # ERT_40 = rep(FALSE, 11*24))
#
# for (func_index in unique(data_filter$funcId)) {
#   for (dim_index in c(2,3)) {
#     for (algo_index in optimizer_interest) {
#     data_inner <- data_filter %>%
#       filter(dimension == dim_index) %>%
#       filter(funcId == func_index) %>%
#       filter(algorithm %in% algo_index)
#     precision <- min(unique(data_inner$precision))
#     data_inner <- data_inner %>% filter(precision == precision)
#     row_full_res <- intersect(which(full_res$optimizer == algo_index),
#                               which(full_res$funcId == func_index))
#     full_res[row_full_res, paste0("ERT_", dim_index)] <- data_inner[1, "ert"]
#     }
#   }
# }
# funcId <- as.factor(seq(1, 24))
#
# # NA means that the optimizer did never reach the necessary precision, thus we set
# # them all to Inf (infinity)
# colSums(is.na(full_res))
# full_res[is.na(full_res)] <- Inf
# colSums(is.na(full_res))
# dim(full_res)
# unique(full_res[which(is.na(full_res$ERT_2)), ]$optimizer)
# unique(full_res[which(is.na(full_res$ERT_5)), ]$optimizer)






### Version 5
### We use the data given by https://dl.p-value.net/2013-ecj_benchmarking/
### corresponding to Mersmann, 0. etal (2015): Analyzing the BBOB results by
### means of benchmarking concepts, Evolutionary Compuation
###
### The selected algorithms correspond to these selected in Section 4.3 of the upper
### paper
###
### Looking at Version 4 we obtain that for two funcID there exist two optimizers where all ERT criteria
### equal Inf as they never reached the necessary
### Thus, setting these two to incomparable analogously as we define it, is not
### necessary meaningful. Thus, we drop these functIds (that are c(2, 5, 7, 13, 16, 17, 18, 23, 24))
### in the following and discuss the same analysis as above
# version_computation <- "_5"
# name_plots <- "5"
#
# load("bbob_ranking.Rdata")
# # View(data)
#
# # Step 1: filter those data needed
# unique(data$algorithm)
# optimizer_interest <- as.factor(c("BFGS", "(1+2_m^s) CMA-ES", "BIPOP-CMA-ES", "MOS", "PSO",
#                                   "RANDOMSEARCH", "FULLNEWUOA", "Nelder-Doerr", "iAMALGAM",
#                                   "IPOP-CMA-ES", "G3PCX"))
# data_filter <- data %>% filter(algorithm %in% optimizer_interest)
# full_res <- data.frame(funcId = sort(rep(seq(1,24)[-c(2, 5, 7, 13, 16, 17, 18, 23, 24)], 11)),
#                        optimizer = rep(optimizer_interest, 15),
#                        ERT_2 = rep(FALSE, 11*15),
#                          ERT_3 = rep(FALSE, 11*15))
#                        # ERT_5 = rep(FALSE, 11*15),
#                        # ERT_10 = rep(FALSE, 11*15),
#                        # ERT_20 = rep(FALSE, 11*15),
#                        # ERT_40 = rep(FALSE, 11*15))
#
# for (func_index in unique(data_filter$funcId)[-c(18,24)]) {
#   for (dim_index in unique(data_filter$dimension)) {
#     for (algo_index in optimizer_interest) {
#       data_inner <- data_filter %>%
#         filter(dimension == dim_index) %>%
#         filter(funcId == func_index) %>%
#         filter(algorithm %in% algo_index)
#       precision <- max(unique(data_inner$precision))
#       data_inner <- data_inner %>% filter(precision == precision)
#       row_full_res <- intersect(which(full_res$optimizer == algo_index),
#                                 which(full_res$funcId == func_index))
#       full_res[row_full_res, paste0("ERT_", dim_index)] <- data_inner[1, "ert"]
#     }
#   }
# }
# funcId <- as.factor(seq(1, 24)[-c(2, 5, 7, 13, 16, 17, 18, 23, 24)])
#
# # NA means that the optimizer did never reach the necessary precision, thus we set
# # them all to Inf (infinity)
# colSums(is.na(full_res))
# full_res[is.na(full_res)] <- Inf
# colSums(is.na(full_res))
# dim(full_res)
# unique(full_res[which(is.na(full_res$ERT_2)), ]$optimizer)
# unique(full_res[which(is.na(full_res$ERT_5)), ]$optimizer)





# TODO
# ACHTUNG!!!! MACHT DAS ÜBERHAUPT SINN
### Version 6
### We use the data given by https://dl.p-value.net/2013-ecj_benchmarking/
### corresponding to Mersmann, 0. etal (2015): Analyzing the BBOB results by
### means of benchmarking concepts, Evolutionary Compuation
###
### The selected algorithms correspond to these selected in Section 4.3 of the upper
### paper.
###
### We computed the precision values by hand and added this as a performance evaluation on dimnesion 2.
###
### The used performance measures are all possible ERT values of different dimensions and the precision value
# version_computation <- "_6"
# name_plots <- "6"
# load("bbob_ranking.Rdata")
#
# unique(data$algorithm)
# optimizer_interest <- as.factor(c("BFGS", "(1+2_m^s) CMA-ES", "BIPOP-CMA-ES", "MOS", "PSO",
#                                   "RANDOMSEARCH", "FULLNEWUOA", "Nelder-Doerr", "iAMALGAM",
#                                   "IPOP-CMA-ES", "G3PCX"))
# data_filter <- data %>% filter(algorithm %in% optimizer_interest)
# full_res <- data.frame(funcId = sort(rep(seq(1,24), 11)),
#                        optimizer = rep(optimizer_interest, 24),
#                        ERT_2 = rep(FALSE, 11*24),
#                        ERT_3 = rep(FALSE, 11*24),
#                        ERT_5 = rep(FALSE, 11*24),
#                        ERT_10 = rep(FALSE, 11*24),
#                        ERT_20 = rep(FALSE, 11*24),
#                        ERT_40 = rep(FALSE, 11*24))
#
# for (func_index in unique(data_filter$funcId)) {
#   for (dim_index in unique(data_filter$dimension)) {
#     for (algo_index in optimizer_interest) {
#       data_inner <- data_filter %>%
#         filter(dimension == dim_index) %>%
#         filter(funcId == func_index) %>%
#         filter(algorithm %in% algo_index)
#       smallest_precision <- min(unique(data_inner$precision))
#       data_inner <- data_inner %>% filter(precision == smallest_precision)
#       row_full_res <- intersect(which(full_res$optimizer == algo_index),
#                                 which(full_res$funcId == func_index))
#       full_res[row_full_res, paste0("ERT_", dim_index)] <- data_inner[1, "ert"]
#     }
#   }
# }
#
# colSums(is.na(full_res))
# full_res[is.na(full_res)] <- Inf
# colSums(is.na(full_res))
#
# # begin edit for function value
# # compute best function value for given time (1000)
#
# ## lowest precision that an algo has achieved at time 1000
# # ERTs higher than 1000 time units
# full_res_sub <- full_res[,-c(1,2)]
# full_res_higher <- apply(full_res_sub, 2,  function(x) x < 1000)
# full_res_fv_bin <- cbind(full_res[,c(1,2)], full_res_higher)
# # replace NAs (~infinite runtime) by False
# full_res_fv_bin[is.na(full_res_fv_bin) == TRUE] <- FALSE
# # now count precision (function values) for which ERT < 1000
# function_value = apply(full_res_fv_bin[,-c(1,2)], 1, sum)
# # and append to full_res data frame
# full_res$function_value = function_value
#
# full_res <- full_res[, c(1,2,3,9)]
#
# ## Interpretation:
# # The variable "function value" describes the number of precision values (thresholds)
# # achieved by an optimizer within a fixed time budget of 1000 time units (seconds?)
#
# funcId <- as.factor(seq(1, 24))



### Version 7
###
### We use the data given by https://dl.p-value.net/2013-ecj_benchmarking/
### corresponding to Mersmann, 0. etal (2015): Analyzing the BBOB results by
### means of benchmarking concepts, Evolutionary Computaion
###
### The selected algorithms correspond to these selected in Section 4.3 of the upper
### paper.
###
### We only use the functions with 2 dimensions.
### We computed the precision values by hand and added this as a performance evaluation
###
### The used performance measures are ERT based on all precision values and the precision value
# version_computation <- "_7"
# name_plots <- "7"
# load("bbob_ranking.Rdata")
# unique(data$algorithm)
# unique(data$dimension)
# optimizer_interest <- as.factor(c("BFGS", "(1+2_m^s) CMA-ES", "BIPOP-CMA-ES", "MOS", "PSO",
#                                   "RANDOMSEARCH", "FULLNEWUOA", "Nelder-Doerr", "iAMALGAM",
#                                   "IPOP-CMA-ES", "G3PCX"))
# data_filter <- data %>% filter(algorithm %in% optimizer_interest)
# full_res <- data.frame(funcId = sort(rep(seq(1,24), 11)),
#                        optimizer = rep(optimizer_interest, 24)
# )
# for (func_index in unique(data_filter$funcId)) {
#   for (prec_index in unique(data_filter$precision)) {
#     for (algo_index in optimizer_interest) {
#       data_inner <- data_filter %>%
#         filter(precision == prec_index) %>%
#         filter(funcId == func_index) %>%
#         filter(algorithm %in% algo_index)
#       smallest_dimension <- 2 #min(unique(data_inner$dimension))
#       data_inner <- data_inner %>% filter(dimension == smallest_dimension)
#       row_full_res <- intersect(which(full_res$optimizer == algo_index),
#                                 which(full_res$funcId == func_index))
#       full_res[row_full_res, paste0("ERT_", prec_index)] <- data_inner[1, "ert"]
#     }
#   }
# }
#
# View(full_res)
# colSums(is.na(full_res))
# dim(full_res)
#
# full_res <- full_res %>% filter(optimizer %in% optimizer_interest)
#
# # begin edit for function value
# # compute best function value for given time (5000)
# # lowest precision that an algo has achieved at time 5000
# # ERTs higher than 1000 time units
# full_res_sub <- full_res[,-c(1,2)]
# full_res_higher <- apply(full_res_sub, 2,  function(x) x < 5000)
# full_res_fv_bin <- cbind(full_res[,c(1,2)], full_res_higher)
# # replace NAs (~infinite runtime) by False
# full_res_fv_bin[is.na(full_res_fv_bin) == TRUE] <- FALSE
# # now count precision (function values) for which ERT < 5000
# function_value = apply(full_res_fv_bin[,-c(1,2)], 1, sum)
# # and append to full_res data frame
# full_res$function_value = function_value
#
# colSums(is.na(full_res))
# full_res[is.na(full_res)] <- Inf
# colSums(is.na(full_res))
#
# funcId <- as.factor(seq(1, 24))



### Version 8
###
### We use the data given by https://dl.p-value.net/2013-ecj_benchmarking/
### corresponding to Mersmann, 0. etal (2015): Analyzing the BBOB results by
### means of benchmarking concepts, Evolutionary Computaion
###
### The selected algorithms correspond to these selected in Section 4.3 of the upper
### paper.
###
### We only use the functions with 2 dimensions.
### We computed the precision values by hand and added this as a performance evaluation
###
### The used performance measures are ERT based on precision 0.001 and the precision value
# version_computation <- "_8"
# name_plots <- "8"
# load("bbob_ranking.Rdata")
# unique(data$algorithm)
# unique(data$dimension)
# optimizer_interest <- as.factor(c("BFGS", "(1+2_m^s) CMA-ES", "BIPOP-CMA-ES", "MOS", "PSO",
#                                   "RANDOMSEARCH", "FULLNEWUOA", "Nelder-Doerr", "iAMALGAM",
#                                   "IPOP-CMA-ES", "G3PCX"))
# data_filter <- data %>% filter(algorithm %in% optimizer_interest)
# full_res <- data.frame(funcId = sort(rep(seq(1,24), 11)),
#                        optimizer = rep(optimizer_interest, 24)
# )
# for (func_index in unique(data_filter$funcId)) {
#   for (prec_index in unique(data_filter$precision)) {
#     for (algo_index in optimizer_interest) {
#       data_inner <- data_filter %>%
#         filter(precision == prec_index) %>%
#         filter(funcId == func_index) %>%
#         filter(algorithm %in% algo_index)
#       smallest_dimension <- 2 #min(unique(data_inner$dimension))
#       data_inner <- data_inner %>% filter(dimension == smallest_dimension)
#       row_full_res <- intersect(which(full_res$optimizer == algo_index),
#                                 which(full_res$funcId == func_index))
#       full_res[row_full_res, paste0("ERT_", prec_index)] <- data_inner[1, "ert"]
#     }
#   }
# }
#
# View(full_res)
# colSums(is.na(full_res))
# dim(full_res)
#
# full_res <- full_res %>% filter(optimizer %in% optimizer_interest)
#
# # begin edit for function value
# # compute best function value for given time (5000)
# # lowest precision that an algo has achieved at time 5000
# # ERTs higher than 1000 time units
# full_res_sub <- full_res[,-c(1,2)]
# full_res_higher <- apply(full_res_sub, 2,  function(x) x < 5000)
# full_res_fv_bin <- cbind(full_res[,c(1,2)], full_res_higher)
# # replace NAs (~infinite runtime) by False
# full_res_fv_bin[is.na(full_res_fv_bin) == TRUE] <- FALSE
# # now count precision (function values) for which ERT < 5000
# function_value = apply(full_res_fv_bin[,-c(1,2)], 1, sum)
# # and append to full_res data frame
# full_res$function_value = function_value
#
# colSums(is.na(full_res))
# full_res[is.na(full_res)] <- Inf
# colSums(is.na(full_res))
#
# full_res <- full_res[, -seq(4,8)]
#
# funcId <- as.factor(seq(1, 24))


### Version 9
### using the evaluations by Table 2 of
### Frank Schneider, Lukas Balles & Philipp Hennig: (2019): DEEPOBS: A DEEP LEARNING OPTIMIZER BENCHMARK SUITE
# version_computation <- "_9"
# name_plots <- "9"
# version_computation <- "_9"
# data <- read.csv(file = "deepops_Schneideretal_data.csv",sep = ";", header = TRUE)
# data <- data[seq(1, 24), ]
#
# full_res <- data
# full_res[seq(1, 24), 3] <- as.numeric(gsub(",", ".",full_res[seq(1, 24), 3]))
# full_res[seq(1, 24), 4] <- as.numeric(gsub(",", ".",full_res[seq(1, 24), 4]))
# # as.numeric(as.character(full_res[seq(2, 24), 4]))
#
# # the measure speed: the smaller the better
# # the measure loss: the smaller the better
# # the measure accuracy: the larger the better --> 1 - accuracy again, the smaller the better
# # for P3, P4, P6, P7, P8 the accuracy
# full_res[c(seq(7, 12), seq(16,24)), 3] <- 1 - as.numeric(full_res[c(seq(7, 12), seq(16,24)), 3])
# colnames(full_res)[c(1, 2)] <- c("funcId", "optimizer")
#
# funcId <- unique(full_res[, 1])


### Version 10
### using the evaluations by Table 2 and Figure 2(for the order of the losses which are not given by the tables) of
### Frank Schneider, Lukas Balles & Philipp Hennig: (2019): DEEPOBS: A DEEP LEARNING OPTIMIZER BENCHMARK SUITE

version_computation <- "_10"
name_plots <- "10"
data <- read.csv(file = "deepops_Schneideretal_data_test_loss.csv",sep = ";", header = TRUE)
data <- data[seq(1, 24), ]

full_res <- data
full_res[seq(1, 24), 3] <- as.numeric(gsub(",", ".",full_res[seq(1, 24), 3]))
full_res[seq(1, 24), 4] <- as.numeric(gsub(",", ".",full_res[seq(1, 24), 4]))
# as.numeric(as.character(full_res[seq(2, 24), 4]))

# the measure speed: the smaller the better
# the measure loss: the smaller the better
colnames(full_res)[c(1, 2)] <- c("funcId", "optimizer")

funcId <- unique(full_res[, 1])








#### the following procedure is the same for for all upper Versions
setwd(paste0("result", version_computation, "/"))

# Step 2: convert the data frame into a list of partial orders (posets)
list_graph <- list()


for (id in funcId) {
  print(paste0("Now at functionId ", id))
  single_data_eval <- filter(full_res, funcId == id)
  single_data_eval <- subset(single_data_eval, select = -c(funcId))
  # to apply the convert_to_matrix function it is important that there exists a
  # column with name optimizer which lists the optimizer of interest
  list_graph <- append(list_graph, list(convert_to_matrix(single_data_eval)))
}


for (graph in list_graph) {
  if (!ddandrda::test_if_porder(graph)) {
    print("Attention: one is only a preoder and not a poset!")
  }
}






### Comments:
# to Version 1: here we have no problem that two times only Inf vs Inf are compared
#     thus, we do do not have the question how to include such an information
#     And all comparisons are indeed a poset
# to Version 2: for function id 18 and 24 we have the problematic how to deal with
#     Inf compared to Inf (for optimizers BFGS, RANDOMSEARCH; FULLNEWUOA)
# to Version 3: like version 2 just without funcid 18 and 24
# to Version 4: for function id 2, 5, 7, 13, 16, 17, 18, 23, 24 we have the problematic how to deal with
#     Inf compared to Inf (for optimizers BFGS, RANDOMSEARCH; FULLNEWUOA; IPOP-CMA-ES;
#     G3PXC)
# to Version 5: like version 4 just without funcid 2, 5, 7, 13, 16, 17, 18, 23, 24
# to Version 6: --> neither computation version 1 nor computation version 2 worked
# to Version 7: --> neither computation version 1 nor computation version 2 worked
# to Version 8: for function id 1, 18 and 24 we have the problematic how to deal with
#     the same performance measures values for both optimizers (for optimizers BFGS, RANDOMSEARCH; IPOP-CMA-ES; G3PCX; FULLNEWUOA)
# to Version 9: no problems :)
################################################################################
#
# PART 1: FIRST IMPRESSION
#
################################################################################

number_optimizer <- dim(list_graph[[1]])[1]
optimizer_interest <- colnames(list_graph[[1]])

### Which edge exists
length(list_graph) # 80
length(unique(list_graph)) # 58
Reduce("|", list_graph)
Reduce("&", list_graph)
Reduce("+", list_graph)

duplicated(list_graph) # no duplications exist


pdf(paste0(name_plots, "_all_observed.pdf"), onefile = TRUE)
for (i in 1:length(list_graph)) {
  mat <- matrix(as.logical(list_graph[[i]]), ncol = number_optimizer)
  colnames(mat) <- rownames(mat) <- colnames(list_graph[[i]])
  hasse(t(mat), parameters = list(arrow = "backward", shape = "roundrect"))
}
dev.off()


# Heatmap
edges <- Reduce("+", list_graph)
colnames(edges) <- rownames(edges) <- optimizer_interest
df_edge_exist <- melt(edges)
df_edge_exist <- df_edge_exist[df_edge_exist$value != 0, ]


jpeg(file = paste0(name_plots, "_heatmap.jpeg"))
ggplot(df_edge_exist, aes(x = Var1, y = Var2)) +
  geom_raster(aes(fill = value)) +
  scale_fill_gradient(low = "lightcyan1", high = "darkcyan") +
  labs(x = "is below", y = "is above") +
  geom_text(aes(label = value)) +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 0.3),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 13, vjust = -1),
        axis.title.y = element_text(size = 13),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_blank())
dev.off()
# note that the plotting is in some sense counterintutitve since it switches the
# x and y axis! Compare to edges





################################################################################
#
# PART 2: Computation and Evaluation UFG Depth
#
################################################################################

item_number <- dim(list_graph[[1]])[1]

start_time <- Sys.time()
depth_premises <- ddandrda::compute_ufg_depth_porder(list_graph,
                                                     print_progress_text = TRUE, # if progress should not be shown, set to FALSE
                                                     save_ufg_premises = FALSE) # eventuelly to many -> thus not saved
total_time <- Sys.time() - start_time
total_time

saveRDS(depth_premises, "depth_premises.rds")
saveRDS(total_time, "total_time_depth_premises.rds")


plot_result(names_columns =  optimizer_interest,
            item_number = number_optimizer,
            depth_value =  depth_premises$depth_ufg,
            list_mat_porders_ml = list_graph,
            file_name_add = name_plots)





###### DAS UNTEN LÖSCHEN!!!!! ##################################################

################################################################################
# version 2
################################################################################
#
# ### Compute the VC dimension
# # Formal context given by the partial orders in list_mat
fc_ml_porder <- ddandrda::compute_conceptual_scaling(input_porder = list_graph)
# ml_porder_model <- oofos::compute_extent_vc_dimension(fc_ml_porder)
# vc_fc_ml_porder <- gurobi::gurobi(ml_porder_model)
# vc <- vc_fc_ml_porder$objval # 8



### Compute the ufg-depth
# Preparation of the computation, needed as input of
porder_all <- ddandrda::compute_all_partial_orders(item_number, list = FALSE, complemented = TRUE)
list_porder_all <- ddandrda::compute_all_partial_orders(item_number, list = TRUE, complemented = FALSE)

data_context <- get_weighted_representation(fc_ml_porder) # duplication
n_row_context <- nrow(data_context$x_weighted)
count_dup <- data_context$counts
number_obs <- sum(data_context$counts)

list_ml_porder_unique <- ddandrda::convert_context_to_list(data_context$x_weighted[ ,(1:item_number * item_number)],  complemented = FALSE)

whole_context <- rbind(data_context$x_weighted, porder_all) # context of all posets
index <- which(!duplicated(whole_context))
whole_context <- whole_context[index,]



# Computation of S, see article (1)
start_time <- Sys.time()
ufg_premises <- oofos::enumerate_ufg_premises(whole_context, n_row_context) # das ist seltsam
total_time <- Sys.time() - start_time

# saveRDS(total_time, "total_time.rds")
# saveRDS(ufg_premises, "ufg_premises.rds")
# length(ufg_premises)


# ufg depth computation
emp_prob <- count_dup / number_obs
depth_ufg <- rep(0, length(list_ml_porder_unique))
constant_c <- 0

for (i in 1:length(ufg_premises)) {
  # print(paste0("Iteration ", i,  " of ", dim(ufg_premises)[1]))
  index_premise <- ufg_premises[[i]]
  if (length(index_premise) < 2) {
    print(paste0("cardinaltiy ufg_premise is ", length(index_premise)))
  }

  prod_emp_ufg <- prod(emp_prob[index_premise])
  concl_ufg <- test_porder_in_concl(list_ml_porder_unique[index_premise], list_ml_porder_unique) * 1

  depth_ufg <- depth_ufg + concl_ufg * prod_emp_ufg
  constant_c <- constant_c + prod_emp_ufg
}

depth_value <- depth_ufg / constant_c


# Adding duplicate values
depth_value_all <- c()
list_data_all <- vector("list", sum(count_dup))
saving <- 1
for (i in 1:length(depth_value)) {
  for (j in 1:count_dup[i]) {
    list_data_all[[saving]] <- list_ml_porder_unique[[i]]
    saving <- saving + 1
  }
  depth_value_all <- append(depth_value_all, rep(depth_value[i], count_dup[i]))

}


# saveRDS(constant_c, "constant_c.rds")
# saveRDS(depth_ufg, "ufg_depth.rds")
# saveRDS(vc, "vc.rds")
# saveRDS(depth_value_all, "depth_values.rds")


###########
prep_ufg_premises <- prepare_ufg_premises(list_graph, number_items = number_classifiers)
start_time <- Sys.time()
ufg_premises <- oofos::enumerate_ufg_premises(prep_ufg_premises$whole_context, prep_ufg_premises$n_row_context)
total_time <- Sys.time() - start_time
depth_value <- compute_ufg_exist_premises(poset_interest = prep_ufg_premises$list_porder_premises,
                                          ufg_premises,
                                          prep_ufg_premises)
