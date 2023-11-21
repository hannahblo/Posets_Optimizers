# blablaba

################################################################################
# R Session
################################################################################

### Information about packages ddandrda. This package is under developement on
# git. Installation can be done by:
# remove.packages("ddandrda")
# install.packages("devtools")
# devtools::install_github("hannahblo/ddandrda")
library(ddandrda)

### All the other R-packages are on CRAN packages (06.10.2023)
library(reshape)
library(dplyr)
library(data.table)
library(ggplot2)
library(forcats)
library(gridExtra)
library(stargazer)
library(prefmod)
library(reshape2)
library(utils)
library(hasseDiagram)

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
      if (all(learner_base < learner_comp)) {
        graph_mat[opt_i, opt_j] <- 1
      }
    }
  }
  return(graph_mat)
}



# this plot function returns summaries of the result
plot_result <- function(names_columns, item_number, depth_value,
                        list_mat_porders_ml,
                        max_plot_number = 100,
                        file_name_add) {

  max_plot_number <- min(max_plot_number, length(list_mat_porders_ml))

  sink(file = paste0("part", file_name_add, "_summery.txt"))
  print(paste0("The minimal value is ", min(depth_value)))
  print(paste0("The maximal value is ", max(depth_value)))
  print(paste0("The mean value is ", mean(depth_value)))
  print(paste0("The standard deviation is ", sd(depth_value)))
  print(paste0("The median is ", median(depth_value)))
  print(paste0("The number of depth value duplicates (reduced by duplicates given by the data) are ", length(depth_value) -
                 length(unique(depth_value))))
  sink(file = NULL)

  ### Distribution of Depth Values
  pdf(paste0("part", file_name_add,"_boxplot_depth.pdf"), onefile = TRUE)
  boxplot(depth_value, main = "Boxplot of the depth values")
  dev.off()



  ## partial orders
  max_depth_index <- sort(depth_value, index.return = TRUE, decreasing = TRUE)$ix
  pdf(paste0("part", file_name_add, "_plots_from_highest_to_lowest.pdf"), onefile = TRUE)
  for (i in max_depth_index[seq(1, max_plot_number)]) {
    mat <- matrix(as.logical(list_mat_porders_ml[[i]]), ncol = item_number)
    colnames(mat) <- rownames(mat) <- names_columns
    hasse(t(mat), parameters = list(arrow = "backward", shape = "roundrect"))
  }
  dev.off()



  ## Intersections (high to low)
  max_depth_index <- sort(depth_value, index.return = TRUE, decreasing = TRUE)$ix
  pdf(paste0("part", file_name_add, "_plots_intersect_from_highest_to_lowest.pdf"), onefile = TRUE)
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
  pdf(paste0("part", file_name_add,"_plots_intersect_from_lowest_to_highest.pdf"), onefile = TRUE)
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
### OLD Version
# # Step 1: Aim to obtain a data frame which following columns
# # - Function id
# # - Optimizer
# # - next columns are the values of the performance measures of interest
#
# # fixed target results
# TR <- read.csv("ERT_Table_Multi.csv")
# # fixed budget results
# BR <- read.csv("FV_Table_Multi.csv")
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
# full_res[is.na(full_res)] <- 0


load("bbob_ranking.Rdata")
#View(data)

# Step 1: filter those data needed
unique(data$algorithm)
optimizer_interest <- as.factor(c("BFGS", "(1+2_m^s) CMA-ES", "BIPOP-CMA-ES", "MOS", "PSO",
                        "RANDOMSEARCH", "FULLNEWUOA", "Nelder-Doerr", "iAMALGAM",
                        "IPOP-CMA-ES", "G3PCX"))
data_filter <- data %>% filter(algorithm %in% optimizer_interest)
full_res <- data.frame(funcId = sort(rep(seq(1,24), 11)),
                       optimizer = rep(optimizer_interest, 24),
                       ERT_2 = rep(FALSE, 11*24),
                       ERT_3 = rep(FALSE, 11*24),
                       ERT_5 = rep(FALSE, 11*24),
                       ERT_10 = rep(FALSE, 11*24),
                       ERT_20 = rep(FALSE, 11*24),
                       ERT_40 = rep(FALSE, 11*24))

for (func_index in unique(data_filter$funcId)) {
  for (dim_index in unique(data_filter$dimension)) {
    for (algo_index in optimizer_interest) {
    data_inner <- data_filter %>%
      filter(dimension == dim_index) %>%
      filter(funcId == func_index) %>%
      filter(algorithm %in% algo_index)
    smallest_precision <- min(unique(data_inner$precision))
    data_inner <- data_inner %>% filter(precision == smallest_precision)
    row_full_res <- intersect(which(full_res$optimizer == algo_index),
                              which(full_res$funcId == func_index))
    full_res[row_full_res, paste0("ERT_", dim_index)] <- data_inner[1, "ert"]
    }
  }
}

colSums(is.na(full_res))
dim(full_res)
unique(full_res[which(is.na(full_res$ERT_2)), ]$optimizer)
unique(full_res[which(is.na(full_res$ERT_5)), ]$optimizer)
# we do not use the algorithm RANDOMSEARCH, BFGS, FULLNEWUOA, G3PCX, PSO and IPOP-CMA-ES
# and (1+2_m^s) CMA-ES
# the comparisons are only done on dimension 2 and 3

optimizer_interest_subset <- as.factor(c("BIPOP-CMA-ES", "MOS", "Nelder-Doerr", "iAMALGAM"))

full_res <- full_res %>% filter(optimizer %in% optimizer_interest_subset)

# Step 2: convert the data frame into a list of partial orders (posets)
funcId <- as.factor(seq(1, 24))
list_graph <- list()


for (id in funcId) {
  single_data_eval <- filter(full_res, funcId == id)
  single_data_eval <- subset(single_data_eval, select = -c(funcId))
  # to apply the convert_to_matrix function it is important that there exists a
  # column with name optimizer which lists the optimizer of interest
  list_graph <- append(list_graph, list(convert_to_matrix(single_data_eval)))
}



################################################################################
#
# PART 1: FIRST IMPRESSION
#
################################################################################
setwd("results/")
number_optimizer <- dim(list_graph[[i]])[1]

### Which edge exists
length(list_graph) # 80
length(unique(list_graph)) # 58
Reduce("|", list_graph)
Reduce("&", list_graph)
Reduce("+", list_graph)

duplicated(list_graph) # no duplications exist


pdf("all_observed.pdf", onefile = TRUE)
for (i in 1:length(list_graph)) {
  mat <- matrix(as.logical(list_graph[[i]]), ncol = number_optimizer)
  colnames(mat) <- rownames(mat) <- colnames(list_graph[[i]])
  hasse(t(mat), parameters = list(arrow = "backward", shape = "roundrect"))
}
dev.off()


# Heatmap
edges <- Reduce("+", list_graph)
colnames(edges) <- rownames(edges) <- c("BS", "CART", "EN", "GBM", "GLM", "LASSO", "RF", "RIDGE")
df_edge_exist <- melt(edges)
df_edge_exist <- df_edge_exist[df_edge_exist$value != 0, ]


jpeg(file = "heatmap_UCI.jpeg")
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

item_number <- 8
names_columns <- colnames(list_graph[[1]])

start_time <- Sys.time()
depth_premises <- ddandrda::compute_ufg_depth_porder(list_graph,
                                                     print_progress_text = TRUE,
                                                     save_ufg_premises = TRUE)
total_time <- Sys.time() - start_time
total_time

# saveRDS(depth_premises, "depth_premises.rds")
# saveRDS(total_time, "total_time_depth_premises.rds")

# TODOOO
plot_result(names_columns =  c("multinom", "ranger", "rpart", "glmnet", "kknn"),
            item_number = number_classifiers,
            depth_value =  depth_value_dupl_22,
            list_mat_porders_ml = list_mat_porders_ml_22,
            file_name_add = "2a_23_obs")
