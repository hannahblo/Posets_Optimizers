# Analysis of the DeepOPS Data given by
# Frank Schneider, Lukas Balles & Philipp Hennig (2019):  DEEPOBS: A DEEP
# LEARNING OPTIMIZER BENCHMARK SUITE
# Table 2 and Figure 2(for the order of the losses which are not given by the tables)

# As performance measures we use
# - the measure speed: the smaller the better
# - the measure loss: the smaller the better


# This code is (partially) copied from
# Hannah Blocher, Georg Schollmeyer, Christoph Jansen, and Malte Nalenz. Depth functions for
# partial orders with a descriptive analysis of machine learning algorithms. In Enrique Miranda,
# Ignacio Montes, Erik Quaeghebeur, and Barbara Vantaggi (eds.), Proceedings of the Thirteenth
# International Symposium on Imprecise Probability: Theories and Applications, volume 215 of
# Proceedings of Machine Learning Research, pp. 59–71. PMLR, 11–14 Jul 2023.
# and the code provided by
# https://github.com/hannahblo/23_Performance_Analysis_ML_Algorithms (accessed: 07.11.2023)


################################################################################
# R Session
################################################################################
library(ddandrda)
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
        print(paste0("Optimizer ", list_optimizer[opt_i], " and Optimizer ",
                     list_optimizer[opt_j], " have both in all components Inf."))
      } else if (all(difference == 0)) {
        print(paste0("Optimizer ", list_optimizer[opt_i], " and Optimizer ",
                     list_optimizer[opt_j], " have both in all components same values."))
      } else if (!any(difference < 0)) {
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

  sink(file = paste0(file_name_add, "_summary.txt"))
  print(paste0("The minimal value is ", min(depth_value)))
  print(paste0("The maximal value is ", max(depth_value)))
  print(paste0("The mean value is ", mean(depth_value)))
  print(paste0("The standard deviation is ", sd(depth_value)))
  print(paste0("The median is ", median(depth_value)))
  print(paste0("The number of depth value duplicates (reduced by duplicates given by the data) are ",
               length(depth_value) - length(unique(depth_value))))
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

  ## Nonedges (high to low)
  max_depth_index <- sort(depth_value, index.return = TRUE, decreasing = TRUE)$ix
  sink(file = paste0(file_name_add, "_matrix_nonedges_intersect_from_highest_to_lowest.txt")) # hasse diagrams cannot be plotted -> not necessarily antisymmetric
  for (i in 1:max_plot_number) {
    intersect <- matrix(rep(TRUE, item_number * item_number), ncol = item_number)
    for (j in seq(1, i)) {
      intersect <- intersect & matrix(!as.logical(list_mat_porders_ml[[max_depth_index[j]]]), ncol = item_number)
    }
    colnames(intersect) <- rownames(intersect) <- names_columns
    print(paste0("\n . \n .\n The ", i, " deepest depth values have these nonedges in common."))
    print(intersect)
    # hasse(t(intersect), parameters = list(arrow = "backward", shape = "roundrect"))
  }
  sink(file = NULL)


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
# we resort the data sets such that each row represents one algorithm,
# one test function and lists in the row all performance values
# the expected running time for precision 0.001 is given, we will compute later
# he number of precision levels achieved overall
#
################################################################################

### Step 0: state the prefix name of the saved documents
name_plots <- "DeepOPS_"

### Step 1: Read the data and convert them to the correct data type
data <- read.csv(file = "deepops_Schneideretal_data_test_loss.csv",sep = ";", header = TRUE)
data <- data[seq(1, 24), ]

full_res <- data
full_res[seq(1, 24), 3] <- as.numeric(gsub(",", ".",full_res[seq(1, 24), 3]))
full_res[seq(1, 24), 4] <- as.numeric(gsub(",", ".",full_res[seq(1, 24), 4]))
colnames(full_res)[c(1, 2)] <- c("funcId", "optimizer")


### Step 2
# convert the data frame into a list of partial orders (posets)
# we say that an optimizer i is better than optimizer j if and only if at least one
# performance evaluation states that optimizer i is better and every other performance
# measure states that optimizer i is not worse than optimizer j.
funcId <- unique(full_res[, 1])

list_graph <- list()

for (id in funcId) {
  print(paste0("Now at functionId ", id))
  single_data_eval <- filter(full_res, funcId == id)
  single_data_eval <- subset(single_data_eval, select = -c(funcId))
  # to apply the convert_to_matrix function it is important that there exists a
  # column with name optimizer which lists the optimizer of interest
  list_graph <- append(list_graph, list(convert_to_matrix(single_data_eval)))
}

# check if everythink is correct and we observe posets
for (graph in list_graph) {
  if (!ddandrda::test_if_porder(graph)) {
    print("Attention: one is only a preoder and not a poset!")
  }
}



### Step 3: First impression of the observed 21 posets
number_optimizer <- dim(list_graph[[1]])[1]
optimizer_interest <- colnames(list_graph[[1]])

# Which edge exists
length(list_graph)
length(unique(list_graph))
Reduce("|", list_graph)
Reduce("&", list_graph)
Reduce("+", list_graph)

duplicated(list_graph) # no duplications exist


pdf(paste0(name_plots, "_all_observed.pdf"), onefile = TRUE)
for (i in 1:length(list_graph)) {
  # print(i)
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



#### Step 4: Computation of the ufg depth
# Hannah Blocher, Georg Schollmeyer, Christoph Jansen, and Malte Nalenz. Depth functions for
# partial orders with a descriptive analysis of machine learning algorithms. In Enrique Miranda,
# Ignacio Montes, Erik Quaeghebeur, and Barbara Vantaggi (eds.), Proceedings of the Thirteenth
# International Symposium on Imprecise Probability: Theories and Applications, volume 215 of
# Proceedings of Machine Learning Research, pp. 59–71. PMLR, 11–14 Jul 2023.
# and the code provided by
# https://github.com/hannahblo/23_Performance_Analysis_ML_Algorithms (accessed: 07.11.2023)
item_number <- dim(list_graph[[1]])[1]

start_time <- Sys.time()
depth_premises <- ddandrda::compute_ufg_depth_porder(list_graph,
                                                     print_progress_text = TRUE, # if progress should not be shown, set to FALSE
                                                     save_ufg_premises = FALSE) # eventually to many -> thus not saved
total_time <- Sys.time() - start_time
total_time

saveRDS(depth_premises, paste0(name_plots, "depth_premises.rds"))
saveRDS(total_time, paste0(name_plots, "total_time_depth_premises.rds"))


plot_result(names_columns =  optimizer_interest,
            item_number = number_optimizer,
            depth_value =  depth_premises$depth_ufg,
            list_mat_porders_ml = list_graph,
            file_name_add = name_plots)
