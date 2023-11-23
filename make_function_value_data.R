# This script produces a data frame like full_res in main.R with function value (fv)
# as an additional metric (:= best function value obtained after at least 1000 time units)

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


load("bbob_ranking.Rdata")

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


optimizer_interest <- as.factor(c("BFGS", "(1+2_m^s) CMA-ES", "BIPOP-CMA-ES", "MOS", "PSO",
                                  "RANDOMSEARCH", "FULLNEWUOA", "Nelder-Doerr", "iAMALGAM",
                                  "IPOP-CMA-ES", "G3PCX"))
full_res <- full_res %>% filter(optimizer %in% optimizer_interest)





######## begin edit for function value

## compute best function value for given time (1000)

## lowest precision that an algo has achieved at time 1000
# ERTs higher than 1000 time units
full_res_sub <- full_res[,-c(1,2)]
full_res_higher <- apply(full_res_sub, 2,  function(x) x < 1000)
full_res_fv_bin <- cbind(full_res[,c(1,2)], full_res_higher) 
# replace NAs (~infinite runtime) by False
full_res_fv_bin[is.na(full_res_fv_bin) == TRUE] <- FALSE
# now count precision (function values) for which ERT < 1000
function_value = apply(full_res_fv_bin[,-c(1,2)], 1, sum)
# and append to full_res data frame
full_res$function_value = function_value

## Interpretation: 
# The variable "function value" describes the number of precision values (thresholds) 
# achieved by an optimizer within a fixed time budget of 1000 time units (seconds?)


