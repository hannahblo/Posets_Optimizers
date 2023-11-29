library(readr)
library(tidyr)
library(tidyverse)
wuetal_data <- read_delim("3_stages_DMOP_Wuetal_data.csv", 
                          delim = ";", escape_double = FALSE, trim_ws = TRUE)
#View(wuetal_data)

wuetal_data = wuetal_data[,1:6]
# remove standard deviation
wuetal_data = wuetal_data %>% mutate(MIGD_total = trimws(str_replace(MIGD_total, "\\(.*?\\)", "")))
wuetal_data = wuetal_data %>% mutate(MIGD_1st_stage = trimws(str_replace(MIGD_1st_stage, "\\(.*?\\)", "")))
wuetal_data = wuetal_data %>% mutate(MIGD_2nd_stage = trimws(str_replace(MIGD_2nd_stage, "\\(.*?\\)", "")))
wuetal_data = wuetal_data %>% mutate(MIGD_3rd_stage = trimws(str_replace(MIGD_3rd_stage, "\\(.*?\\)", "")))
#wuetal_data

write.csv(wuetal_data, file = "3_stages_DMOP_Wuetal_data_no_std_err.csv")
