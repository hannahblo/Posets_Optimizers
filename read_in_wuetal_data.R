library(readr)
library(tidyr)
library(tidyverse)
wuetal_data <- read_delim("DMOP_Wuetal_data.csv", 
                               delim = ";", escape_double = FALSE, trim_ws = TRUE)
#View(wuetal_data)

wuetal_data = wuetal_data[,1:4]
# remove standard deviation
wuetal_data = wuetal_data %>% mutate(MIGD = trimws(str_replace(MIGD, "\\(.*?\\)", "")))
wuetal_data = wuetal_data %>% mutate(MHVD = trimws(str_replace(MHVD, "\\(.*?\\)", "")))
wuetal_data

write.csv(wuetal_data, file = "DMOP_Wuetal_data_no_std_err.csv")
?save
