##no source
no_source()

sxtTools::setwd_project()
source("R/tools.R")

setwd("data/lipid20201020/internal_standard/")
library(tidyverse)

####internal standard
####positive
extract_is(polarity = "positive",
           path = ".",
           is_table_name = "IS_information_20201020.xlsx")

####negative
extract_is(polarity = "negative",
           path = ".",
           is_table_name = "IS_information_20201020.xlsx")



# ##ChE
# ####positive
# extract_is(polarity = "positive",
#            path = ".",
#            is_table_name = "IS_informaiton_ChE.xlsx")
# 
# ####negative
# extract_is(polarity = "negative",
#            path = ".",
#            is_table_name = "IS_informaiton_ChE.xlsx")
