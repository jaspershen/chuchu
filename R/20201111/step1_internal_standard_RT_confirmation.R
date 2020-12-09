##no source
no_source()

sxtTools::setwd_project()

source("R/20201111/step0_parameter_setting.R")

ls()

###set work directory
setwd(file.path("data", file_path, "internal_standard/"))
library(tidyverse)

###extract IS peaks
####positive
extract_is(polarity = "positive",
           path = ".",
           is_table_name = is_table_name)

####negative
extract_is(polarity = "negative",
           path = ".",
           is_table_name = is_table_name)
