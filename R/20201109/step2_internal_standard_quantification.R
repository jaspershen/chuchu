##no source
no_source()

sxtTools::setwd_project()
source("R/tools.R")

setwd("data/lipid20201020/IS_quantification/POS/")
library(tidyverse)

###positive
trans_is_table(is_table_name = "IS_information_20201020.xlsx", polarity = "positive")

extract_targeted_peaks(
  path = ".",
  targeted_targeted_peak_table_name = "is_table.xlsx",
  sample_numer_show = 3,
  fit.gaussian = TRUE,
  integrate_xcms = FALSE
)


###negative
sxtTools::setwd_project()
source("R/tools.R")

setwd("data/lipid20201020/IS_quantification/NEG/")
library(tidyverse)

trans_is_table(is_table_name = "IS_information_20201020.xlsx", polarity = "negative")

extract_targeted_peaks(
  path = ".",
  targeted_targeted_peak_table_name = "is_table.xlsx",
  sample_numer_show = 3,
  fit.gaussian = TRUE,
  integrate_xcms = FALSE
)
