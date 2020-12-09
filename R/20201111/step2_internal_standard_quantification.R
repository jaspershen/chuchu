##no source
no_source()

sxtTools::setwd_project()

source("R/20201111/step0_parameter_setting.R")

file.copy(
  from = file.path("data", file_path, "internal_standard", is_table_name),
  to = file.path("data", file_path, "IS_quantification/POS"),
  overwrite = TRUE,
  recursive = TRUE
)

setwd(file.path("data", file_path, "IS_quantification/POS"))

###positive
trans_is_table(is_table_name = is_table_name,
               polarity = "positive")

extract_targeted_peaks(
  path = ".",
  targeted_targeted_peak_table_name = "is_table.xlsx",
  forced_targeted_peak_table_name = NULL,
  fit.gaussian = TRUE,
  integrate_xcms = TRUE,
  output_eic = TRUE,
  output_integrate = TRUE,
  ppm = 40,
  rt.tolerance = 180
)

###stop
###manual check peak integration
extract_targeted_peaks(
  path = ".",
  targeted_targeted_peak_table_name = "is_table.xlsx",
  forced_targeted_peak_table_name = "forced_targeted.xlsx",
  fit.gaussian = TRUE,
  integrate_xcms = TRUE,
  output_eic = TRUE,
  ppm = 40,
  rt.tolerance = 180,
  output_integrate = TRUE
)

##combine raw lipid table and our quantification table
quantification_data <-
  readxl::read_xlsx("quantification_table.xlsx")

data_pos <- readxl::read_xlsx(is_table_name) %>%
  dplyr::rename(name = `Compound Name`) %>%
  dplyr::mutate(name = stringr::str_trim(name))

##impute MV
quantification_data$Blank1[is.na(quantification_data$Blank1)] <- 0

quantification_data <-
  process_mv(quantification_data = quantification_data)

quantification_data <-
  quantification_data %>%
  dplyr::left_join(data_pos, by = "name") %>%
  dplyr::select(name:adduct, Index:`RT NEG (second)`, everything())

openxlsx::write.xlsx(quantification_data,
                     file = "quantification_data_final.xlsx",
                     asTable = TRUE)


###negative
sxtTools::setwd_project()
source("R/20201111/step0_parameter_setting.R")

file.copy(
  from = file.path("data", file_path, "internal_standard", is_table_name),
  to = file.path("data", file_path, "IS_quantification/NEG"),
  overwrite = TRUE,
  recursive = TRUE
)

setwd(file.path("data", file_path, "IS_quantification/NEG"))
library(tidyverse)

trans_is_table(is_table_name = is_table_name, polarity = "negative")

extract_targeted_peaks(
  path = ".",
  targeted_targeted_peak_table_name = "is_table.xlsx",
  forced_targeted_peak_table_name = NULL,
  fit.gaussian = TRUE,
  integrate_xcms = TRUE,
  output_eic = TRUE,
  output_integrate = TRUE,
  ppm = 40,
  rt.tolerance = 180
)

###manual check
###manual check peak integration
# extract_targeted_peaks(
#   path = ".",
#   targeted_targeted_peak_table_name = "is_table.xlsx",
#   forced_targeted_peak_table_name = "forced_targeted.xlsx",
#   fit.gaussian = TRUE,
#   integrate_xcms = TRUE,
#   output_eic = TRUE,
#   ppm = 40,
#   rt.tolerance = 180,
#   output_integrate = TRUE
# )

##combine raw lipid table and our quantification table
quantification_data <-
  readxl::read_xlsx("quantification_table.xlsx")

data_neg <- readxl::read_xlsx(is_table_name) %>%
  dplyr::rename(name = `Compound Name`) %>%
  dplyr::mutate(name = stringr::str_trim(name))

##impute MV
quantification_data$Blank1[is.na(quantification_data$Blank1)] <- 0
quantification_data <-
  process_mv(quantification_data)

quantification_data %>%
  dplyr::select(name:adduct) %>%
  is.na() %>%
  sum()

quantification_data %>%
  dplyr::select(name:adduct) %>%
  `<`(0) %>%
  sum()

quantification_data <-
  quantification_data %>%
  dplyr::left_join(data_neg, by = "name") %>%
  dplyr::select(name:adduct, Index:`RT NEG (second)`, everything())

openxlsx::write.xlsx(quantification_data,
                     file = "quantification_data_final.xlsx",
                     asTable = TRUE)
