##no source
no_source()

sxtTools::setwd_project()

source("R/20201111/step0_parameter_setting.R")

file.copy(
  from = file.path("data", file_path, "internal_standard", is_table_name),
  to = file.path("data", file_path, "lipid_search_quantification/POS"),
  overwrite = TRUE,
  recursive = TRUE
)

##positive
dir.create(file.path("data", file_path, "lipid_search_quantification/POS"))
setwd(file.path("data", file_path, "lipid_search_quantification/POS"))

###Copy mzXML data from IS_quantification folder
from_file <- dir("../../IS_quantification/POS", full.names = TRUE)
from_file <- grep("mzXML" , from_file, value = TRUE)
file.copy(
  from = from_file,
  to = ".",
  overwrite = TRUE,
  recursive = TRUE
)

data_pos <-
  tidy_lipidsearch_data(file = lipid_data_pos_name,
                        path = ".",
                        polarity = "positive")

data_pos <-
  clean_lipid_data(x = data_pos)

###add Chol
chol_matrix <-
  data.frame(
    peak_name = "peak_pos_chol",
    mz = 369.3521,
    rt = chol_rt / 60,
    name = "Cholesterol",
    lipid_raw_name = NA,
    adduct = "NH4",
    polarity = "positive",
    Class = "Chol",
    FattyAcid = NA,
    IonFormula = NA,
    FA1 = NA,
    FA2 = NA,
    FA3 = NA,
    FA4 = NA,
    mean.int = 100000
  )

sample_matrix <-
  matrix(NA,
         nrow = 1,
         ncol = data_pos %>% dplyr::select(-c(peak_name:mean.int)) %>% ncol()) %>%
  as.data.frame()

colnames(sample_matrix) <-
  data_pos %>% dplyr::select(-c(peak_name:mean.int)) %>% colnames()

data_pos <-
  rbind(data_pos,
        cbind(chol_matrix, sample_matrix))

peak_table <-
  data_pos %>%
  dplyr::select(name = peak_name, mz, rt, adduct, Class, compound_name = name) %>%
  dplyr::mutate(rt = rt * 60)

openxlsx::write.xlsx(peak_table, "peak_table.xlsx", asTable = TRUE)
dim(peak_table)

extract_targeted_peaks(
  path = ".",
  targeted_targeted_peak_table_name = "peak_table.xlsx",
  forced_targeted_peak_table_name = NULL,
  fit.gaussian = TRUE,
  integrate_xcms = TRUE,
  output_eic = TRUE,
  output_integrate = TRUE,
  ppm = 40,
  rt.tolerance = 180,
  from_lipid_search = TRUE
)

###manual checking
extract_targeted_peaks(
  path = ".",
  targeted_targeted_peak_table_name = "peak_table.xlsx",
  forced_targeted_peak_table_name = "forced_targeted.xlsx",
  fit.gaussian = TRUE,
  integrate_xcms = TRUE,
  output_eic = TRUE,
  ppm = 40,
  rt.tolerance = 180,
  output_integrate = TRUE, 
  from_lipid_search = TRUE
)

##combine lipid search and our quantification table
quantification_data <-
  readxl::read_xlsx("quantification_table.xlsx")

##impute MV
###remove peaks with more than 80% missing values
quantification_data$Blank1[is.na(quantification_data$Blank1)] <- 0
quantification_data <-
  process_mv(quantification_data, na.tolerance = 0.5)

colnames(quantification_data)[1] <- "peak_name"

quantification_data <-
  quantification_data %>%
  dplyr::select(-c(mz, rt, adduct)) %>%
  dplyr::left_join(data_pos %>%
                     dplyr::select(peak_name:mean.int), by = "peak_name") %>%
  dplyr::select(peak_name, mz:mean.int, everything())

quantification_data$mean.int <-
  quantification_data %>%
  dplyr::select(-c(peak_name:mean.int)) %>%
  apply(1, function(x) {
    mean(x, na.rm = TRUE)
  })

openxlsx::write.xlsx(quantification_data,
                     file = "quantification_data_final.xlsx",
                     asTable = TRUE)

##delete mzXML data
file_mzxml <- c(grep("mzXML",dir("."), value = TRUE), "peak_data", "raw_data")

unlink(x = file_mzxml, recursive = TRUE, force = TRUE)





##negative---------------------------------------------------------------------
sxtTools::setwd_project()
file.copy(
  from = file.path("data", file_path, "internal_standard", is_table_name),
  to = file.path("data", file_path, "lipid_search_quantification/NEG"),
  overwrite = TRUE,
  recursive = TRUE
)
setwd(file.path("data", file_path, "lipid_search_quantification/NEG"))
library(tidyverse)

###Copy mzXML data from IS_quantification folder
from_file <- dir("../../IS_quantification/NEG", full.names = TRUE)
from_file <- grep("mzXML" , from_file, value = TRUE)
file.copy(
  from = from_file,
  to = ".",
  overwrite = TRUE,
  recursive = TRUE
)

data_neg <-
  tidy_lipidsearch_data(file = lipid_data_neg_name,
                        path = ".",
                        polarity = "negative")

data_neg <-
  clean_lipid_data(x = data_neg)

peak_table <-
  data_neg %>%
  dplyr::select(name = peak_name, mz, rt, adduct, Class, compound_name = name) %>%
  dplyr::mutate(rt = rt * 60)

openxlsx::write.xlsx(peak_table, "peak_table.xlsx", asTable = TRUE)

extract_targeted_peaks(
  path = ".",
  targeted_targeted_peak_table_name = "peak_table.xlsx",
  forced_targeted_peak_table_name = NULL,
  from_lipid_search = TRUE,
  fit.gaussian = FALSE,
  integrate_xcms = TRUE,
  output_eic = TRUE,
  output_integrate = TRUE,
  ppm = 40,
  rt.tolerance = 180
)

###manual checking
extract_targeted_peaks(
  path = ".",
  targeted_targeted_peak_table_name = "peak_table.xlsx",
  forced_targeted_peak_table_name = "forced_targeted.xlsx",
  fit.gaussian = TRUE,
  integrate_xcms = TRUE,
  output_eic = TRUE,
  ppm = 40,
  rt.tolerance = 180,
  output_integrate = TRUE, 
  from_lipid_search = TRUE
)

##combine lipid search and our quantification table
quantification_data <-
  readxl::read_xlsx("quantification_table.xlsx")

##impute MV
quantification_data$Blank1[is.na(quantification_data$Blank1)] <- 0
quantification_data <-
  process_mv(quantification_data, na.tolerance = 0.5)

colnames(quantification_data)[1] <- "peak_name"

quantification_data <-
  quantification_data %>%
  dplyr::select(-c(mz, rt, adduct)) %>%
  dplyr::left_join(data_neg %>%
                     dplyr::select(peak_name:mean.int), by = "peak_name") %>%
  dplyr::select(peak_name, mz:mean.int, everything())

quantification_data$mean.int <-
  quantification_data %>%
  dplyr::select(-c(peak_name:mean.int)) %>%
  apply(1, function(x) {
    mean(x, na.rm = TRUE)
  })

openxlsx::write.xlsx(quantification_data,
                     file = "quantification_data_final.xlsx",
                     asTable = TRUE)

##delete mzXML data
file_mzxml <- c(grep("mzXML",dir("."), value = TRUE), "peak_data", "raw_data")
unlink(x = file_mzxml, recursive = TRUE, force = TRUE)
