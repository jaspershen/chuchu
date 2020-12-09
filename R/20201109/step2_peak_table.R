##
no_source()

sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")
source("R/data_preparation.R")
library(tidyverse)

###specific folder and file
file_path <- "lipid20201020"
is_table_name <- "IS_information_20201020.xlsx"

setwd(file.path("data", file_path, "internal_standard"))
is_table <- readxl::read_xlsx(is_table_name)

is_table_pos <-
  clean_is_table(x = is_table, polarity = "positive")

sxtTools::setwd_project()
setwd(file.path("data", file_path, "peak_table/POS"))

######samples in different files
file <- dir(".")
remove_idx <- grep("pdf", file)
if(length(remove_idx) > 0) {
  file <- file[-remove_idx]
}

for (i in file) {
  cat(i, " ")
  peak_table_pos <-
    readr::read_csv(file.path(i, "Peak_table_for_cleaning.csv"))
  variable_info_pos <-
    peak_table_pos %>%
    dplyr::select(c(name, mz, rt))
  
  expression_data_pos <-
    peak_table_pos %>%
    dplyr::select(-c(name, mz, rt))
  
  ###missing value imputation
  library(impute)
  expression_data_pos <-
    impute.knn(
      data = as.matrix(expression_data_pos),
      k = 10,
      rowmax = 0.9,
      colmax = 0.9
    )
  
  expression_data_pos <-
    expression_data_pos$data
  
  peak_table_pos <-
    data.frame(variable_info_pos, expression_data_pos, stringsAsFactors = FALSE)
  
  ####get IS from peak table
  setwd(i)
  get_is_quantify_table(
    is_table = is_table_pos,
    peak_table = list(peak_table_pos),
    # sample_info = sample_info,
    mz.tol = 25,
    rt.tol = 90,
    figure = TRUE, 
    polarity = "positive"
  )
  
  setwd("..")
  write.csv(peak_table_pos,
            file.path(i, "peak_table_pos.csv"),
            row.names = FALSE)
}


###-----------------------------------------------------------------------------
####negative mode
sxtTools::setwd_project()
source("R/tools.R")
library(tidyverse)
setwd(file.path("data", file_path, "internal_standard"))
is_table <- readxl::read_xlsx(is_table_name)

is_table_neg <-
  clean_is_table(x = is_table, polarity = "negative")

is_table_neg

sxtTools::setwd_project()
setwd(file.path("data", file_path, "peak_table/NEG"))

######samples in different files
file <- dir(".")

remove_idx <- 
  grep("PDF", file)

if (length(remove_idx) > 0) {
  file <- file[-remove_idx]
}

for (i in file) {
  cat(i, " ")
  peak_table_neg <-
    readr::read_csv(file.path(i, "Peak_table_for_cleaning.csv"))
  
  variable_info_neg <-
    peak_table_neg %>%
    dplyr::select(c(name, mz, rt))
  
  expression_data_neg <-
    peak_table_neg %>%
    dplyr::select(-c(name, mz, rt))
  
  ###missing value imputation
  library(impute)
  expression_data_neg <-
    impute.knn(
      data = as.matrix(expression_data_neg),
      k = 10,
      rowmax = 0.9,
      colmax = 0.9
    )
  
  expression_data_neg <-
    expression_data_neg$data
  
  peak_table_neg <-
    data.frame(variable_info_neg, expression_data_neg, stringsAsFactors = FALSE)
  
  ####get IS from peak table
  setwd(i)
  get_is_quantify_table(
    is_table = is_table_neg,
    peak_table = list(peak_table_neg),
    # sample_info = sample_info,
    mz.tol = 25,
    rt.tol = 90,
    figure = TRUE, 
    polarity = "negative"
  )
  setwd("..")
  write.csv(peak_table_neg,
            file.path(i, "peak_table_neg.csv"),
            row.names = FALSE)
}

