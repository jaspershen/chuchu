##avoid source
no_function()

sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")

file_path <- "lipid20201020"
setwd(file.path("data", file_path, "Result"))

####extract peaks for each lipid in some samples
lipid_data <- readr::read_csv("lipid_data_um.csv")


match_item_pos <-
  list(
    "AEA" = NA,
    "Cer" = "d18:1 (d7)-15:0 Cer",
    "ChE" = c("18:1(d7) Chol Ester", "Cholesterol (d7)"),
    "Chol" = "Cholesterol (d7)",
    "CL" = NA,
    "Co" = NA,
    "DG" = "15:0-18:1(d7) DAG",
    "Hex1Cer" = NA,
    "Hex2Cer" = NA,
    "LPC" = "18:1(d7) Lyso PC",
    "LPE" = "18:1(d7) Lyso PE",
    "MG" = "18:1 (d7) MG",
    "PA" = "15:0-18:1(d7) PA (Na Salt)",
    "PC" = "15:0-18:1(d7) PC",
    "PE" = "15:0-18:1(d7) PE",
    "PG" = "15:0-18:1(d7) PG (Na Salt)",
    "PI" = "15:0-18:1(d7) PI (NH4 Salt)",
    "PPE" = "C18(Plasm)-18:1(d9) PE",
    "PS" = "15:0-18:1(d7) PS (Na Salt)",
    "SM" = "d18:1-18:1(d9) SM",
    "ST" = NA,
    # "SPH" = NA,
    # "SPHP" = NA,
    "TG" = "15:0-18:1(d7)-15:0 TAG",
    "ZyE" = NA
  )

match_item_neg <-
  list(
    # "AEA" = NA,
    "Cer" = "d18:1 (d7)-15:0 Cer",
    "Chol" = "Cholesterol (d7)",
    # "ChE" = c("18:1(d7) Chol Ester", "Cholesterol (d7)"),
    # "CL" = NA,
    # "Co" = NA,
    # "DG" = "15:0-18:1(d7) DAG",
    "Hex1Cer" = NA,
    # "Hex2Cer" = NA,
    "LPC" = "18:1(d7) Lyso PC",
    "LPE" = "18:1(d7) Lyso PE",
    "LPI" = NA,
    # "MG" = NA,
    "PA" = NA,
    "PC" = "15:0-18:1(d7) PC",
    "PE" = "15:0-18:1(d7) PE",
    # "PG" = "15:0-18:1(d7) PG (Na Salt)",
    "PI" = "15:0-18:1(d7) PI (NH4 Salt)",
    # "PPE" = "C18(Plasm)-18:1(d9) PE",
    "PS" = "15:0-18:1(d7) PS (Na Salt)",
    "SM" = "d18:1-18:1(d9) SM"
    # "ST" = NA,
    # "SPH" = NA,
    # "SPHP" = NA,
    # "TG" = "15:0-18:1(d7)-15:0 TAG",
    # "ZyE" = NA
  )


###positive
extract_peak(path = ".", 
             peak_table = lipid_data %>% dplyr::filter(polarity == "positive"),
             match_item = match_item_pos,
             polarity = "positive")


###negative
extract_peak(path = ".", 
             peak_table = lipid_data %>% dplyr::filter(polarity == "negative"),
             match_item = match_item_neg,
             polarity = "negative")



