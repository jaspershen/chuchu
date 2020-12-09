sxtTools::setwd_project()
rm(list = ls())
source("R/20201111/tools.R")

library(xcms)
library(MSnbase)
library(tidyverse)
###where is the data
file_path <- "lipid20201025"
##the name of the IS table
is_table_name = "IS_information_20201025.xlsx"

lipid_data_pos_name = "102820_Align_Pos.xlsx"
lipid_data_neg_name = "102820_Align_Neg.xlsx"

chol_rt <- 1169

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
    "ChE" = c("18:1(d7) Chol Ester", "Cholesterol (d7)"),
    # "CL" = NA,
    # "Co" = NA,
    # "DG" = "15:0-18:1(d7) DAG",
    "Hex1Cer" = NA,
    # "Hex2Cer" = NA,
    "LPC" = "18:1(d7) Lyso PC",
    "LPE" = "18:1(d7) Lyso PE",
    "LPI" = NA,
    "MG" = NA,
    "PA" = NA,
    "PC" = "15:0-18:1(d7) PC",
    "PE" = "15:0-18:1(d7) PE",
    "PG" = "15:0-18:1(d7) PG (Na Salt)",
    "PI" = "15:0-18:1(d7) PI (NH4 Salt)",
    "PPE" = "C18(Plasm)-18:1(d9) PE",
    "PS" = "15:0-18:1(d7) PS (Na Salt)",
    "SM" = "d18:1-18:1(d9) SM"
    # "ST" = NA,
    # "SPH" = NA,
    # "SPHP" = NA,
    # "TG" = "15:0-18:1(d7)-15:0 TAG",
    # "ZyE" = NA
  )

###check files
###internal standard

if (all(dir("data") != file.path(file_path))) {
  stop("No ", file_path, " folder\n")
} else{
  ####internal_standard
  if (all(dir(file.path("data", file_path)) != "internal_standard")) {
    stop("No internal_standard folder\n")
  } else{
    ##is_table_name
    if (all(dir(file.path(
      "data", file_path, "internal_standard"
    ))
    != is_table_name)) {
      stop("No ",
           is_table_name,
           " in data/",
           file_path,
           "/internal_standard")
    }
    
    if (all(dir(file.path(
      "data", file_path, "internal_standard"
    ))
    != "mzXML")) {
      stop("No mzXML folder in data/",
           file_path,
           "/internal_standard")
    }
  }
  
  ##IS_quantification
  if (all(dir(file.path("data", file_path)) != "IS_quantification")) {
    stop("No IS_quantification folder\n")
  } else{
    if (all(dir(file.path(
      "data", file_path, "IS_quantification"
    )) != "POS")) {
      stop("No POS folder in IS_quantification")
    }
    
    if (all(dir(file.path(
      "data", file_path, "IS_quantification"
    )) != "NEG")) {
      stop("No POS folder in IS_quantification")
    }
  }
  
  
  ##lipid_search_quantification
  if (all(dir(file.path("data", file_path)) != "lipid_search_quantification")) {
    stop("No lipid_search_quantification folder\n")
  } else{
    if (all(dir(
      file.path("data", file_path, "lipid_search_quantification")
    ) != "POS")) {
      stop("No POS folder in lipid_search_quantification")
    } else{
      if (all(dir(
        file.path("data", file_path, "lipid_search_quantification/POS")
      ) != lipid_data_pos_name)) {
        stop("No ",
             lipid_data_pos_name,
             " in lipid_search_quantification/POS")
      }
      
      if (all(dir(
        file.path("data", file_path, "lipid_search_quantification/NEG")
      ) != lipid_data_neg_name)) {
        stop("No ",
             lipid_data_neg_name,
             " in lipid_search_quantification/NEG")
      }
      
    }
    
    if (all(dir(
      file.path("data", file_path, "lipid_search_quantification")
    ) != "NEG")) {
      stop("No POS folder in lipid_search_quantification")
    }
  }
  
  ###absolute_quantification
  if (all(dir(file.path("data", file_path)) != "absolute_quantification")) {
    stop("No absolute_quantification folder\n")
  } else{
    if (all(dir(file.path(
      "data", file_path, "absolute_quantification"
    )) != "POS")) {
      stop("No POS folder in absolute_quantification")
    }
    
    if (all(dir(file.path(
      "data", file_path, "absolute_quantification"
    )) != "NEG")) {
      stop("No POS folder in absolute_quantification")
    }
  }
}

