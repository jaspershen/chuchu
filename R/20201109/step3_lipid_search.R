##avoid source
no_function()

sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")

file_path <- "lipid20201020"

######set work directory
setwd(file.path("data", file_path))

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

#####get concentration
sample_info <-
  data.frame(sample_id = c("Blank",
                           paste("test_sample", 1:4, sep = "")),
             class = c(rep("test", 5)))

sample_info

path = "."
sample_info = sample_info
lipid_data_pos_name = "102420_Align_Pos.xlsx"
lipid_data_neg_name = "102420_Align_Neg.xlsx"
is_table_name =
  "IS_information_20201020.xlsx"
##no change
peak_table_pos_name =
  "peak_table_pos.csv"
peak_table_neg_name =
  "peak_table_neg.csv"

match_item_pos = match_item_pos
match_item_neg = match_item_neg

get_quantification_data(
  path = ".",
  sample_info = sample_info,
  lipid_data_pos_name = lipid_data_pos_name,
  lipid_data_neg_name = lipid_data_neg_name,
  is_table_name = is_table_name,
  peak_table_pos_name = peak_table_pos_name,
  peak_table_neg_name = peak_table_neg_name,
  match_item_pos = match_item_pos,
  match_item_neg = match_item_neg
)



