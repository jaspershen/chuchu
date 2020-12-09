##no source
no_source()

sxtTools::setwd_project()
source("R/20201111/step0_parameter_setting.R")
setwd(file.path("data", file_path))

###positive mode
is_data_pos <- 
  readxl::read_xlsx("IS_quantification/POS/quantification_data_final.xlsx")

lipid_data_pos <- 
  readxl::read_xlsx("lipid_search_quantification/POS/quantification_data_final.xlsx")

get_quantification_data2(
  path = "absolute_quantification/POS/",
  lipid_data = lipid_data_pos,
  is_data = is_data_pos,
  match_item = match_item_pos
)

###negative mode
is_data_neg <- 
  readxl::read_xlsx("IS_quantification/NEG/quantification_data_final.xlsx")

lipid_data_neg <- 
  readxl::read_xlsx("lipid_search_quantification/NEG/quantification_data_final.xlsx")

get_quantification_data2(
  path = "absolute_quantification/NEG/",
  lipid_data = lipid_data_neg,
  is_data = is_data_neg,
  match_item = match_item_neg
)


