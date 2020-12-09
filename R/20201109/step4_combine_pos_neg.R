##avoid source
no_function()

sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")

file_path <- "lipid20201020"
setwd(file.path("data", file_path))

###combine pos and neg
load("lipid_search/POS/variable_info_abs_pos")
load("lipid_search/POS/express_data_abs_pos1")
load("lipid_search/POS/express_data_abs_pos2")

load("lipid_search/NEG/variable_info_abs_neg")
load("lipid_search/NEG/express_data_abs_neg1")
load("lipid_search/NEG/express_data_abs_neg2")

colnames(variable_info_abs_pos)
colnames(variable_info_abs_neg)

intersect_name <-
  intersect(colnames(express_data_abs_pos1),
            colnames(express_data_abs_neg1))

express_data_abs_pos1 <-
  express_data_abs_pos1[, intersect_name]

express_data_abs_pos2 <-
  express_data_abs_pos2[, intersect_name]

express_data_abs_neg1 <-
  express_data_abs_neg1[, intersect_name]

express_data_abs_neg2 <-
  express_data_abs_neg2[, intersect_name]

colnames(express_data_abs_pos1)
colnames(express_data_abs_neg1)

colnames(express_data_abs_pos2)
colnames(express_data_abs_neg2)

dim(express_data_abs_pos1)
dim(express_data_abs_pos2)
dim(express_data_abs_neg1)
dim(express_data_abs_neg2)

combine_pos_neg_quantification(
  path = "Result",
  express_data_abs_pos1 = express_data_abs_pos1,
  express_data_abs_pos2 = express_data_abs_pos2,
  variable_info_abs_pos = variable_info_abs_pos,
  express_data_abs_neg1 = express_data_abs_neg1,
  express_data_abs_neg2 = express_data_abs_neg2,
  variable_info_abs_neg = variable_info_abs_neg
)


