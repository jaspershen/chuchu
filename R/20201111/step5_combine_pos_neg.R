##no source
no_source()

sxtTools::setwd_project()
source("R/20201111/step0_parameter_setting.R")
setwd(file.path("data", file_path, "absolute_quantification"))

###combine pos and neg
load("POS/variable_info_abs")
variable_info_abs_pos <- variable_info_abs

load("POS/express_data_abs_ug_ml")
express_data_abs_ug_ml_pos <-
  express_data_abs_ug_ml

load("POS/express_data_abs_um")
express_data_abs_um_pos <-
  express_data_abs_um

load("NEG/variable_info_abs")
variable_info_abs_neg <- variable_info_abs

load("NEG/express_data_abs_ug_ml")
express_data_abs_ug_ml_neg <-
  express_data_abs_ug_ml

load("NEG/express_data_abs_um")
express_data_abs_um_neg <-
  express_data_abs_um

colnames(variable_info_abs_pos)
colnames(variable_info_abs_neg)

intersect_name <-
  intersect(colnames(express_data_abs_um_pos),
            colnames(express_data_abs_um_neg))

express_data_abs_ug_ml_pos <-
  express_data_abs_ug_ml_pos[, intersect_name]

express_data_abs_um_pos <-
  express_data_abs_um_pos[, intersect_name]

express_data_abs_ug_ml_neg <-
  express_data_abs_ug_ml_neg[, intersect_name]

express_data_abs_um_neg <-
  express_data_abs_um_neg[, intersect_name]

colnames(express_data_abs_um_pos)
colnames(express_data_abs_um_neg)

dim(express_data_abs_um_pos)
dim(express_data_abs_um_neg)

combine_pos_neg_quantification(
  path = "Result",
  express_data_abs_ug_ml_pos = express_data_abs_ug_ml_pos,
  express_data_abs_um_pos = express_data_abs_um_pos,
  variable_info_abs_pos = variable_info_abs_pos,
  express_data_abs_ug_ml_neg = express_data_abs_ug_ml_neg,
  express_data_abs_um_neg = express_data_abs_um_neg,
  variable_info_abs_neg = variable_info_abs_neg
)

