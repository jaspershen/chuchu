sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")

# load data
load("data/lipid20200727/lipidsearch_result/lipid_table")

is_table_pos <- 
  readr::read_csv("data/lipid20200727/internal_standard_table/is_table_pos.csv")

is_table_neg <- 
  readr::read_csv("data/lipid20200727/internal_standard_table/is_table_neg.csv")

setwd("data/lipid20200727/data_preparation_for_analysis/")

library(tidyverse)



##first, we should find the is in the peak table
data1 <- is_table_pos[,c("mz", "rt")]

data2 <- peak_table_pos[,c("mz", "rt")]

match_result_pos <- 
  sxtTools::sxtMTmatch(data1 = as.matrix(data1), 
                       data2 = as.matrix(data2), 
                       mz.tol = 25,
                       rt.tol = 30,
                       rt.error.type = "abs")


data1 <- is_table_neg[,c("mz", "rt")]

data2 <- peak_table_neg[,c("mz", "rt")]

match_result_neg <- 
  sxtTools::sxtMTmatch(data1 = as.matrix(data1), 
                       data2 = as.matrix(data2), 
                       mz.tol = 25,
                       rt.tol = 30,
                       rt.error.type = "abs")


colnames(lipid_table_pos)
colnames(lipid_table_neg)

is_tag_pos <- 
is_table_pos[match_result_pos[,1],]

is_data_pos <- peak_table_pos[match_result_pos[,2],] %>% 
  dplyr::select(contains('X'))

is_data_pos <- 
  cbind(is_tag_pos, is_data_pos)

write.csv(is_data_pos, "is_data_pos.csv", row.names = FALSE)


is_tag_neg <- 
  is_table_neg[match_result_neg[,1],]

is_data_neg <- peak_table_neg[match_result_neg[,2],] %>% 
  dplyr::select(contains('X'))

is_data_neg <- 
  cbind(is_tag_neg, is_data_neg)

write.csv(is_data_neg, "is_data_neg.csv", row.names = FALSE)

save(is_data_pos, file = "is_data_pos")
save(is_data_neg, file = "is_data_neg")



##output relative
expression_data_relative <- 
  rbind(lipid_table_pos, lipid_table_neg)

variable_info_relative <- 
  expression_data_relative[,c(1:16)]

expression_data_relative <- 
  expression_data_relative[,-c(1:16)]


expression_data_relative <- 
  expression_data_relative %>% 
  dplyr::select(-contains("QC")) %>% 
  dplyr::select(-contains("Blk"))

variable_info_relative <- 
  variable_info_relative %>% 
  dplyr::select(peak_name, name, everything())

save(expression_data_relative, file = "expression_data_relative")
save(variable_info_relative, file = "variable_info_relative")

write.csv(expression_data_relative, "expression_data_relative.csv", row.names = FALSE)
write.csv(variable_info_relative, "variable_info_relative.csv", row.names = FALSE)


###calculate concentration for each samples and lipid
##positive
is_tag_pos <- 
  is_data_pos[,1:9]

is_table_pos <- 
  is_data_pos[,-c(1:9)] %>% 
  dplyr::select(contains("X"))

lipid_tag_pos <- 
  lipid_table_pos[,c(1:16)]

lipid_table_pos <- 
  lipid_table_pos[,-c(1:16)] %>% 
  dplyr::select(contains("X"))

colnames(lipid_table_pos) == colnames(is_table_pos)

lipid_tag_pos$Class %>% table()

is_tag_pos$name

match_item_pos <- 
list(
  "AEA" = NA,
  "Cer" = NA,
  "ChE" = c(9, 15),
  "CL" = NA,
  "DG" = 11,
  "HexlCer" = NA,
  "LPC" = 7,
  "LPE" = 8,
  "MG" = NA,
  "PC" = 1,
  "PE" = 2,
  "PI" = 5,
  "PS" = 3,
  "SM" = 13,
  "SPH" = NA,
  "SPHP" = NA,
  "ST" = NA,
  "TG" = 12
)


lipid_tag_pos$rt <- 
  lipid_tag_pos$rt * 60

expresion_data_abs_pos <- 
  cal_abs(lipid_tag = lipid_tag_pos, 
          lipid_table = lipid_table_pos, 
          is_tag = is_tag_pos, 
          is_table = is_table_pos, 
          match_item = match_item_pos)

express_data_abs_pos1 <- expresion_data_abs_pos$express_data_abs1
express_data_abs_pos2 <- expresion_data_abs_pos$express_data_abs2
variable_info_abs_pos <- expresion_data_abs_pos$variable_info_abs


##negative
is_tag_neg <- 
  is_data_neg[,1:9]

is_table_neg <- 
  is_data_neg[,-c(1:9)] %>% 
  dplyr::select(contains("X"))

lipid_tag_neg <- 
  lipid_table_neg[,c(1:16)]

lipid_table_neg <- 
  lipid_table_neg[,-c(1:16)] %>% 
  dplyr::select(contains("X"))

colnames(lipid_table_neg) == colnames(is_table_neg)

lipid_tag_neg$Class %>% table()

is_tag_neg$name

match_item_neg <- 
  list(
    "AEA" = NA,
    "Cer" = NA,
    "ChE" = NA,
    "CL" = NA,
    "DG" = NA,
    "HexlCer" = NA,
    "LPC" =NA,
    "LPE" = NA,
    "MG" = NA,
    "PC" = NA,
    "PE" = c(1,3),
    "PI" = NA,
    "PS" = NA,
    "SM" = NA,
    "SPH" = NA,
    "SPHP" = NA,
    "ST" = NA,
    "TG" = NA
  )


lipid_tag_neg$rt <- 
  lipid_tag_neg$rt * 60

expresion_data_abs_neg <- 
  cal_abs(lipid_tag = lipid_tag_neg, 
          lipid_table = lipid_table_neg, 
          is_tag = is_tag_neg, 
          is_table = is_table_neg, 
          match_item = match_item_neg)


express_data_abs_neg1 <- expresion_data_abs_neg$express_data_abs1
express_data_abs_neg2 <- expresion_data_abs_neg$express_data_abs2
variable_info_abs_neg <- expresion_data_abs_neg$variable_info_abs


###combine positive and negative
express_data_abs1 <- rbind(express_data_abs_pos1, express_data_abs_neg1)
express_data_abs2 <- rbind(express_data_abs_pos2, express_data_abs_neg2)
variable_info_abs <- rbind(variable_info_abs_pos, variable_info_abs_neg)

mass <- 
stringr::str_replace_all(variable_info_abs$IonFormula, " ", "") %>% 
  purrr::map(.f = function(x){
  x %>% 
      Rdisop::getMolecule() %>% 
      Rdisop::getMass()
  })
  

mass <- unlist(mass)

variable_info_abs <- data.frame(variable_info_abs, mass, stringsAsFactors = FALSE)

express_data_abs1 <- 
  express_data_abs1 %>% 
  as.matrix() %>% 
  `/`(mass)
  

save(express_data_abs1, file = "express_data_abs1")
save(express_data_abs2, file = "express_data_abs2")
save(variable_info_abs, file = "variable_info_abs")

write.csv(express_data_abs1, "expression_data_abs1.csv", row.names = FALSE)
write.csv(express_data_abs2, "expression_data_abs2.csv", row.names = FALSE)
write.csv(variable_info_abs, "variable_info_abs.csv", row.names = FALSE)






