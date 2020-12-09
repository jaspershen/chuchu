##no source
no_source()

###positive mode
sxtTools::setwd_project()
source("R/20201111/step0_parameter_setting.R")
lipid_data_um <-
  readr::read_csv(file.path("data", file_path, "absolute_quantification/Result/lipid_data_um.csv"))

setwd(file.path("data", file_path, "lipid_search_quantification/POS/"))

match_item <-
  match_item_pos

remove_idx <-
  lapply(match_item, function(x) {
    all(is.na(unique(x)))
  }) %>%
  unlist() %>%
  which()

if (length(remove_idx) > 0) {
  match_item <-
    match_item[-remove_idx]
}

for (i in 1:length(match_item)) {
  cat(i, "\n")
  dir.create(file.path("peak_shape/", names(match_item)[i]), showWarnings = FALSE)
  peak_name <-
    lipid_data_um %>%
    dplyr::filter(Class == names(match_item)[i]) %>%
    dplyr::pull(peak_name)
  
  lipid_name <-
    lipid_data_um %>%
    dplyr::filter(Class == names(match_item)[i]) %>%
    dplyr::pull(name) %>%
    stringr::str_replace_all("\\/", "_")
  
  if (length(peak_name) > 0) {
    file_name1 <- paste(peak_name, ".html", sep = "")

    file.copy(
      from = file.path("peak_shape", file_name1),
      to = file.path("peak_shape", names(match_item)[i]),
      overwrite = TRUE, recursive = TRUE
    )

    file.rename(
      from = file.path("peak_shape", names(match_item)[i], file_name1),
      to = file.path(
        "peak_shape",
        names(match_item)[i],
        paste(lipid_name, ".html", sep = "")
      )
    )
    
    unlink(file.path("peak_shape", file_name1), recursive = TRUE, force = TRUE)
  }
}



###negative mode----------------------------------------------------------------
sxtTools::setwd_project()
source("R/20201111/step0_parameter_setting.R")
lipid_data_um <-
  readr::read_csv(file.path("data", file_path, "absolute_quantification/Result/lipid_data_um.csv"))

setwd(file.path("data", file_path, "lipid_search_quantification/NEG/"))

match_item <-
  match_item_neg

remove_idx <-
  lapply(match_item, function(x) {
    all(is.na(unique(x)))
  }) %>%
  unlist() %>%
  which()

if (length(remove_idx) > 0) {
  match_item <-
    match_item[-remove_idx]
}

for (i in 1:length(match_item)) {
  cat(i, "\n")
  dir.create(file.path("peak_shape/", names(match_item)[i]), showWarnings = FALSE)
  peak_name <-
    lipid_data_um %>%
    dplyr::filter(Class == names(match_item)[i]) %>%
    dplyr::pull(peak_name)
  
  lipid_name <-
    lipid_data_um %>%
    dplyr::filter(Class == names(match_item)[i]) %>%
    dplyr::pull(name) %>%
    stringr::str_replace_all("\\/", "_")
  
  if (length(peak_name) > 0) {
    file_name1 <- paste(peak_name, ".html", sep = "")
    
    file.copy(
      from = file.path("peak_shape", file_name1),
      to = file.path("peak_shape", names(match_item)[i]),
      overwrite = TRUE,
      recursive = TRUE
    )

    file.rename(
      from = file.path("peak_shape", names(match_item)[i], file_name1),
      to = file.path(
        "peak_shape",
        names(match_item)[i],
        paste(lipid_name, ".html", sep = "")
      )
    )
    
    unlink(file.path("peak_shape", file_name1),
           recursive = TRUE,
           force = TRUE)
    
  }
}
