# lipid_tag <- lipid_tag_pos
# lipid_table <- lipid_table_pos
# is_tag <- is_tag_pos
# is_table <- is_table_pos
# match_item <- match_item_pos
cal_abs <- function(lipid_tag,
                    lipid_table,
                    is_tag,
                    is_table,
                    match_item) {
  ##remove the class which can bot be use for absolute quantificaion
  remain_idx <- lapply(match_item, function(x) {
    !is.na(x)[1]
  }) %>%
    unlist() %>%
    which()
  
  match_item <- match_item[remain_idx]
  
  ####get the final match_item and final is_table and is_table
  intersect_name <-
    intersect(is_tag$name, unique(unlist(match_item)))
  
  remain_idx <- which(is_tag$name %in% intersect_name)
  
  is_tag <-
    is_tag[remain_idx, ]
  
  is_table <-
    is_table[remain_idx, ]
  
  match_item <-
    lapply(match_item, function(x) {
      x[(x %in% intersect_name)]
    })
  
  remain_idx <-
    lapply(match_item, length) %>%
    unlist() %>%
    `!=`(0)
  
  match_item <- match_item[remain_idx]
  
  ########--------------------------------------
  remain_idx <-
    which(lipid_tag$Class %in% names(match_item))
  
  lipid_table <- lipid_table[remain_idx,]
  lipid_tag <- lipid_tag[remain_idx, ]
  
  express_data_abs_ug_ml <-
    express_data_abs_um <-
    lipid_table %>% 
    as.data.frame()
  
  for (i in 1:nrow(lipid_tag)) {
    cat(i, " ")
    temp_class <- lipid_tag$Class[i]
    temp_idx <- match_item[[temp_class]] %>%
      match(., is_tag$name)
    temp_rt <- lipid_tag$rt[i]
    
    temp_idx <-
      temp_idx[which.min(abs(temp_rt - is_tag$rt[temp_idx]))]
    
    ratio <-
      unlist(lipid_table[i, , drop = TRUE]) / unlist(is_table[temp_idx, , drop = TRUE])
    
    # ratio[is.na(ratio)] <- 0
    # ratio[is.infinite(ratio)] <- max(ratio[!is.infinite(ratio)])
    ratio[is.infinite(ratio)] <- 0
    ratio[is.na(ratio)] <- 0
    
    concentration1 <-
      is_tag$`μg/mL`[temp_idx] * ratio
    
    concentration2 <-
      is_tag$μM[temp_idx] * ratio
    
    express_data_abs_ug_ml[i, ] <- concentration1
    express_data_abs_um[i, ] <- concentration2
  }
  
  
  return_result <- list(
    express_data_abs_ug_ml = express_data_abs_ug_ml,
    express_data_abs_um = express_data_abs_um,
    variable_info_abs = lipid_tag
  )
  
  return(return_result)
  
}



###clean lipid data
clean_lipid_data <-
  function(x) {
    ##rename class
    x$Class[x$Class == "PC" &
              stringr::str_detect(x$name, "p")] <- "PPC"
    x$Class[x$Class == "PE" &
              stringr::str_detect(x$name, "p")] <- "PPE"
    
    ##remove some wrong annotation
    x <-
      x %>%
      dplyr::filter(rt > 1)
    
    dim(x)
    
    remain_idx <-
      as.data.frame(t(x)) %>%
      purrr::map(
        .f = function(x) {
          ##RT
          if (x[7] == "LPC" & as.numeric(x[4]) > 15) {
            return(FALSE)
          }
          
          if (x[7] == "LPE" & as.numeric(x[4]) > 15) {
            return(FALSE)
          }
          
          if (x[7] == "PC" & as.numeric(x[4]) < 15) {
            return(FALSE)
          }
          
          if (x[7] == "PE" & as.numeric(x[4]) < 15) {
            return(FALSE)
          }
          
          if (x[7] == "MG" & as.numeric(x[4]) > 10) {
            return(FALSE)
          }
          
          if (x[7] == "DG" & as.numeric(x[4]) < 10) {
            return(FALSE)
          }
          
          if (x[7] == "TG" & as.numeric(x[4]) < 15) {
            return(FALSE)
          }
          
          return(TRUE)
        }
      ) %>%
      unlist() %>%
      which()
    
    x <-
      x[remain_idx,]
    x
  }



clean_is_table <- function(x,
                           polarity = c("positive", "negative")) {
  polarity <- match.arg(polarity)
  
  x <-
    x %>%
    dplyr::arrange(Index) %>%
    dplyr::mutate(`Compound Name` = stringr::str_trim(`Compound Name`))
  
  if (any(colnames(x) == "RT POS (second)")) {
    x$`RT POS (min)` <- as.numeric(x$`RT POS (second)`) / 60
    x$`RT NEG (min)` <- as.numeric(x$`RT NEG (second)`) / 60
    x <-
      x %>%
      dplyr::select(-c(`RT POS (second)`, `RT NEG (second)`))
  }
  
  
  if (polarity == "positive") {
    x_pos <-
      x %>%
      dplyr::select(-c(`RT NEG (min)`, `Adduct NEG`)) %>%
      dplyr::mutate(`RT POS (min)` = as.numeric(`RT POS (min)`)) %>%
      dplyr::filter(!is.na(`RT POS (min)`)) %>%
      dplyr::rename(rt = `RT POS (min)`) %>%
      dplyr::mutate(rt = rt * 60)
    
    cat(nrow(x_pos),
        "out of",
        nrow(x),
        "internal standards have RT.\n")
    
    library(Rdisop)
    
    adduct <-
      x_pos$`Adduct POS`
    
    shift_mz <-
      purrr::map(
        .x = adduct,
        .f = function(x) {
          if (x == "M+H") {
            return(getMass(getMolecule(formula = "H")))
          }
          
          if (x == "M+H-H2O") {
            return(-getMass(getMolecule(formula = "HO")))
          }
          
          if (x == "M+NH4") {
            return(getMass(getMolecule(formula = "NH4")))
          }
          
          if (x == "M+NH4-H2O") {
            return(getMass(getMolecule(formula = "NH2")) - getMass(getMolecule(formula = "O")))
          }
          
          if (x == "M+Na") {
            return(getMass(getMolecule(formula = "Na")))
          }
          
          return(NA)
        }
      ) %>%
      unlist()
    
    mz <-
      purrr::map2(
        .x = as.numeric(x_pos$`Exact Mass`),
        .y = shift_mz,
        .f = function(x, y) {
          x + y
        }
      ) %>% unlist()
    
    x_pos$mz <- mz
    return(x_pos)
  }
  
  
  if (polarity == "negative") {
    x_neg <-
      x %>%
      dplyr::select(-c(`RT POS (min)`, `Adduct POS`)) %>%
      dplyr::mutate(`RT NEG (min)` = as.numeric(`RT NEG (min)`)) %>%
      dplyr::filter(!is.na(`RT NEG (min)`)) %>%
      dplyr::rename(rt = `RT NEG (min)`) %>%
      dplyr::mutate(rt = rt * 60)
    
    cat(nrow(x_neg),
        "out of",
        nrow(x),
        "internal standards have RT.\n")
    
    library(Rdisop)
    
    adduct <-
      x_neg$`Adduct NEG`
    
    shift_mz <-
      purrr::map(
        .x = adduct,
        .f = function(x) {
          if (x == "M-H") {
            return(-getMass(getMolecule(formula = "H")))
          }
          
          if (x == "M+CH3COO") {
            return(getMass(getMolecule(formula = "CH3COO")))
          }
          
          if (x == "M+HCOO") {
            return(getMass(getMolecule(formula = "HCOO")))
          }
          return(NA)
        }
      ) %>%
      unlist()
    
    mz <-
      purrr::map2(
        .x = as.numeric(x_neg$`Exact Mass`),
        .y = shift_mz,
        .f = function(x, y) {
          x + y
        }
      ) %>% unlist()
    
    x_neg$mz <- mz
    return(x_neg)
  }
  
}


match_sample_peak_table_lipid_search <- function(peak_table,
                                                 lipid_search) {
  colnames(peak_table) <-
    stringr::str_replace(colnames(peak_table), "^X", "")
  
  last_idx_pos <- which(colnames(lipid_search) == "mean.int")
  
  name1 <- setdiff(colnames(lipid_search)[-c(1:last_idx_pos)],
                   colnames(peak_table)[-c(1:3)])
  
  name2 <- setdiff(colnames(peak_table)[-c(1:3)],
                   colnames(lipid_search)[-c(1:last_idx_pos)])
  
  if (length(name1) > 0) {
    lipid_search <-
      lipid_search %>%
      dplyr::select(-name1)
  }
  
  if (length(name2) > 0) {
    peak_table <-
      peak_table %>%
      dplyr::select(-name2)
  }
  
  return(list(peak_table, lipid_search))
  
}




get_is_quantify_table <- function(is_table,
                                  peak_table,
                                  # sample_info,
                                  mz.tol = 25,
                                  rt.tol = 90,
                                  figure = FALSE,
                                  polarity = c("positive", "negative")) {
  polarity <- match.arg(polarity)
  data_is <- is_table[, c("mz", "rt")]
  
  is_data <-
    purrr::map(
      .x = peak_table,
      .f = function(x) {
        data_table <- x[, c("mz", "rt")]
        match_result <-
          sxtTools::sxtMTmatch(
            data1 = as.matrix(data_is),
            data2 = as.matrix(data_table),
            mz.tol = mz.tol,
            rt.tol = rt.tol,
            rt.error.type = "abs"
          )
        
        match_result <-
          match_result %>%
          as.data.frame() %>%
          dplyr::group_by(Index1) %>%
          dplyr::filter(`rt error` == min(`rt error`)) %>%
          dplyr::ungroup()
        
        cat(
          nrow(match_result),
          "out of",
          nrow(is_table),
          "internal standards are found in peak table.\n"
        )
        
        is_tag <-
          is_table[match_result$Index1,]
        
        is_data <- x[match_result$Index2, ]  %>%
          dplyr::select(-c(mz, rt))
        
        is_data <-
          cbind(is_tag, is_data)
        is_data
      }
    )
  
  is_data <-
    purrr::map2(
      .x = is_data,
      .y = 1:length(is_data),
      .f = function(x, y) {
        colnames(x)[11] <-
          paste("name", y, sep = "_")
        x
      }
    )
  
  if (length(is_data) == 1) {
    temp_data <- is_data[[1]]
  } else{
    temp_data <- is_data[[1]]
    for (i in 2:length(is_data)) {
      temp_data <-
        temp_data %>%
        dplyr::full_join(is_data[[i]],
                         by = intersect(colnames(temp_data), colnames(is_data[[i]])))
    }
  }
  
  write.csv(
    temp_data,
    file = paste(polarity, "is_data.csv", sep = "_"),
    row.names = FALSE
  )
  
  temp_data <-
    temp_data %>%
    dplyr::select(-contains("name_"))
  
  if (figure) {
    dir.create(path = "IS_figure", showWarnings = FALSE)
    name <- colnames(temp_data)[-c(1:10)]
    purrr::walk(
      as.data.frame(t(temp_data)),
      .f = function(x) {
        plot <-
          data.frame(name = name, value = as.numeric(x[-c(1:10)]))  %>%
          ggplot(aes(name, value)) +
          geom_point(shape = 21,
                     size = 4,
                     fill = "red") +
          labs(x = "", y = "Response") +
          theme_bw()
        ggsave(
          plot,
          filename = file.path("IS_figure", paste(x[2], ".pdf", sep = "")),
          width = 7,
          height = 7
        )
      }
    )
  }
  
  temp_data
  
}



get_quantification_data <- function(path = ".",
                                    sample_info,
                                    lipid_data_pos_name = "100720_Align_Pos.xlsx",
                                    lipid_data_neg_name = "100720_Align_Neg.xlsx",
                                    is_table_name =
                                      "IS_information_092920.xlsx",
                                    peak_table_pos_name =
                                      "peak_table_pos.csv",
                                    peak_table_neg_name =
                                      "peak_table_neg.csv",
                                    match_item_pos = match_item_pos,
                                    match_item_neg = match_item_neg) {
  ######################positive quantification-----------------------------------
  #######read lipid search data
  data_pos <-
    tidy_lipidsearch_data(
      file = file.path("lipid_search/POS", lipid_data_pos_name),
      path = ".",
      polarity = "positive"
    )
  
  data_pos <-
    clean_lipid_data(x = data_pos)
  
  ###remove samples which are not in sample_info
  last_idx <- which(colnames(data_pos) == "mean.int")
  data_pos <-
    data_pos[, c(c(1:last_idx), which(colnames(data_pos) %in% sample_info$sample_id))]
  
  ###read IS table
  is_table <-
    readxl::read_xlsx(file.path("internal_standard", is_table_name))
  
  is_table_pos <-
    clean_is_table(x = is_table, polarity = "positive")
  
  #####################read peak table-------------------------------------------
  peak_table_pos <-
    purrr::map(
      unique(sample_info$class),
      .f = function(x) {
        x <-
          readr::read_csv(file.path("peak_table/POS", x, peak_table_pos_name))
        colnames(x) <-
          stringr::str_replace(colnames(x), "^X", "")
        x
      }
    )
  
  
  ###match samples in lipid search and peak table---------------------------------
  #################find internal standards in peak table-------------------------
  is_data_pos <-
    get_is_quantify_table(
      is_table = is_table_pos,
      peak_table = peak_table_pos,
      mz.tol = 25,
      rt.tol = 90
    )
  
  ###get Cholesterol from peak table and add it to lipid search table
  temp_idx <-
    lapply(peak_table_pos, function(x) {
      which(abs(x$mz - 369.3521) * 10 ^ 6 / 369.3521 <= 25 &
              abs(x$rt - is_table_pos$rt[is_table_pos$`Compound Name` == "Cholesterol (d7)"]) <= 90)
    })
  
  if (length(unique(unlist(temp_idx))) >= 1) {
    temp_data <-
      purrr::map2(
        .x = temp_idx,
        .y = peak_table_pos,
        .f = function(x, y) {
          x <-
            x[which.min(abs(y$rt[x] - is_table_pos$rt[is_table_pos$`Compound Name` == "Cholesterol (d7)"]))]
          
          temp_data <-
            y[x, , drop = FALSE] %>%
            dplyr::mutate(peak_name = name) %>%
            dplyr::mutate(name = "Cholesterol") %>%
            dplyr::mutate(Class = "Chol") %>%
            dplyr::mutate(rt = rt / 60) %>%
            dplyr::mutate(mean.int = 10000)
          
          temp_data <-
            temp_data %>%
            dplyr::select(one_of(colnames(data_pos)))
          temp_data
        }
      )
    
    peak_name <-
      lapply(temp_data, function(x) {
        x$peak_name
      }) %>%
      unlist() %>%
      `[`(1)
    
    rt <-
      lapply(temp_data, function(x) {
        x$rt
      }) %>%
      unlist() %>%
      `[`(1)
    
    temp_data <-
      lapply(temp_data, function(x) {
        x$peak_name <- peak_name
        x$rt <- rt
        x
      })
    
    if (length(temp_data) == 1) {
      new_data <- temp_data[[1]]
    } else{
      new_data <- temp_data[[1]]
      for (x in temp_data[-1]) {
        new_data <-
          new_data %>%
          dplyr::left_join(x, by = intersect(colnames(x), colnames(new_data)))
      }
    }
    
    data_pos <-
      data_pos %>%
      dplyr::full_join(new_data, by = intersect(colnames(data_pos), colnames(new_data)))
    
    data_pos$adduct[is.na(data_pos$adduct)] <- "H-H2O"
    data_pos$polarity[is.na(data_pos$polarity)] <- "positive"
    data_pos$IonFormula[is.na(data_pos$IonFormula)] <- "C27 H45 O0"
  }
  
  ###calculate concentration for each samples and lipid--------------------------
  is_tag_pos <-
    is_data_pos[, 1:10]
  
  is_table_pos <-
    is_data_pos[, -c(1:10)]
  
  last_idx_pos <- which(colnames(data_pos) == "mean.int")
  
  lipid_tag_pos <-
    data_pos[, c(1:last_idx_pos)]
  
  lipid_table_pos <-
    data_pos[, -c(1:last_idx_pos)]
  
  lipid_table_pos <-
    lipid_table_pos[, sort(colnames(lipid_table_pos))]
  
  is_table_pos <-
    is_table_pos[, sort(colnames(is_table_pos))]
  
  lipid_tag_pos$rt <-
    lipid_tag_pos$rt * 60
  
  lipid_table_pos2 <-
    get_lipid_quantify_table(
      lipid_table = data.frame(
        peak_name = lipid_tag_pos$peak_name,
        mz = lipid_tag_pos$mz,
        rt = lipid_tag_pos$rt,
        stringsAsFactors = FALSE
      ),
      peak_table = peak_table_pos,
      mz.tol = 25,
      rt.tol = 90,
      figure = FALSE
    )
  
  lipid_table_pos2 <-
    lipid_tag_pos[, "peak_name", drop = FALSE] %>%
    dplyr::left_join(lipid_table_pos2 %>% dplyr::select(-c(mz, rt)), by = "peak_name") %>%
    dplyr::select(-peak_name)
  
  expresion_data_abs_pos <-
    cal_abs(
      lipid_tag = lipid_tag_pos,
      lipid_table = lipid_table_pos,
      is_tag = is_tag_pos,
      is_table = is_table_pos,
      match_item = match_item_pos
    )
  
  express_data_abs_ug_ml_pos <- expresion_data_abs_pos$express_data_abs_ug_ml
  express_data_abs_um_pos <- expresion_data_abs_pos$express_data_abs_um
  variable_info_abs_pos <- expresion_data_abs_pos$variable_info_abs
  
  save(variable_info_abs_pos, file = "lipid_search/POS/variable_info_abs_pos")
  save(express_data_abs_ug_ml_pos, file = "lipid_search/POS/express_data_abs_ug_ml_pos")
  save(express_data_abs_um_pos, file = "lipid_search/POS/express_data_abs_um_pos")
  save(lipid_tag_pos, file = "lipid_search/POS/lipid_tag_pos")
  save(lipid_table_pos, file = "lipid_search/POS/lipid_table_pos")
  save(is_tag_pos, file = "lipid_search/POS/is_tag_pos")
  save(is_table_pos, file = "lipid_search/POS/is_table_pos")
  
  
  ######################negative quantification-----------------------------------
  #######read lipid search data
  data_neg <-
    tidy_lipidsearch_data(
      file = file.path("lipid_search/NEG", lipid_data_neg_name),
      path = ".",
      polarity = "negative"
    )
  
  data_neg <-
    clean_lipid_data(x = data_neg)
  
  ###remove samples which are not in sample_info
  last_idx <- which(colnames(data_neg) == "mean.int")
  data_neg <-
    data_neg[, c(c(1:last_idx), which(colnames(data_neg) %in% sample_info$sample_id))]
  
  ###read IS table
  is_table <-
    readxl::read_xlsx(file.path("internal_standard", is_table_name))
  
  is_table_neg <-
    clean_is_table(x = is_table, polarity = "negative")
  
  #####################read peak table-------------------------------------------
  peak_table_neg <-
    purrr::map(
      unique(sample_info$class),
      .f = function(x) {
        x <-
          readr::read_csv(file.path("peak_table/NEG", x, peak_table_neg_name))
        colnames(x) <-
          stringr::str_replace(colnames(x), "^X", "")
        x
      }
    )
  
  
  ###match samples in lipid search and peak table---------------------------------
  #################find internal standards in peak table-------------------------
  is_data_neg <-
    get_is_quantify_table(is_table = is_table_neg,
                          peak_table = peak_table_neg)
  
  ###get ChE from peak table and add it to lipid search table
  temp_idx <-
    lapply(peak_table_neg, function(x) {
      which(abs(x$mz - 369.3521) * 10 ^ 6 / 369.3521 <= 25 &
              abs(x$rt - is_table_neg$rt[is_table_neg$`Compound Name` == "Cholesterol (d7)"]) <= 90)
    })
  
  if (length(unique(unlist(temp_idx))) >= 1) {
    temp_data <-
      purrr::map2(
        .x = temp_idx,
        .y = peak_table_neg,
        .f = function(x, y) {
          x <-
            x[which.min(abs(y$rt[x] - is_table_neg$rt[is_table_neg$`Compound Name` == "Cholesterol (d7)"]))]
          
          temp_data <-
            y[x, , drop = FALSE] %>%
            dplyr::mutate(peak_name = name) %>%
            dplyr::mutate(name = "Cholesterol") %>%
            dplyr::mutate(Class = "Chol") %>%
            dplyr::mutate(rt = rt / 60) %>%
            dplyr::mutate(mean.int = 10000)
          
          temp_data <-
            temp_data %>%
            dplyr::select(one_of(colnames(data_neg)))
          temp_data
        }
      )
    
    peak_name <-
      lapply(temp_data, function(x) {
        x$peak_name
      }) %>%
      unlist() %>%
      `[`(1)
    
    rt <-
      lapply(temp_data, function(x) {
        x$rt
      }) %>%
      unlist() %>%
      `[`(1)
    
    temp_data <-
      lapply(temp_data, function(x) {
        x$peak_name <- peak_name
        x$rt <- rt
        x
      })
    
    if (length(temp_data) == 1) {
      new_data <- temp_data[[1]]
    } else{
      new_data <- temp_data[[1]]
      for (x in temp_data[-1]) {
        new_data <-
          new_data %>%
          dplyr::left_join(x, by = intersect(colnames(x), colnames(new_data)))
      }
    }
    
    data_neg <-
      data_neg %>%
      dplyr::full_join(new_data, by = intersect(colnames(data_neg), colnames(new_data)))
  }
  
  ###calculate concentration for each samples and lipid--------------------------
  is_tag_neg <-
    is_data_neg[, 1:10]
  
  is_table_neg <-
    is_data_neg[, -c(1:10)]
  
  last_idx_neg <- which(colnames(data_neg) == "mean.int")
  
  lipid_tag_neg <-
    data_neg[, c(1:last_idx_neg)]
  
  lipid_table_neg <-
    data_neg[, -c(1:last_idx_neg)]
  
  lipid_table_neg <-
    lipid_table_neg[, sort(colnames(lipid_table_neg))]
  
  is_table_neg <-
    is_table_neg[, sort(colnames(is_table_neg))]
  
  lipid_tag_neg$rt <-
    lipid_tag_neg$rt * 60
  
  expresion_data_abs_neg <-
    cal_abs(
      lipid_tag = lipid_tag_neg,
      lipid_table = lipid_table_neg,
      is_tag = is_tag_neg,
      is_table = is_table_neg,
      match_item = match_item_neg
    )
  
  express_data_abs_ug_ml_neg <- expresion_data_abs_neg$express_data_abs_ug_ml
  express_data_abs_um_neg <- expresion_data_abs_neg$express_data_abs_um
  variable_info_abs_neg <- expresion_data_abs_neg$variable_info_abs
  
  save(variable_info_abs_neg, file = "lipid_search/NEG/variable_info_abs_neg")
  save(express_data_abs_ug_ml_neg, file = "lipid_search/NEG/express_data_abs_ug_ml_neg")
  save(express_data_abs_um_neg, file = "lipid_search/NEG/express_data_abs_um_neg")
  
  save(lipid_tag_neg, file = "lipid_search/NEG/lipid_tag_neg")
  save(lipid_table_neg, file = "lipid_search/NEG/lipid_table_neg")
  save(is_tag_neg, file = "lipid_search/NEG/is_tag_neg")
  save(is_table_neg, file = "lipid_search/NEG/is_table_neg")
}



combine_pos_neg_quantification <- function(path = ".",
                                           express_data_abs_ug_ml_pos,
                                           express_data_abs_um_pos,
                                           variable_info_abs_pos,
                                           express_data_abs_ug_ml_neg,
                                           express_data_abs_um_neg,
                                           variable_info_abs_neg) {
  dir.create(path, showWarnings = FALSE)
  # ####remove duplicated
  # ###positive mode
  # library(plyr)
  # express_data_abs_ug_ml_pos <-
  #   express_data_abs_ug_ml_pos %>%
  #   data.frame(variable_info_abs_pos, ., check.names = FALSE) %>%
  #   plyr::dlply(.variables = .(name)) %>%
  #   purrr::map(
  #     .f = function(x) {
  #       x <-
  #         x %>%
  #         dplyr::filter(mean.int == max(mean.int)) %>%
  #         head(1)
  #       # dplyr::select(colnames(express_data_abs_ug_ml))
  #     }
  #   ) %>%
  #   do.call(rbind, .) %>%
  #   as.data.frame()
  # 
  # rownames(express_data_abs_ug_ml_pos) <- NULL
  # 
  # express_data_abs_ug_ml_pos <-
  #   express_data_abs_ug_ml_pos %>%
  #   tibble::column_to_rownames(var = "peak_name") %>%
  #   dplyr::select(-colnames(variable_info_abs_pos)[-1])
  # 
  # express_data_abs_um_pos <-
  #   express_data_abs_um_pos %>%
  #   data.frame(variable_info_abs_pos, ., check.names = FALSE) %>%
  #   plyr::dlply(.variables = .(name)) %>%
  #   purrr::map(
  #     .f = function(x) {
  #       x <-
  #         x %>%
  #         dplyr::filter(mean.int == max(mean.int)) %>%
  #         head(1)
  #     }
  #   ) %>%
  #   do.call(rbind, .) %>%
  #   as.data.frame()
  # 
  # rownames(express_data_abs_um_pos)   <- NULL
  # 
  # express_data_abs_um_pos <-
  #   express_data_abs_um_pos %>%
  #   tibble::column_to_rownames(var = "peak_name") %>%
  #   dplyr::select(-colnames(variable_info_abs_pos)[-1])
  # 
  # variable_info_abs_pos <-
  #   variable_info_abs_pos[match(rownames(express_data_abs_ug_ml_pos), variable_info_abs_pos$peak_name),]
  # 
  # ###negative mode
  # library(plyr)
  # express_data_abs_ug_ml_neg <-
  #   express_data_abs_ug_ml_neg %>%
  #   data.frame(variable_info_abs_neg, ., check.names = FALSE) %>%
  #   plyr::dlply(.variables = .(name)) %>%
  #   purrr::map(
  #     .f = function(x) {
  #       x <-
  #         x %>%
  #         dplyr::filter(mean.int == max(mean.int)) %>%
  #         head(1)
  #       # dplyr::select(colnames(express_data_abs_ug_ml))
  #     }
  #   ) %>%
  #   do.call(rbind, .) %>%
  #   as.data.frame()
  # 
  # rownames(express_data_abs_ug_ml_neg) <- NULL
  # 
  # express_data_abs_ug_ml_neg <-
  #   express_data_abs_ug_ml_neg %>%
  #   tibble::column_to_rownames(var = "peak_name") %>%
  #   dplyr::select(-colnames(variable_info_abs_neg)[-1])
  # 
  # express_data_abs_um_neg <-
  #   express_data_abs_um_neg %>%
  #   data.frame(variable_info_abs_neg, ., check.names = FALSE) %>%
  #   plyr::dlply(.variables = .(name)) %>%
  #   purrr::map(
  #     .f = function(x) {
  #       x <-
  #         x %>%
  #         dplyr::filter(mean.int == max(mean.int)) %>%
  #         head(1)
  #     }
  #   ) %>%
  #   do.call(rbind, .) %>%
  #   as.data.frame()
  # 
  # rownames(express_data_abs_um_neg)   <- NULL
  # 
  # express_data_abs_um_neg <-
  #   express_data_abs_um_neg %>%
  #   tibble::column_to_rownames(var = "peak_name") %>%
  #   dplyr::select(-colnames(variable_info_abs_neg)[-1])
  # 
  # variable_info_abs_neg <-
  #   variable_info_abs_neg[match(rownames(express_data_abs_ug_ml_neg), variable_info_abs_neg$peak_name),]
  # 
  # #####for PC, PS and PI, if they are detected in pos and neg, only keep in neg
  # remove_idx <-
  #   which(
  #     variable_info_abs_pos$name %in% variable_info_abs_neg$name &
  #       (
  #         variable_info_abs_pos$Class == "PC" |
  #           variable_info_abs_pos$Class == "PI" |
  #           variable_info_abs_pos$Class == "PS" |
  #           variable_info_abs_pos$Class == "PG" |
  #           variable_info_abs_pos$Class == "LPE" |
  #           variable_info_abs_pos$Class == "PPC" |
  #           variable_info_abs_pos$Class == "PPE"
  #       )
  #   )
  # 
  # if(length(remove_idx) > 0){
  #   variable_info_abs_pos <-
  #     variable_info_abs_pos[-remove_idx,]
  #   express_data_abs_ug_ml_pos <- 
  #     express_data_abs_ug_ml_pos[-remove_idx,]
  #   express_data_abs_um_pos <- 
  #     express_data_abs_um_pos[-remove_idx,]
  # }
  
  ##combine positive and negative
  rownames(express_data_abs_ug_ml_pos) <-
    rownames(express_data_abs_um_pos) <-
    variable_info_abs_pos$peak_name
  
  rownames(express_data_abs_ug_ml_neg) <-
    rownames(express_data_abs_um_neg) <-
    variable_info_abs_neg$peak_name
  
  express_data_abs_ug_ml <-
    rbind(express_data_abs_ug_ml_pos,
          express_data_abs_ug_ml_neg)
  
  colnames(express_data_abs_um_pos) <-
    colnames(express_data_abs_um_neg)
  
  express_data_abs_um <-
    rbind(express_data_abs_um_pos,
          express_data_abs_um_neg)
  
  variable_info_abs <-
    variable_info_abs_pos %>%
    dplyr::full_join(variable_info_abs_neg, by = intersect(
      colnames(variable_info_abs_pos),
      colnames(variable_info_abs_neg)
    ))
  
  ####remove duplicated
  library(plyr)
  express_data_abs_ug_ml <-
    express_data_abs_ug_ml %>%
    data.frame(variable_info_abs, ., check.names = FALSE) %>%
    plyr::dlply(.variables = .(name)) %>%
    purrr::map(
      .f = function(x) {
        x <-
          x %>%
          dplyr::filter(mean.int == max(mean.int)) %>%
          head(1)
        # dplyr::select(colnames(express_data_abs_ug_ml))
      }
    ) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  rownames(express_data_abs_ug_ml)  <- NULL
  
  express_data_abs_ug_ml <-
    express_data_abs_ug_ml %>%
    tibble::column_to_rownames(var = "peak_name") %>%
    dplyr::select(-colnames(variable_info_abs)[-1])
  
  library(plyr)
  express_data_abs_um <-
    express_data_abs_um %>%
    data.frame(variable_info_abs, ., check.names = FALSE) %>%
    plyr::dlply(.variables = .(name)) %>%
    purrr::map(
      .f = function(x) {
        x <-
          x %>%
          dplyr::filter(mean.int == max(mean.int)) %>%
          head(1)
      }
    ) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  rownames(express_data_abs_um)   <- NULL
  
  express_data_abs_um <-
    express_data_abs_um %>%
    tibble::column_to_rownames(var = "peak_name") %>%
    dplyr::select(-colnames(variable_info_abs)[-1])
  
  variable_info_abs <-
    variable_info_abs[match(rownames(express_data_abs_ug_ml),
                            variable_info_abs$peak_name), ]
  
  library(plyr)
  class_data_ug_ml <-
    express_data_abs_ug_ml %>%
    data.frame(
      class = variable_info_abs$Class,
      .,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) %>%
    plyr::dlply(.variables = .(class)) %>%
    purrr::map(
      .f = function(x) {
        apply(x[, -1], 2, function(y) {sum(y, na.rm = TRUE)})
      }
    ) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  write.csv(class_data_ug_ml, file.path(path, "class_data_ug_ml"))
  
  library(plyr)
  class_data_um <-
    express_data_abs_um %>%
    data.frame(
      class = variable_info_abs$Class,
      .,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) %>%
    plyr::dlply(.variables = .(class)) %>%
    purrr::map(
      .f = function(x) {
        apply(x[, -1], 2, function(y) {sum(y, na.rm = TRUE)})
      }
    ) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  write.csv(class_data_um, file.path(path, "class_data_um"))
  
  plot1 <-
    class_data_ug_ml %>%
    tibble::rownames_to_column(var = "class") %>%
    tidyr::pivot_longer(cols = -class,
                        names_to = "sample_id",
                        values_to = "value") %>%
    ggplot(aes(sample_id, value)) +
    geom_bar(stat = "identity", position = "fill", aes(fill = class)) +
    theme_bw() +
    labs(x = "") +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    scale_x_discrete(expand = expansion(mult = c(0, 0))) +
    theme(axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ))
  
  ggsave(
    plot1,
    filename = file.path(path, "plot_ug.pdf"),
    width = 14,
    height = 7
  )
  
  
  plot2 <-
    class_data_um %>%
    tibble::rownames_to_column(var = "class") %>%
    tidyr::pivot_longer(cols = -class,
                        names_to = "sample_id",
                        values_to = "value") %>%
    ggplot(aes(sample_id, value)) +
    geom_bar(stat = "identity",
             position = "fill",
             aes(fill = class)) +
    theme_bw() +
    labs(x = "") +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    scale_x_discrete(expand = expansion(mult = c(0, 0))) +
    theme(axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ))
  
  ggsave(
    plot2,
    filename = file.path(path, "plot_um.pdf"),
    width = 14,
    height = 7
  )
  
  express_data_abs_um_per <-
    express_data_abs_um %>%
    apply(2, function(x) {
      x * 100 / sum(x)
    }) %>%
    as.data.frame()
  
  class_data_um_per <-
    class_data_um %>%
    apply(2, function(x) {
      x * 100 / sum(x)
    }) %>%
    as.data.frame()
  
  write.csv(
    cbind(variable_info_abs, express_data_abs_ug_ml),
    file.path(path, "lipid_data_ug_ml.csv"),
    row.names = FALSE
  )
  
  write.csv(
    cbind(variable_info_abs, express_data_abs_um),
    file.path(path, "lipid_data_um.csv"),
    row.names = FALSE
  )
  
  write.csv(
    cbind(variable_info_abs, express_data_abs_um_per),
    file.path(path, "lipid_data_um_per.csv"),
    row.names = FALSE
  )
  
  write.csv(class_data_ug_ml,
            file.path(path, "lipid_data_class_ug_ml.csv"),
            row.names = TRUE)
  
  write.csv(class_data_um,
            file.path(path, "lipid_data_class_um.csv"),
            row.names = TRUE)
  
  write.csv(class_data_um_per,
            file.path(path, "lipid_data_class_um_per.csv"),
            row.names = TRUE)
  
}










# sxtTools::setwd_project()
# # rm(list = ls())
# setwd("data/lipid20201020/lipid_search/")
#
# library(tidyverse)
#
# file <- "102420_Align_Pos.xlsx"
# path <- "POS/"
#
# data_pos <- tidy_lipidsearch_data(file = file,
#                                    path = path,
#                                    polarity = "positive")
#
# data_neg <- tidy_lipidsearch_data(file = "CW_091520_AlignTestNeg.xlsx",
#                                    path = ".",
#                                    polarity = "negative")

tidy_lipidsearch_data <-
  function(file,
           path = ".",
           polarity = c("positive", "negative")) {
    require(tidyverse)
    format <- stringr::str_split(file, "\\.")[[1]] %>%
      tail(1)
    polarity <- match.arg(polarity)
    ##read lipid search result
    lipid_table <-
      readxl::read_xlsx(path = file.path(path, file),
                        sheet = 1,
                        col_names = FALSE)
    
    idx <- which(is.na(lipid_table$...1))[1]
    rename <- lipid_table[1:(idx - 1), 1]
    
    rename <-
      rename$...1 %>%
      stringr::str_split("\\:") %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::rename(new_name = V1, old_name = V2) %>%
      dplyr::filter(stringr::str_detect(old_name, "raw"))
    
    rename$new_name <-
      rename$new_name %>%
      stringr::str_replace("#", "")
    
    rename$old_name <-
      rename$old_name %>%
      stringr::str_replace("\\.raw", "")
    
    lipid_table <-
      lipid_table[-c(1:idx), ]
    
    colnames(lipid_table) <- as.character(lipid_table[1, ])
    lipid_table <- lipid_table[-1, ]
    
    lipid_table <-
      lipid_table %>%
      dplyr::rename(rentention_time = `GroupTopPos[c]`)
    
    mz <-
      lipid_table %>%
      dplyr::select(contains("ObsMz")) %>%
      apply(1, function(x) {
        mean(as.numeric(x), na.rm = TRUE)
      })
    
    lipid_table <-
      lipid_table %>%
      dplyr::select(-contains("Height")) %>%
      dplyr::select(-contains("Score")) %>%
      dplyr::select(-contains("Norm")) %>%
      dplyr::select(-contains("Top")) %>%
      dplyr::select(-contains("Hwhm")) %>%
      dplyr::select(-contains("DataId")) %>%
      dplyr::select(-contains("Scan")) %>%
      dplyr::select(-contains("ObsMz")) %>%
      dplyr::select(-contains("Rt")) %>%
      dplyr::select(-contains("Delta")) %>%
      dplyr::select(-contains("z")) %>%
      dplyr::select(-contains("It"))
    
    colnames(lipid_table) <-
      colnames(lipid_table) %>%
      stringr::str_replace("Area", "")
    
    expression_data <-
      lipid_table[, rename$new_name] %>%
      apply(2, as.numeric) %>%
      as.data.frame()
    
    colnames(expression_data) <- rename$old_name
    
    variable_info <-
      lipid_table %>%
      dplyr::select(-rename$new_name)
    
    adduct <-
      variable_info$LipidIon %>%
      purrr::map(
        .f = function(x) {
          temp <- stringr::str_split(x, "\\)\\+", n = 2)[[1]]
          if (length(temp) == 1) {
            temp <- stringr::str_split(x, "\\)\\-", n = 2)[[1]]
          }
          temp
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::rename(name = V1, adduct = V2) %>%
      dplyr::mutate(name = paste(name, ")", sep = ""))
    
    adduct <-
      data.frame(
        lipid_raw_name = variable_info$LipidIon,
        adduct,
        polarity,
        stringsAsFactors = FALSE
      )
    
    variable_info <-
      data.frame(adduct, variable_info, stringsAsFactors = FALSE)
    
    mean.int <- apply(expression_data, 1, median)
    
    variable_info <-
      variable_info %>%
      dplyr::select(-c(ARatio.0.:HDiff.0.)) %>%
      dplyr::select(
        name,
        lipid_raw_name,
        rt = rentention_time,
        adduct,
        polarity,
        Class,
        FattyAcid,
        IonFormula,
        contains("FA")
        # mean.int = Average.
      ) %>%
      data.frame(., mean.int, stringsAsFactors = FALSE)
    
    
    variable_info$rt <- as.numeric(variable_info$rt)
    data <- cbind(variable_info, expression_data)
    
    ##give unique name for each lipid
    if (polarity == "positive") {
      data$peak_name <-
        paste("peak_pos", 1:nrow(data), sep = "")
    } else{
      data$peak_name <-
        paste("peak_neg", 1:nrow(data), sep = "")
    }
    
    rownames(data) <- data$peak_name
    
    data <-
      data %>%
      dplyr::mutate(mz = mz) %>%
      dplyr::select(peak_name, mz, rt, dplyr::everything())
    
    return(data)
  }

# data <- combine_pos_neg(data_pos, data_neg)

combine_pos_neg <- function(data_pos, data_neg) {
  ##combine same lipid
  data_pos <-
    data_pos %>%
    dplyr::arrange(name)
  
  # library(plyr)
  #
  # data_pos <-
  #   data_pos %>%
  #   plyr::dlply(.variables = .(name))
  
  # data_pos <-
  #   lapply(data_pos, function(x) {
  #     if (nrow(x) == 1) {
  #       return(x)
  #     }
  #     if (any(as.numeric(x$rt) > 1)) {
  #       x <-
  #         x %>%
  #         dplyr::filter(as.numeric(rt) > 1)
  #
  #       if (nrow(x) == 1) {
  #         return(x)
  #       }
  #     }
  #     x <-
  #       x %>%
  #       dplyr::filter(as.numeric(x$mean.int) == max(as.numeric(x$mean.int)))
  #   })
  
  data_neg <-
    data_neg %>%
    dplyr::arrange(name)
  
  # library(plyr)
  #
  # data_neg <-
  #   data_neg %>%
  #   plyr::dlply(.variables = .(name))
  #
  # data_neg <-
  #   lapply(data_neg, function(x) {
  #     if (nrow(x) == 1) {
  #       return(x)
  #     }
  #     if (any(as.numeric(x$rt) > 1)) {
  #       x <-
  #         x %>%
  #         dplyr::filter(as.numeric(rt) > 1)
  #
  #       if (nrow(x) == 1) {
  #         return(x)
  #       }
  #     }
  #     x <-
  #       x %>%
  #       dplyr::filter(as.numeric(x$mean.int) == max(as.numeric(x$mean.int)))
  #   })
  #
  #
  # data_pos <-
  #   data_pos %>%
  #   do.call(rbind, .)
  #
  # data_neg <-
  #   data_neg %>%
  #   do.call(rbind, .)
  
  data <-
    data_pos %>%
    dplyr::full_join(data_neg, by = intersect(colnames(data_pos), colnames(data_neg)))
  
  return(data)
}



###filter variables
filter_variable <- function(data,
                            rt.cutoff = 1) {
  data <-
    data %>%
    dplyr::filter(rt > rt.cutoff)
}




extract_is <- function(polarity = c("positive", "negative"),
                       path = ".",
                       is_table_name = "IS_information_092920_new.xlsx") {
  is_table <- readxl::read_xlsx(is_table_name)
  mass <- as.numeric(is_table$`Exact Mass`)
  library(Rdisop)
  
  mz_pos_h <- mass + Rdisop::getMass(molecule = getMolecule("H"))
  mz_pos_na <- mass + Rdisop::getMass(molecule = getMolecule("Na"))
  mz_pos_nh4 <-
    mass + Rdisop::getMass(molecule = getMolecule("NH4"))
  mz_pos_h_h20 <-
    mass - Rdisop::getMass(molecule = getMolecule("HO"))
  mz_pos_nh4_h20 <-
    mass + Rdisop::getMass(molecule = getMolecule("NH4")) -
    Rdisop::getMass(molecule = getMolecule("H2O"))
  
  mz_neg_h <- mass - Rdisop::getMass(molecule = getMolecule("H"))
  mz_neg_ch3cOO <-
    mass + Rdisop::getMass(molecule = getMolecule("CH3COO"))
  mz_neg_hcooh <-
    mass + Rdisop::getMass(molecule = getMolecule("HCOOH"))
  
  is_table <- data.frame(
    is_table,
    mz_pos_h,
    mz_pos_na,
    mz_pos_nh4,
    mz_pos_h_h20,
    mz_pos_nh4_h20,
    mz_neg_h,
    mz_neg_ch3cOO,
    stringsAsFactors = FALSE
  )
  
  write.csv(is_table, "is_table.csv", row.names = FALSE)
  
  if (polarity == "positive") {
    library(metflow2)
    new_path <- "mzXML/POS"
    
    is_table_pos <-
      rbind(
        data.frame(
          name = is_table[, 2],
          mz = mz_pos_h,
          adduct = "H"
        ),
        cbind(
          name = is_table[, 2],
          mz = mz_pos_na,
          adduct = 'Na'
        ),
        cbind(
          name = is_table[, 2],
          mz = mz_pos_nh4,
          adduct = "NH4"
        ),
        cbind(
          name = is_table[, 2],
          mz = mz_pos_h_h20,
          adduct = "+H-H2O"
        ),
        cbind(
          name = is_table[, 2],
          mz = mz_pos_nh4_h20,
          adduct = "+NH4-H2O"
        )
      ) %>%
      as.data.frame() %>%
      dplyr::arrange(name, adduct)
    
    is_table_pos$name <-
      stringr::str_trim(is_table_pos$name)
    
    is_table_pos$mz <- as.numeric(is_table_pos$mz)
    
    xlsx::write.xlsx(is_table_pos,
                     file.path(new_path, "is_table_pos.xlsx"),
                     row.names = FALSE)
    
    result_pos <-
      metflow2::extractPeaks(
        path = new_path,
        ppm = 25,
        threads = 3,
        is.table = "is_table_pos.xlsx"
      )
    
    plot_pos_path <- paste(polarity, "plot", sep = "_")
    dir.create(plot_pos_path, showWarnings = FALSE)
    for (i in 1:nrow(is_table_pos)) {
      cat(i, " ")
      x <- as.character(is_table_pos[i,])
      plot_path <- file.path(plot_pos_path, x[1])
      dir.create(plot_path, showWarnings = FALSE)
      name <- paste(x[1], x[3], sep = "_")
      
      plot <-
        metflow2::showPeak(object = result_pos,
                           peak.index = i,
                           title = name)
      
      library(htmlwidgets)
      
      saveWidget(
        widget = plot,
        file = file.path(getwd(), plot_path, paste(name, "html", sep = ".")),
        selfcontained = TRUE,
        libdir = "lib"
      )
    }
  }
  
  
  if (polarity == "negative") {
    new_path <- "mzXML/neg"
    
    is_table_neg <-
      rbind(
        data.frame(
          name = is_table[, 2],
          mz = mz_neg_h,
          adduct = "H"
        ),
        cbind(
          name = is_table[, 2],
          mz = mz_neg_ch3cOO,
          adduct = 'CH3COO'
        ),
        cbind(
          name = is_table[, 2],
          mz = mz_neg_hcooh,
          adduct = "HCOOH"
        )
      ) %>%
      as.data.frame() %>%
      dplyr::arrange(name, adduct)
    
    is_table_neg$name <-
      stringr::str_trim(is_table_neg$name)
    
    is_table_neg$mz <- as.numeric(is_table_neg$mz)
    
    xlsx::write.xlsx(is_table_neg,
                     file.path(new_path, "is_table_neg.xlsx"),
                     row.names = FALSE)
    
    result_neg <-
      metflow2::extractPeaks(
        path = new_path,
        ppm = 25,
        threads = 3,
        is.table = "is_table_neg.xlsx"
      )
    
    plot_neg_path <- paste(polarity, "plot", sep = "_")
    dir.create(plot_neg_path, showWarnings = FALSE)
    for (i in 1:nrow(is_table_neg)) {
      cat(i, " ")
      x <- as.character(is_table_neg[i,])
      plot_path <- file.path(plot_neg_path, x[1])
      dir.create(plot_path, showWarnings = FALSE)
      name <- paste(x[1], x[3], sep = "_")
      
      plot <-
        metflow2::showPeak(object = result_neg,
                           peak.index = i,
                           title = name)
      
      library(htmlwidgets)
      
      saveWidget(
        widget = plot,
        file = file.path(getwd(), plot_path, paste(name, "html", sep = ".")),
        selfcontained = TRUE,
        libdir = "lib"
      )
    }
  }
}















extract_peak <- function(path = ".",
                         peak_table,
                         match_item,
                         polarity = c("positive", "negative")) {
  polarity <- match.arg(polarity)
  
  library(metflow2)
  if (polarity == "positive") {
    new_path <- "mzXML/POS"
  } else{
    new_path <- "mzXML/NEG"
  }
  
  result <-
    metflow2::extractPeaks(
      path = new_path,
      ppm = 25,
      threads = 3,
      mz = peak_table$mz,
      rt = peak_table$rt,
      rt.tolerance = 180
    )
  
  plot_path <- "intensity_plot"
  dir.create(plot_path, showWarnings = FALSE)
  
  for (i in names(match_item)) {
    is_name <- match_item[[i]]
    if (is.na(is_name)) {
      next()
    }
    
    new_path <- file.path(plot_path, i)
    
    dir.create(new_path, showWarnings = FALSE)
    
    temp_idx <- which(peak_table$Class == i)
    
    if (length(temp_idx) > 0) {
      for (j in temp_idx) {
        x <- as.character(peak_table[j, ])
        name <- x[1] %>%
          stringr::str_replace_all("\\/", "_") %>%
          paste("name-", ., "_mz-", x[2], "_rt-", x[3], sep = "")
        
        plot <-
          metflow2::showPeak(object = result,
                             peak.index = j,
                             title = name)
        
        library(htmlwidgets)
        try(expr = saveWidget(
          widget = plot,
          file = file.path(
            getwd(),
            new_path,
            paste(
              x[4] %>% stringr::str_replace_all("\\/", "_"),
              "html",
              sep = "."
            )
          ),
          selfcontained = TRUE,
          libdir = "lib"
        ),
        silent = TRUE)
      }
    }
  }
  
}






# trans_is_table(is_table_name = "IS_information_20201020.xlsx", polarity = "positive")

trans_is_table <-
  function(is_table_name = "IS_information_20201020.xlsx",
           polarity = c("positive", "negative")) {
    polarity <- match.arg(polarity)
    is_table <- readxl::read_xlsx(is_table_name)
    
    is_table$`Compound Name` <-
      is_table$`Compound Name` %>%
      stringr::str_trim()
    
    is_table_old <- is_table
    mass <- as.numeric(is_table$`Exact Mass`)
    library(Rdisop)
    
    if (polarity == "positive") {
      mz_pos_h <- mass + Rdisop::getMass(molecule = getMolecule("H"))
      mz_pos_na <-
        mass + Rdisop::getMass(molecule = getMolecule("Na"))
      mz_pos_nh4 <-
        mass + Rdisop::getMass(molecule = getMolecule("NH4"))
      mz_pos_h_h20 <-
        mass - Rdisop::getMass(molecule = getMolecule("HO"))
      mz_pos_nh4_h20 <-
        mass + Rdisop::getMass(molecule = getMolecule("NH4")) -
        Rdisop::getMass(molecule = getMolecule("H2O"))
      
      is_table <-
        rbind(
          cbind(
            name = is_table$`Compound Name`,
            mz = mz_pos_h,
            rt = is_table$`RT POS (second)`,
            adduct = "M+H"
          ),
          cbind(
            name = is_table$`Compound Name`,
            mz = mz_pos_na,
            rt = is_table$`RT POS (second)`,
            adduct = 'M+Na'
          ),
          cbind(
            name = is_table$`Compound Name`,
            mz = mz_pos_nh4,
            rt = is_table$`RT POS (second)`,
            adduct = "M+NH4"
          ),
          cbind(
            name = is_table$`Compound Name`,
            mz = mz_pos_h_h20,
            rt = is_table$`RT POS (second)`,
            adduct = "M+H-H2O"
          ),
          cbind(
            name = is_table$`Compound Name`,
            mz = mz_pos_nh4_h20,
            rt = is_table$`RT POS (second)`,
            adduct = "M+NH4-H2O"
          )
        ) %>%
        as.data.frame() %>%
        dplyr::arrange(name, adduct)
      
      is_table$name <-
        stringr::str_trim(is_table$name)
      is_table$mz <- as.numeric(is_table$mz)
      is_table$rt <- as.numeric(is_table$rt)
      
      is_table <-
        is_table %>%
        dplyr::filter(!is.na(rt))
      
      library(plyr)
      is_table <-
        is_table %>%
        plyr::dlply(.variables = .(name)) %>%
        purrr::map(
          .f = function(x) {
            temp <-
              is_table_old %>%
              dplyr::filter(`Compound Name` == x$name[1])
            x <-
              x %>%
              dplyr::filter(adduct == temp$`Adduct POS`)
            x
          }
        ) %>%
        do.call(rbind, .) %>%
        as.data.frame()
      
      
      xlsx::write.xlsx(is_table,
                       file.path("is_table.xlsx"),
                       row.names = FALSE)
    } else{
      mz_neg_h <- mass - Rdisop::getMass(molecule = getMolecule("H"))
      mz_neg_ch3cOO <-
        mass + Rdisop::getMass(molecule = getMolecule("CH3COO"))
      mz_neg_hcooh <-
        mass + Rdisop::getMass(molecule = getMolecule("HCOOH"))
      
      is_table <-
        rbind(
          cbind(
            name = is_table$`Compound Name`,
            mz = mz_neg_h,
            rt = is_table$`RT NEG (second)`,
            adduct = "M-H"
          ),
          cbind(
            name = is_table$`Compound Name`,
            mz = mz_neg_ch3cOO,
            rt = is_table$`RT NEG (second)`,
            adduct = 'M+CH3COO'
          ),
          cbind(
            name = is_table$`Compound Name`,
            mz = mz_neg_hcooh,
            rt = is_table$`RT NEG (second)`,
            adduct = "M+HCOO"
          )
        ) %>%
        as.data.frame() %>%
        dplyr::arrange(name, adduct)
      
      is_table$name <-
        stringr::str_trim(is_table$name)
      is_table$mz <- as.numeric(is_table$mz)
      is_table$rt <- as.numeric(is_table$rt)
      
      is_table <-
        is_table %>%
        dplyr::filter(!is.na(rt)) %>%
        dplyr::filter(rt != "")
      
      library(plyr)
      is_table <-
        is_table %>%
        plyr::dlply(.variables = .(name)) %>%
        purrr::map(
          .f = function(x) {
            temp <-
              is_table_old %>%
              dplyr::filter(`Compound Name` == x$name[1])
            x <-
              x %>%
              dplyr::filter(adduct == temp$`Adduct NEG`)
            x
          }
        ) %>%
        do.call(rbind, .) %>%
        as.data.frame()
      
      xlsx::write.xlsx(is_table,
                       file.path("is_table.xlsx"),
                       row.names = FALSE)
    }
  }


# extract_targeted_peaks(fit.gaussian = TRUE, integrate_xcms = TRUE)

extract_targeted_peaks <-
  function(path = ".",
           targeted_targeted_peak_table_name = "is_table.xlsx",
           forced_targeted_peak_table_name = "forced_table.xlsx",
           from_lipid_search = FALSE,
           sample_numer_show = 5,
           fit.gaussian = TRUE,
           integrate_xcms = FALSE,
           output_eic = TRUE,
           output_integrate = TRUE,
           ppm = 25,
           rt.tolerance = 90) {
    
    peak_table <-
      readxl::read_xlsx(file.path(path , targeted_targeted_peak_table_name))
    
    library(tidyverse)
    
    if(any(dir(path) == "result_raw")){
      load(file.path(path, "result_raw"))
    }else{
      result_raw <-
        extractPeaks2(
          path = path,
          ppm = ppm,
          threads = 3,
          is.table = targeted_targeted_peak_table_name,
          rt.tolerance = rt.tolerance,
          msLevel = 1L
          # filled = TRUE,
          # include = "any"
        )
      
      save(result_raw, file = file.path(path, "result_raw"))
      
    }
    
    dir.create(file.path(path, "peak_shape"))
    
    # load("raw_data")
    # 
    # raw_data %>%
    #   filterRt(rt = as.numeric(result_raw@featureData@data[1,3:4])) %>%
    #   filterMz(mz = as.numeric(result_raw@featureData@data[1,1:2])) %>%
    #   plot(type = "XIC")
    # x <- 
    #   raw_data@featureData@data
    # 
    # x <- 
    # x %>% 
    #   dplyr::filter(fileIdx == 1)
    
    library(xcms)
    
    if(any(dir(path) == "result")){
      load(file.path(path, "result"))
    }else{
      result <- findChromPeaks(
        object = result_raw,
        param = CentWaveParam(
          peakwidth = c(5, 30),
          snthresh = 2,
          ppm = ppm,
          fitgauss = FALSE,
          noise = 500,
          prefilter = c(3, 500),
          integrate = 1
          # mzdiff = -0.01
        )
      )
      
      # plot(result, col = "black", type = "b")
      
      mpp <- MergeNeighboringPeaksParam(expandRt = 3)
      
      result <- refineChromPeaks(result, param = mpp)
      # plot(result, col = "black", type = "b")  
      save(result, file = file.path(path, file.path(path, "result")))  
    }
    
    peak_value <-
      chromPeaks(result) %>%
      as.data.frame()
    
    library(plyr)
    
    peak_value <- 
      peak_value %>% 
      plyr::dlply(.variables = .(row, column)) %>% 
      purrr::map(function(x){
        x <-
          x %>% 
          dplyr::filter(into == max(into))
        x
      }) %>% 
      do.call(rbind, .) %>% 
      as.data.frame()
    
    peak_value$name <-
      peak_table$name[peak_value$row]
    
    peak_value$sample <-
      result@phenoData@data$sample_group[peak_value$column]
    
    peak_value <-
      peak_value %>%
      dplyr::select(-c(row, column)) %>%
      dplyr::select(name, sample, everything())
    
    raw_info <- result@.Data
    
    rownames(raw_info) <- peak_table$name
    colnames(raw_info) <- result@phenoData@data$sample_group
    # raw_info_old <- raw_info
    
    ##output quantification information
    output_quantification_table <- 
      peak_value %>% 
      dplyr::select(name, sample, into, rt) %>% 
      dplyr::group_by(name) %>%
      dplyr::mutate(rt2 = median(rt)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-rt) %>%
      dplyr::rename(rt= rt2) %>%
      dplyr::distinct(name, sample, .keep_all = TRUE) %>% 
      tidyr::pivot_wider(names_from = sample, values_from = "into")
    
    diff_name <- 
      setdiff(
        result@phenoData@data$sample_group,
        colnames(output_quantification_table)
      )
    
    if(length(diff_name) > 0) {
      diff_matrix <-
        matrix(NA,
               nrow = nrow(output_quantification_table),
               ncol = length(diff_name)) %>%
        as.data.frame()
      colnames(diff_matrix) <- 
        diff_name
      output_quantification_table <- 
        cbind(output_quantification_table, diff_matrix)
    }
    
    output_quantification_table <-
      as.data.frame(output_quantification_table)
    
    output_quantification_table_old <- output_quantification_table
    
    if(!integrate_xcms){
      output_quantification_table[,-c(1:2)] <- NA
    }
    
    ###check missing value
    if(fit.gaussian){
      for(temp_name in output_quantification_table$name){
        cat(temp_name, "\n")
        for(temp_sample in colnames(output_quantification_table)[-c(1:2)]){
          # cat(temp_sample, " ")
          value <- 
            output_quantification_table[output_quantification_table$name == temp_name, temp_sample] %>% 
            as.numeric()
          
          if(!is.na(value)){
            next()
          }else{
            temp_raw_info <-
              raw_info[temp_name, temp_sample][[1]] 
            
            xy <- 
              data.frame(x = temp_raw_info@rtime, 
                         y = temp_raw_info@intensity, 
                         stringsAsFactors = FALSE) %>% 
              dplyr::filter(!is.na(y))
            
            if(sum(!is.na(xy$y)) <= 6){
              value <- 0
            }else{
              fit_result <-
                try(
                  expr = fit_gaussian(x = xy$x, y = xy$y), 
                  silent = TRUE
                )
              
              if(class(fit_result$y) == "try-error"){
                xy <-
                  data.frame(xy, z = 0)
              } else{
                xy <-
                  data.frame(xy, z = fit_result$y)                
              }
              
              
              # plot(xy$x, xy$y, type = "b")
              # points(xy$x, xy$z, type = "b", col = "red")
              # 
              # abline(v = fit_result$center)
              # abline(v = fit_result$center - fit_result$width * 2)
              # abline(v = fit_result$center + fit_result$width * 2)
              # abline(h = 0)
              
              l1 <- loess(y ~ x, xy, control = loess.control(surface = "direct"))
              # l2 <-
              #   loess(z ~ x, xy, control = loess.control(surface = "direct"))
              f1 <- function(x){
                predict(l1, newdata = x)
              }
              # f1 <- approxfun(xy$x, xy$y)
              # f2 <- function(x)
              #   predict(l2, newdata = x)
              # f2 <- approxfun(xy$x, xy$z)
              
              value1 <- try(expr = integrate(f = f1,
                                             lower = min(xy$x),
                                             upper = max(xy$x))$value,
                            silent = TRUE
              )
              # 
              # value2 <- try(integrate(f = f2,
              #                         lower = min(xy$x),
              #                         upper = max(xy$x))$value, silent = TRUE
              # )
              
              if(class(value1) == "try-error"){
                value1 <- NA
              }
              # 
              # if(class(value2) == "try-error"){
              #   value2 <- NA
              # }
              
              fit_error <- 
                abs(as.numeric(fit_result$residual))/xy$y
              
              max_idx <- which.max(xy$z)
              
              item1 <- 
                try(expr = shapiro.test(x = xy$z)$p.value > 0.05, silent = TRUE)
              
              if(class(item1) == "try-error"){
              item1 <- FALSE  
              }
                
              item2 <- 
                (sum(fit_error < 0.3) >= nrow(xy)/2) &
                (abs(max_idx - nrow(xy)/2) < nrow(xy)*0.3/2)
              
              if(item1 | item2){
                value <- value1
                ##add information to raw_info
                raw_info[temp_name, temp_sample][[1]]@chromPeaks <- 
                  data.frame(
                    rt = fit_result$center,
                    rtmin = min(xy$x),
                    rtmax = max(xy$x),
                    into = value,
                    maxo = max(xy$y),
                    sn = NA
                  ) %>% 
                  as.matrix()
                
                raw_info[temp_name, temp_sample][[1]]@rtime <-
                  xy$x
                raw_info[temp_name, temp_sample][[1]]@intensity <-
                  xy$y
                
                raw_info[temp_name, temp_sample][[1]]@productMz <-
                  xy$z
              }else{
                value <- value  
              }
            }
            output_quantification_table[output_quantification_table$name == temp_name, temp_sample] <-
              value
          }
        }
      }  
    }
    
    ###add mz and adduct
    output_quantification_table <-
      output_quantification_table %>% 
      dplyr::left_join(peak_table[,c("name", "mz", "adduct")], by = "name") %>% 
      dplyr::select(name, mz, rt, adduct, everything())
    
    if(from_lipid_search){
      output_quantification_table <- 
        output_quantification_table %>% 
        dplyr::left_join(peak_table[,c("name", "compound_name")],
                         by = "name")
      
      ##remove duplicated
      output_quantification_table <-
        output_quantification_table %>%
        plyr::dlply(.variables = .(compound_name)) %>%
        purrr::map(function(y) {
          if (nrow(y) == 1) {
            return(y)
          }
          
          mean.int <- 
            y %>% 
            dplyr::select(-c(name:adduct, compound_name)) %>% 
            apply(1, function(z){
              mean(z, na.rm = TRUE)
            })
          y <- y %>% 
            dplyr::mutate(mean.int = mean.int)
          
          ###remain the peaks with largest mean.int          
          if(nrow(y) > 1) {
            y <- 
              y %>% 
              dplyr::filter(mean.int == max(mean.int)) %>% 
              `[`(.,1,) %>% 
              as.data.frame()
            
            return(y %>% dplyr::select(-mean.int)) 
          }
        }) %>% 
        do.call(rbind, .) %>% 
        as.data.frame()
      
      output_quantification_table <- 
        output_quantification_table %>% 
        dplyr::select(-compound_name)
      
    } 
    
    ###manual check
    # browser()
    if(!is.null(forced_targeted_peak_table_name)){
      cat("Manual check..\n")
      forced_targeted_peak_table <- 
        readxl::read_xlsx(file.path(path, forced_targeted_peak_table_name)) %>% 
        dplyr::filter(!is.na(begin_rt)) %>% 
        dplyr::filter(begin_rt != "")
      
      if(nrow(forced_targeted_peak_table) > 0){
        for(i in 1:nrow(forced_targeted_peak_table)){
          cat(forced_targeted_peak_table$name[i], "\n")  
          temp_name <- forced_targeted_peak_table$name[i]
          temp_sample <- forced_targeted_peak_table$sample[i]
          
          temp_raw_info <-
            raw_info[temp_name, temp_sample][[1]] 
          
          xy <- 
            data.frame(x = temp_raw_info@rtime, 
                       y = temp_raw_info@intensity, 
                       stringsAsFactors = FALSE) %>% 
            dplyr::filter(!is.na(y)) %>% 
            dplyr::filter(x >= as.numeric(forced_targeted_peak_table$begin_rt[i]) &
                            x <= as.numeric(forced_targeted_peak_table$end_rt[i]))
          
          # fit_result <-
          #   try(expr = fit_gaussian(x = xy$x, y = xy$y), silent = TRUE
          #   )
          # 
          # xy <-
          #   data.frame(xy, z = fit_result$y)
          
          l1 <- loess(y ~ x, xy, 
                      control = loess.control(surface = "direct"))
          
          f1 <- function(x){
            predict(l1, newdata = x)
          }
          
          value1 <- try(expr = integrate(f = f1,
                                         lower = min(xy$x),
                                         upper = max(xy$x))$value,
                        silent = TRUE
          )
          
          if(class(value1) == "try-error"){
            value1 <- NA
          }
          
          value <- value1
          ##add information to raw_info
          raw_info[temp_name, temp_sample][[1]]@chromPeaks <- 
            data.frame(
              rt = fit_result$center,
              rtmin = min(xy$x),
              rtmax = max(xy$x),
              into = value,
              maxo = max(xy$y),
              sn = NA
            ) %>% 
            as.matrix()
          
          raw_info[temp_name, temp_sample][[1]]@rtime <-
            xy$x
          
          raw_info[temp_name, temp_sample][[1]]@intensity <-
            xy$y
          
          # raw_info[temp_name, temp_sample][[1]]@productMz <-
          #   xy$z
          
          output_quantification_table[output_quantification_table$name == temp_name, temp_sample] <-
            value
          
        }  
      }
      
    }
    
    openxlsx::write.xlsx(x = output_quantification_table,
                         file = file.path(path, "quantification_table.xlsx"),
                         asTable = TRUE)
    
    example_temp <-
      matrix(nrow = 0, ncol = 8) %>% 
      as.data.frame()
    
    colnames(example_temp)  <-
      c("name",
        "rt",
        "rtmin",
        "rtmax",
        "into",
        "intb",
        "maxo",
        "sn")
    
    peak_value_old <- peak_value
    
    peak_value <-
      purrr::map2(
        .x = as.data.frame(raw_info),
        .y = colnames(as.data.frame(raw_info)),
        .f = function(data, sample_name) {
          names(data) <- rownames(raw_info)
          temp <-
            purrr::map2(
              .x = data,
              .y = names(data),
              .f = function(data1, peak_name) {
                if (nrow(as.data.frame(data1@chromPeaks)) == 0) {
                  peak_name <- logical(0)
                }
                test <-
                  data.frame(name = peak_name,
                             as.data.frame(data1@chromPeaks) %>% 
                               dplyr::filter(into == max(into)),
                             stringsAsFactors = FALSE)
                
                if(ncol(test) < ncol(example_temp)){
                  diff_name <- setdiff(colnames(example_temp), colnames(test))
                  add_matrix <- 
                    matrix(nrow = nrow(test), ncol = length(diff_name)) %>% 
                    as.data.frame()
                  colnames(add_matrix) <- diff_name
                  test <- data.frame(test, add_matrix)[,colnames(example_temp)]
                  test
                }
                test
              }
            )  %>%
            do.call(rbind, .) %>%
            as.data.frame()
          
          if (nrow(temp) == 0) {
            sample_name <- logical(0)
          }
          
          data.frame(sample = sample_name, temp, stringsAsFactors = FALSE)
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    rownames(peak_value) <- NULL
    
    peak_table <-
      peak_table %>%
      dplyr::filter(name %in% output_quantification_table$name)
    
    # browser()
    
    forced_targeted_peak_table_temple <- 
      matrix(NA, nrow = nrow(raw_info), ncol = ncol(raw_info)) %>% 
      as.data.frame()
    
    colnames(forced_targeted_peak_table_temple) <- colnames(raw_info)
    rownames(forced_targeted_peak_table_temple) <- rownames(raw_info)
    
    forced_targeted_peak_table_temple <- 
    forced_targeted_peak_table_temple %>% 
      tibble::rownames_to_column(var = "name") %>% 
      tidyr::pivot_longer(cols = -name, names_to = "sample", values_to = "value") %>% 
      dplyr::mutate(begin_rt = NA, end_rt = NA) %>% 
      dplyr::select(-value)
    
    forced_targeted_peak_table_temple <- 
    forced_targeted_peak_table_temple %>% 
      dplyr::filter(name %in% output_quantification_table$name)
    
    openxlsx::write.xlsx(forced_targeted_peak_table_temple, 
                         file.path(path, "forced_targeted_peak_table_temple.xlsx"),
                         asTable = TRUE)
    
    if(output_eic){
      cat("Output peak shapes...\n")
      setwd("peak_shape")
      for (i in 1:nrow(peak_table)) {
       test <- 
        try(expr = {
          if(!is.null(forced_targeted_peak_table_name)){
            if(!peak_table$name[i] %in% forced_targeted_peak_table$name){
              next()
            }
          }
          
          cat(i, " ")
          x <- as.character(peak_table[i,])
          name <- paste(x[1], sep = "_")
          
          plot <-
            showPeak2(object = result_raw,
                      peak.index = match(x[1],rownames(raw_info)),
                      title = name, 
                      interactive = TRUE,
                      area_data = peak_value[peak_value$name == name,],
                      raw_info = raw_info)
          
          library(htmlwidgets)
          
          saveWidget(
            widget = plot,
            file = file.path(path, paste(name, "html", sep = ".")),
            selfcontained = TRUE,
            title = name,
            libdir = "lib"
          )
       })
       
       if(class(test) == "try-error"){
         setwd("..")
         stop("error")
       }
        
      }
      cat("Done\n")
      setwd("..")
    }
  
    
  }

# x <- temp_raw_info@rtime
# y <- temp_raw_info@intensity
# 
# xy <- data.frame(x, y, stringsAsFactors = FALSE) %>% 
#   dplyr::filter(!is.na(y))
# 
# x <- xy$x
# y <- xy$y
# 
# plot(x,y, type = "b")
# 
fit_gaussian <-
  function (x,
            y,
            start.center = NULL,
            start.width = NULL,
            start.height = NULL,
            start.floor = NULL,
            fit.floor = FALSE) {
    who.max <- which.max(y)
    if (is.null(start.center)){
      start.center <- x[who.max]
    }
      
    if (is.null(start.height)){
      start.height <- y[who.max]
    }
      
    if (is.null(start.width)){
      start.width <- sum(y > (start.height / 2), na.rm = TRUE) / 2
    }
      
    controlList <- nls.control(maxiter = 100,
                               minFactor = 1 / 512,
                               warnOnly = TRUE)
    if (!fit.floor) {
      starts <- list(center = start.center,
                     width = start.width,
                     height = start.height)
      nlsAns <- try(nls(y ~ gaussian(x, center, width, height),
                        start = starts,
                        control = controlList))
    } else {
      if (is.null(start.floor)){
        start.floor <- quantile(y, seq(0, 1, 0.1))[2]
      }
      starts <- list(
        center = start.center,
        width = start.width,
        height = start.height,
        floor = start.floor
      )
      nlsAns <- try(nls(
        y ~ gaussian(x, center, width, height,
                     floor),
        start = starts,
        control = controlList
      ))
    }
    if (class(nlsAns) == "try-error") {
      centerAns <- start.center
      widthAns <- start.width
      heightAns <- start.height
      floorAns <- if (fit.floor){
        start.floor
      }else{
        0
      }
      yAns <-
        try(expr = gaussian(x, centerAns, widthAns, heightAns, floorAns),
            silent = TRUE
        )
      if(class(yAns) == "try-error"){
        residualAns <- y - 0
      }else{
        residualAns <- y - yAns  
      }
      
    } else {
      coefs <- coef(nlsAns)
      centerAns <- coefs[1]
      widthAns <- coefs[2]
      heightAns <- coefs[3]
      floorAns <- if (fit.floor){
        coefs[4]
      } else{
        0
      }

      yAns <- fitted(nlsAns)
      residualAns <- residuals(nlsAns)
    }
    widthAns <- abs(widthAns)
    out <-
      list(
        center = centerAns,
        width = widthAns,
        height = heightAns,
        y = yAns,
        residual = residualAns
      )
    if (fit.floor) {
      out <- c(out, floor = floorAns)
    }
    return(out)
  }


find_local_minimal <- function(x, y){
  i.mins  <- which(diff(sign(diff(c(
    Inf, y, Inf
  )))) == 2)
  i.mins
}









setGeneric(
  name = "extractPeaks2",
  def = function(path = ".",
                 ppm = 15,
                 threads = 4,
                 is.table = "is.xlsx",
                 mz = NULL,
                 rt = NULL,
                 rt.tolerance = 40,
                 msLevel = 1,
                 filled = TRUE,
                 include = c("apex_within", "any", "none")
                 ) {
    options(warn = -1)
    output_path <- path
    include <- match.arg(include)
    # dir.create(output_path)
    ##peak detection
    
    f.in <- list.files(
      path = path,
      pattern = '\\.(mz[X]{0,1}ML|cdf)',
      recursive = TRUE,
      full.names = TRUE
    )
    
    sample_group <-
      BiocGenerics::basename(f.in) %>%
      stringr::str_replace("\\.(mz[X]{0,1}ML|cdf)", "")
    
    pd <-
      data.frame(# sample_name = sub(
        basename(f.in),
        #   pattern = ".mzXML",
        #   replacement = "",
        #   fixed = TRUE
        # ),
        sample_group = sample_group,
        stringsAsFactors = FALSE)
    
    # requireNamespace("xcms")
    cat(crayon::green("Reading raw data, it will take a while...\n"))
    
    if (any(dir(path) == "raw_data")) {
      cat(crayon::yellow("Use old data.\n"))
      load(file.path(path, "raw_data"))
    } else{
      raw_data <- MSnbase::readMSData(
        files = f.in,
        pdata = new("NAnnotatedDataFrame", pd),
        mode = "onDisk",
        verbose = TRUE
      )
      
      save(raw_data,
           file = file.path(output_path, "raw_data"),
           compress = "xz")
    }
    
    cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
    
    is.table <-
      try(readxl::read_xlsx(file.path(path, is.table)), silent = TRUE)
    
    if (!is.null(mz) & !is.null(rt)) {
      if (length(mz) != length(rt)) {
        cat(crayon::yellow("Lenght of mz and rt you provied are different.\n"))
      }
      is.table <- data.frame(mz = as.numeric(mz),
                             rt = as.numeric(rt),
                             stringsAsFactors = FALSE)
      is.table$name <- paste("feature", 1:nrow(is.table), sep = "_")
      
      is.table <-
        is.table %>%
        dplyr::select(name, mz, rt)
    }
    
    if (class(is.table)[1] == "try-error") {
      stop(crayon::red('Please provide right is table or mz and rt.\n'))
    }
    
    mz <-
      is.table %>%
      dplyr::pull(2)
    
    mz <- as.numeric(mz)
    
    mz_range <-
      lapply(mz, function(x) {
        c(x - ppm * x / 10 ^ 6, ppm * x / 10 ^ 6 + x)
      })
    
    mz_range <- do.call(rbind, mz_range)
    
    if (any(colnames(is.table) == "rt")) {
      rt <-
        is.table %>%
        dplyr::pull(3) %>%
        as.numeric()
      
      rt_range <-
        lapply(rt, function(x) {
          c(x - rt.tolerance, x + rt.tolerance)
        }) %>%
        do.call(rbind, .)
    } else{
      rt_range <- NA
    }
    
    cat(crayon::green("Extracting peaks, it will take a while..."))
    if (!is.na(rt_range)) {
      
      # raw_data <- xcms::filterRt(raw_data, c(min(rt_range), max(rt_range)))
      peak_data <- xcms::chromatogram(object = raw_data, 
                                      mz = mz_range,
                                      rt = rt_range,
                                      aggregationFun = "sum",
                                      msLevel = msLevel
                                      # filled = filled,
                                      # include = include
                                      )
      
    } else{
      peak_data <- xcms::chromatogram(object = raw_data,
                                      mz = mz_range,
                                      aggregationFun = "sum",
                                      msLevel = msLevel
                                      # filled = filled,
                                      # include = include
                                      )
    }
    cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
    
    save(peak_data, file = file.path(output_path, "peak_data"))
    return(peak_data)
  }
)







setGeneric(
  name = "showPeak2",
  def = function(object,
                 peak.index = 1,
                 title.size = 12,
                 lab.size = 12,
                 axis.text.size = 12,
                 alpha = 0.5,
                 title = "",
                 interactive = FALSE,
                 area_data,
                 raw_info) {
    options(warn = -1)
    info <- object@phenoData@data
    data <- object@.Data
    rm(list = c("object"))
    if (peak.index > nrow(data)) {
      peak.index <- nrow(data)
      cat("peak.index is ", nrow(data), '\n')
    }
    data <- apply(data, 2, function(x) {
      x <- x[[peak.index]]
      x <-
        data.frame(
          "rt" = x@rtime,
          "intensity" = x@intensity,
          stringsAsFactors = FALSE
        )
      list(x)
    })
    
    data <- lapply(data, function(x) {
      x[[1]]
    })
    
    data <- mapply(
      FUN = function(x, y, z) {
        x <- data.frame(
          x,
          "group" = y,
          "sample" = z,
          stringsAsFactors = FALSE
        )
        list(x)
      },
      x = data,
      y = info[, 2],
      z = info[, 1]
    )
    
    data <- do.call(rbind, args = data)
    
    ######
    ##check rt
    area_data <- 
      area_data %>% 
      apply(1, function(z){
        if(as.numeric(z[3]) > as.numeric(z[5]) |
           as.numeric(z[3]) < as.numeric(z[4])){
          z[3] <- mean(c(as.numeric(z[4]), as.numeric(z[5])))
        }
        z
      }) %>% 
      t() %>% 
      as.data.frame()
    
    area_data$rt <- as.numeric(area_data$rt)
    area_data$rtmin <- as.numeric(area_data$rtmin)
    area_data$rtmax <- as.numeric(area_data$rtmax)
    area_data$into <- as.numeric(area_data$into)
    area_data$intb <- as.numeric(area_data$intb)
    area_data$maxo <- as.numeric(area_data$maxo)
    area_data$sn <- as.numeric(area_data$sn)
    
    area_data <-
      area_data %>%
      dplyr::rename(rtmedian = rt) %>% 
      dplyr::arrange(sample) 
    
    area_data_raw_info <-
      purrr::map2(
        .x = area_data$sample,
        .y = raw_info[area_data$name[1], area_data$sample],
        .f = function(x, y) {
          if(length(y@productMz) > 2){
            fit_int <- as.numeric(y@productMz)
          }else{
            fit_int <- NA
          }
          data.frame(
            rt = y@rtime,
            int = y@intensity,
            fit_int,
            sample = x,
            stringsAsFactors = FALSE
          )
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    area_data_raw_info <- 
      area_data_raw_info %>%
      dplyr::left_join(area_data, by = "sample") %>% 
      dplyr::select(-c(rt:fit_int, name)) %>% 
      dplyr::distinct(sample, .keep_all = TRUE)
    

    final_data <-
      data %>% 
      dplyr::select(group, rt, intensity) %>% 
      dplyr::left_join(
        area_data_raw_info, by = c("group" = "sample")
      ) %>% 
      dplyr::mutate(group = paste(group,
                                   format(into, scientific = TRUE),
                                   sep = ":"))
    
    plot <-
      ggplot(final_data) +
      geom_line(aes(rt, intensity, color = group, group = group)) +
      geom_area(
        aes(rt, intensity, fill = group, group = group),
        alpha = 0.3,
        data = final_data %>%
          dplyr::filter(rt > rtmin & rt < rtmax)
      ) +
      geom_rect(
        aes(
          xmin = rtmin,
          xmax = rtmax,
          ymin = 0,
          ymax = maxo,
          col = group
        ),
        fill = "transparent",
        data = final_data[, c("group", "rtmin", "rtmax", "maxo")] %>%
          dplyr::distinct(group, .keep_all = TRUE)
      ) +
      # geom_line(aes(rt, fit_int, group = sample), color = "black") +
      geom_point(aes(rt, intensity),
                 shape = 16,
                 size = 0.5) +
      labs(x = "", y = "") +
      theme_bw() +
      facet_wrap(facets = vars(group)) +
      theme(
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1
        )
      )
    
    if (interactive) {
      plot <- plotly::ggplotly(plot)
    }
    return(plot)
  }
)




get_quantification_data2 <- function(path = ".",
                                     lipid_data,
                                     is_data,
                                     match_item) {
  
  lipid_data$rt <- 
    lipid_data$rt * 60
  # ###get Cholesterol from peak table and add it to lipid search table
  # temp_idx <-
  #   lapply(peak_table, function(x) {
  #     which(abs(x$mz - 369.3521) * 10 ^ 6 / 369.3521 <= 25 &
  #             abs(x$rt - is_table$rt[is_table$`Compound Name` == "Cholesterol (d7)"]) <= 90)
  #   })
  # 
  # if (length(unique(unlist(temp_idx))) >= 1) {
  #   temp_data <-
  #     purrr::map2(
  #       .x = temp_idx,
  #       .y = peak_table,
  #       .f = function(x, y) {
  #         x <-
  #           x[which.min(abs(y$rt[x] - is_table$rt[is_table$`Compound Name` == "Cholesterol (d7)"]))]
  #         
  #         temp_data <-
  #           y[x, , drop = FALSE] %>%
  #           dplyr::mutate(peak_name = name) %>%
  #           dplyr::mutate(name = "Cholesterol") %>%
  #           dplyr::mutate(Class = "Chol") %>%
  #           dplyr::mutate(rt = rt / 60) %>%
  #           dplyr::mutate(mean.int = 10000)
  #         
  #         temp_data <-
  #           temp_data %>%
  #           dplyr::select(one_of(colnames(data)))
  #         temp_data
  #       }
  #     )
  # 
  # 
  # 
  # data <-
  #   data %>%
  #   dplyr::full_join(new_data, by = intersect(colnames(data), colnames(new_data)))
  # 
  # data$adduct[is.na(data$adduct)] <- "H-H2O"
  # data$polarity[is.na(data$polarity)] <- "positive"
  # data$IonFormula[is.na(data$IonFormula)] <- "C27 H45 O0"
  # 
  
  ###calculate concentration for each samples and lipid--------------------------
  lipid_tag <- 
    lipid_data %>% 
    dplyr::select(peak_name : mean.int)
  
  lipid_table <- 
    lipid_data %>% 
    dplyr::select(-c(peak_name : mean.int))
  
  is_tag <- 
    is_data %>% 
    dplyr::select(name : `RT NEG (second)`)
  
  is_table <- 
    is_data %>% 
    dplyr::select(-c(name : `RT NEG (second)`))
  
  diff_name1 <- 
    setdiff(colnames(lipid_table), colnames(is_table))
  
  diff_name2 <- 
    setdiff(colnames(is_table), colnames(lipid_table))
  
  if(length(diff_name1) > 0){
    lipid_table <-
      lipid_table %>% 
      dplyr::select(-diff_name1)
  }
  
  if(length(diff_name2) > 0){
    is_table <-
      is_table %>% 
      dplyr::select(-diff_name2)
  }
  
  expression_data_abs <-
    cal_abs(
      lipid_tag = lipid_tag,
      lipid_table = lipid_table,
      is_tag = is_tag,
      is_table = is_table,
      match_item = match_item
    )
  
  express_data_abs_ug_ml <- expression_data_abs$express_data_abs_ug_ml
  express_data_abs_um <- expression_data_abs$express_data_abs_um
  variable_info_abs <- expression_data_abs$variable_info_abs
  
  rownames(express_data_abs_ug_ml) <-
    rownames(express_data_abs_um) <-
    variable_info_abs$peak_name
  
  save(variable_info_abs, file = file.path(path, "variable_info_abs"))
  save(express_data_abs_ug_ml, file = file.path(path, "express_data_abs_ug_ml"))
  save(express_data_abs_um, file = file.path(path, "express_data_abs_um"))
  save(lipid_tag, file = file.path(path, "lipid_tag"))
  save(lipid_table, file = file.path(path, "lipid_table"))
  save(is_tag, file = file.path(path, "is_tag"))
  save(is_table, file = file.path(path, "is_table"))
}





process_mv <- function(quantification_data,
                       na.tolerance = 0.8) {
  tag <-
    quantification_data %>%
    dplyr::select(c(name:adduct))
  
  data <-
    quantification_data %>%
    dplyr::select(-c(name:adduct))

  remain_idx <- 
  apply(data, 1, function(x){
    sum(is.na(x))/ncol(data)
  }) %>% 
    `<`(na.tolerance) %>% 
    which()

  tag <- tag[remain_idx,]  

  data <- data[remain_idx,]
  
  data <- 
    impute::impute.knn(data = as.matrix(data), k = 10)$data %>% 
    as.data.frame()
  return(cbind(tag, data))    
  }


