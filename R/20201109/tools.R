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
  
  ####get the final match_itme and final is_table and is_table
  intersect_name <-
    intersect(is_tag$`Compound Name`, unique(unlist(match_item)))
  remain_idx <- which(is_tag$`Compound Name` %in% intersect_name)
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
  
  express_data_abs1 <-
    express_data_abs2 <-
    lipid_table
  
  for (i in 1:nrow(lipid_tag)) {
    cat(i, " ")
    temp_class <- lipid_tag$Class[i]
    temp_idx <- match_item[[temp_class]] %>%
      match(., is_tag$`Compound Name`)
    temp_rt <- lipid_tag$rt[i]
    
    temp_idx <-
      temp_idx[which.min(abs(temp_rt - is_tag$rt[temp_idx]))]
    
    ratio <-
      unlist(lipid_table[i, , drop = TRUE]) / unlist(is_table[temp_idx, , drop = TRUE])
    
    ratio[is.na(ratio)] <- 0
    ratio[is.infinite(ratio)] <- max(ratio[!is.infinite(ratio)])
    
    concentration1 <-
      is_tag$`μg/mL`[temp_idx] * ratio
    
    concentration2 <-
      is_tag$μM[temp_idx] * ratio
    
    express_data_abs1[i, ] <- concentration1
    express_data_abs2[i, ] <- concentration2
  }
  
  
  return_result <- list(
    express_data_abs1 = express_data_abs1,
    express_data_abs2 = express_data_abs2,
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
  
  express_data_abs_pos1 <- expresion_data_abs_pos$express_data_abs1
  express_data_abs_pos2 <- expresion_data_abs_pos$express_data_abs2
  variable_info_abs_pos <- expresion_data_abs_pos$variable_info_abs
  
  save(variable_info_abs_pos, file = "lipid_search/POS/variable_info_abs_pos")
  save(express_data_abs_pos1, file = "lipid_search/POS/express_data_abs_pos1")
  save(express_data_abs_pos2, file = "lipid_search/POS/express_data_abs_pos2")
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
  
  express_data_abs_neg1 <- expresion_data_abs_neg$express_data_abs1
  express_data_abs_neg2 <- expresion_data_abs_neg$express_data_abs2
  variable_info_abs_neg <- expresion_data_abs_neg$variable_info_abs
  
  save(variable_info_abs_neg, file = "lipid_search/NEG/variable_info_abs_neg")
  save(express_data_abs_neg1, file = "lipid_search/NEG/express_data_abs_neg1")
  save(express_data_abs_neg2, file = "lipid_search/NEG/express_data_abs_neg2")
  
  save(lipid_tag_neg, file = "lipid_search/NEG/lipid_tag_neg")
  save(lipid_table_neg, file = "lipid_search/NEG/lipid_table_neg")
  save(is_tag_neg, file = "lipid_search/NEG/is_tag_neg")
  save(is_table_neg, file = "lipid_search/NEG/is_table_neg")
}



combine_pos_neg_quantification <- function(path = ".",
                                           express_data_abs_pos1,
                                           express_data_abs_pos2,
                                           variable_info_abs_pos,
                                           express_data_abs_neg1,
                                           express_data_abs_neg2,
                                           variable_info_abs_neg) {
  dir.create(path, showWarnings = FALSE)
  # ####remove duplicated
  # ###positive mode
  # library(plyr)
  # express_data_abs_pos1 <-
  #   express_data_abs_pos1 %>%
  #   data.frame(variable_info_abs_pos, ., check.names = FALSE) %>%
  #   plyr::dlply(.variables = .(name)) %>%
  #   purrr::map(
  #     .f = function(x) {
  #       x <-
  #         x %>%
  #         dplyr::filter(mean.int == max(mean.int)) %>%
  #         head(1)
  #       # dplyr::select(colnames(express_data_abs1))
  #     }
  #   ) %>%
  #   do.call(rbind, .) %>%
  #   as.data.frame()
  # 
  # rownames(express_data_abs_pos1) <- NULL
  # 
  # express_data_abs_pos1 <-
  #   express_data_abs_pos1 %>%
  #   tibble::column_to_rownames(var = "peak_name") %>%
  #   dplyr::select(-colnames(variable_info_abs_pos)[-1])
  # 
  # express_data_abs_pos2 <-
  #   express_data_abs_pos2 %>%
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
  # rownames(express_data_abs_pos2)   <- NULL
  # 
  # express_data_abs_pos2 <-
  #   express_data_abs_pos2 %>%
  #   tibble::column_to_rownames(var = "peak_name") %>%
  #   dplyr::select(-colnames(variable_info_abs_pos)[-1])
  # 
  # variable_info_abs_pos <-
  #   variable_info_abs_pos[match(rownames(express_data_abs_pos1), variable_info_abs_pos$peak_name),]
  # 
  # ###negative mode
  # library(plyr)
  # express_data_abs_neg1 <-
  #   express_data_abs_neg1 %>%
  #   data.frame(variable_info_abs_neg, ., check.names = FALSE) %>%
  #   plyr::dlply(.variables = .(name)) %>%
  #   purrr::map(
  #     .f = function(x) {
  #       x <-
  #         x %>%
  #         dplyr::filter(mean.int == max(mean.int)) %>%
  #         head(1)
  #       # dplyr::select(colnames(express_data_abs1))
  #     }
  #   ) %>%
  #   do.call(rbind, .) %>%
  #   as.data.frame()
  # 
  # rownames(express_data_abs_neg1) <- NULL
  # 
  # express_data_abs_neg1 <-
  #   express_data_abs_neg1 %>%
  #   tibble::column_to_rownames(var = "peak_name") %>%
  #   dplyr::select(-colnames(variable_info_abs_neg)[-1])
  # 
  # express_data_abs_neg2 <-
  #   express_data_abs_neg2 %>%
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
  # rownames(express_data_abs_neg2)   <- NULL
  # 
  # express_data_abs_neg2 <-
  #   express_data_abs_neg2 %>%
  #   tibble::column_to_rownames(var = "peak_name") %>%
  #   dplyr::select(-colnames(variable_info_abs_neg)[-1])
  # 
  # variable_info_abs_neg <-
  #   variable_info_abs_neg[match(rownames(express_data_abs_neg1), variable_info_abs_neg$peak_name),]
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
  #   express_data_abs_pos1 <- 
  #     express_data_abs_pos1[-remove_idx,]
  #   express_data_abs_pos2 <- 
  #     express_data_abs_pos2[-remove_idx,]
  # }
  
  ##combine positive and negative
  rownames(express_data_abs_pos1) <-
    rownames(express_data_abs_pos2) <-
    variable_info_abs_pos$peak_name
  
  rownames(express_data_abs_neg1) <-
    rownames(express_data_abs_neg2) <-
    variable_info_abs_neg$peak_name
  
  express_data_abs1 <-
    rbind(express_data_abs_pos1,
          express_data_abs_neg1)
  
  colnames(express_data_abs_pos2) <-
    colnames(express_data_abs_neg2)
  
  express_data_abs2 <-
    rbind(express_data_abs_pos2,
          express_data_abs_neg2)
  
  variable_info_abs <-
    variable_info_abs_pos %>%
    dplyr::full_join(variable_info_abs_neg, by = intersect(
      colnames(variable_info_abs_pos),
      colnames(variable_info_abs_neg)
    ))
  
  ####remove duplicated
  library(plyr)
  express_data_abs1 <-
    express_data_abs1 %>%
    data.frame(variable_info_abs, ., check.names = FALSE) %>%
    plyr::dlply(.variables = .(name)) %>%
    purrr::map(
      .f = function(x) {
        x <-
          x %>%
          dplyr::filter(mean.int == max(mean.int)) %>%
          head(1)
        # dplyr::select(colnames(express_data_abs1))
      }
    ) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  rownames(express_data_abs1)  <- NULL
  
  express_data_abs1 <-
    express_data_abs1 %>%
    tibble::column_to_rownames(var = "peak_name") %>%
    dplyr::select(-colnames(variable_info_abs)[-1])
  
  library(plyr)
  express_data_abs2 <-
    express_data_abs2 %>%
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
  
  rownames(express_data_abs2)   <- NULL
  
  express_data_abs2 <-
    express_data_abs2 %>%
    tibble::column_to_rownames(var = "peak_name") %>%
    dplyr::select(-colnames(variable_info_abs)[-1])
  
  variable_info_abs <-
    variable_info_abs[match(rownames(express_data_abs1), variable_info_abs$peak_name),]
  
  library(plyr)
  class_data1 <-
    express_data_abs1 %>%
    data.frame(
      class = variable_info_abs$Class,
      .,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) %>%
    plyr::dlply(.variables = .(class)) %>%
    purrr::map(
      .f = function(x) {
        apply(x[, -1], 2, sum)
      }
    ) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  write.csv(class_data1, file.path(path, "class_data1"))
  
  
  library(plyr)
  class_data2 <-
    express_data_abs2 %>%
    data.frame(
      class = variable_info_abs$Class,
      .,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) %>%
    plyr::dlply(.variables = .(class)) %>%
    purrr::map(
      .f = function(x) {
        apply(x[, -1], 2, sum)
      }
    ) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  write.csv(class_data2, file.path(path, "class_data2"))
  
  plot1 <-
    class_data1 %>%
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
    width = 7,
    height = 7
  )
  
  
  plot2 <-
    class_data2 %>%
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
    width = 7,
    height = 7
  )
  
  express_data_abs2_per <-
    express_data_abs2 %>%
    apply(2, function(x) {
      x * 100 / sum(x)
    }) %>%
    as.data.frame()
  
  class_data2_per <-
    class_data2 %>%
    apply(2, function(x) {
      x * 100 / sum(x)
    }) %>%
    as.data.frame()
  
  write.csv(
    cbind(variable_info_abs, express_data_abs1),
    file.path(path, "lipid_data_ug_ml.csv"),
    row.names = FALSE
  )
  
  write.csv(
    cbind(variable_info_abs, express_data_abs2),
    file.path(path, "lipid_data_um.csv"),
    row.names = FALSE
  )
  
  write.csv(
    cbind(variable_info_abs, express_data_abs2_per),
    file.path(path, "lipid_data_um_per.csv"),
    row.names = FALSE
  )
  
  write.csv(class_data1,
            file.path(path, "lipid_data_class_ug_ml.csv"),
            row.names = TRUE)
  
  write.csv(class_data2,
            file.path(path, "lipid_data_class_um.csv"),
            row.names = TRUE)
  
  write.csv(class_data2_per,
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











get_lipid_quantify_table <- function(lipid_table,
                                     peak_table,
                                     # sample_info,
                                     mz.tol = 25,
                                     rt.tol = 90,
                                     figure = FALSE) {
  data_is <- lipid_table[, c("mz", "rt")]
  
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
          nrow(lipid_table),
          "internal standards are found in peak table.\n"
        )
        
        is_tag <-
          lipid_table[match_result$Index1,]
        
        is_data <- x[match_result$Index2, ]  %>%
          dplyr::select(-c(name, mz, rt))
        
        is_data <-
          cbind(is_tag, is_data)
        is_data
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
           sample_numer_show = 3,
           fit.gaussian = TRUE,
           integrate_xcms = FALSE) {
    peak_table <- readxl::read_xlsx(targeted_targeted_peak_table_name)
    library(tidyverse)
    result_raw <-
      metflow2::extractPeaks(
        path = path,
        ppm = 40,
        threads = 3,
        is.table = targeted_targeted_peak_table_name,
        rt.tolerance = 90
      )
    
    # load("raw_data")
    # raw_data %>%
    #   filterRt(rt = as.numeric(result_raw@featureData@data[1,3:4])) %>%
    #   filterMz(mz = as.numeric(result_raw@featureData@data[1,1:2])) %>%
    #   plot(type = "XIC")
    
    library(xcms)
    result <- findChromPeaks(object = result_raw,
                             param = CentWaveParam(peakwidth = c(5, 80),
                                                   snthresh = 0.1, 
                                                   ppm = 40, 
                                                   fitgauss = TRUE, 
                                                   # prefilter = c(3, 500),
                                                   mzdiff = 0.01))
    
    # plot(result)
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
    
    
    output_quantification_table <- as.data.frame(output_quantification_table)
    output_quantification_table_old <- output_quantification_table
    
    if(!integrate_xcms){
      output_quantification_table[,-c(1:2)] <- NA
    }
    
    ###check missing value
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
            data.frame(x = temp_raw_info@rtime, y = temp_raw_info@intensity, stringsAsFactors = FALSE) %>% 
            dplyr::filter(!is.na(y))
          
          if(nrow(xy) <= 5){
            value <- 0
          }else{
            fit_result <-
              fit_gaussian(x = xy$x, y = xy$y)
            xy <-
              data.frame(xy, z = fit_result$y)
            
            # plot(xy$x, xy$y, type = "b")
            # points(xy$x, xy$z, type = "b", col = "red")
            # 
            # abline(v = fit_result$center)
            # abline(v = fit_result$center - fit_result$width * 2)
            # abline(v = fit_result$center + fit_result$width * 2)
            # abline(h = 0)
            
            l1 <- loess(y ~ x, xy, control = loess.control(surface = "direct"))
            l2 <-
              loess(z ~ x, xy, control = loess.control(surface = "direct"))
            # f1 <- function(x)
            #   predict(l1, newdata = x)
            f1 <- approxfun(xy$x, xy$y)
            # f2 <- function(x)
            #   predict(l2, newdata = x)
            f2 <- approxfun(xy$x, xy$z)
            
            value1 <- try(expr = integrate(f = f1,
                                           lower = min(xy$x),
                                           upper = max(xy$x))$value,
                          silent = TRUE
            )
            
            value2 <- try(integrate(f = f2,
                                    lower = min(xy$x),
                                    upper = max(xy$x))$value, silent = TRUE
            )
            
            if(class(value1) == "try-error"){
              value1 <- NA
            }
            
            if(class(value2) == "try-error"){
              value2 <- NA
            }
            
            
            if(fit.gaussian){
              value <- value2
            }else{
              value <- value1
            }
            
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
          }
          output_quantification_table[output_quantification_table$name == temp_name, temp_sample] <-
            value
        }
      }
    }
    
    ###add mz and adduct
    output_quantification_table <-
      output_quantification_table %>% 
      dplyr::left_join(peak_table[,c("name", "mz", "adduct")], by = "name") %>% 
        dplyr::select(name, mz, rt, adduct, everything())
    
    openxlsx::write.xlsx(x = output_quantification_table,
                         file = file.path(path, "quantification_table.xlsx"),
                         asTable = TRUE)
    
    peak_value <-
      purrr::map2(
        as.data.frame(raw_info),
        .y = colnames(as.data.frame(raw_info)),
        .f = function(x, sample_name) {
          temp <- 
            purrr::map2(x, .y = names(x), .f = function(y, peak_name){
              if(nrow(as.data.frame(y@chromPeaks)) == 0){
                peak_name <- logical(0)
              }
              data.frame(name = peak_name,
                         as.data.frame(y@chromPeaks),
                         stringsAsFactors = FALSE)
            })  %>% 
            do.call(rbind, .) %>% 
            as.data.frame()
          
          if(nrow(temp) == 0){
            sample_name <- logical(0)
          }
          
          data.frame(sample = sample_name, temp, stringsAsFactors = FALSE)
        }
      ) %>% 
      do.call(rbind, .) %>% 
      as.data.frame()
    
    rownames(peak_value) <- NULL
    
    ###output peak shape
    dir.create(file.path(path, "peak_shape"))
    purrr::walk(
      .x = unique(peak_value$name),
      .f = function(x) {
        temp <-
          peak_value %>%
          dplyr::filter(name == x)
        
        temp <-
          temp %>%
          dplyr::rename(rtmedian = rt) %>% 
          dplyr::arrange(sample)
        
        if (nrow(temp) > sample_numer_show) {
          temp <- temp[1:sample_numer_show, ]
        }
        
        temp_raw_info <-
          raw_info[x, temp$sample] %>%
          purrr::map2(
            .x = temp$sample,
            .y = .,
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
        
        plot <-
          temp_raw_info %>%
          dplyr::left_join(temp, by = "sample") %>%
          dplyr::filter(rt >= min(rtmin) - 10 &
                          rt <= max(rtmax) + 10) %>%
          ggplot() +
          geom_line(aes(rt, int, color = sample, group = sample)) +
          geom_line(aes(rt, fit_int, group = sample), color = "black") +
          geom_point(aes(rt, int, color = sample),
                     shape = 16,
                     size = 2) +
          ggrepel::geom_text_repel(aes(
            x = rtmedian,
            y = maxo,
            color = sample,
            label = round(into)
          ),
          data = temp) +
          geom_rect(
            mapping = aes(
              xmin = rtmin,
              xmax = rtmax,
              ymin = 0,
              ymax = maxo,
              color = sample
            ),
            # fill = "transparent",
            fill = NA,
            alpha = 0.8,
            data = temp
          ) +
          labs(x = "RT (second)", y = "Intensity") +
          ggsci::scale_color_futurama() +
          ggsci::scale_fill_futurama() +
          theme_bw() +
          theme(panel.grid = element_blank(), legend.position = "bottom")
        
        name <- paste(x, "_peak_shape.pdf", sep = "")
        ggsave(
          plot = plot,
          filename = file.path(path, "peak_shape", name),
          width = 8,
          height = 7
        )
      }
    )
    
    
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

fit_gaussian <-
  function (x,
            y,
            start.center = NULL,
            start.width = NULL,
            start.height = NULL,
            start.floor = NULL,
            fit.floor = FALSE)
  {
    who.max <- which.max(y)
    if (is.null(start.center))
      start.center <- x[who.max]
    if (is.null(start.height))
      start.height <- y[who.max]
    if (is.null(start.width))
      start.width <- sum(y > (start.height / 2), na.rm = TRUE) / 2
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
      if (is.null(start.floor))
        start.floor <- quantile(y, seq(0, 1, 0.1))[2]
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
      yAns <- gaussian(x, centerAns, widthAns, heightAns, floorAns)
      residualAns <- y - yAns
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