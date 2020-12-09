##avoid source
no_function()

sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")

file_path <- "lipid20201020"
setwd(file.path("data", file_path))

load("lipid_search/POS/is_table_pos")
load("lipid_search/POS/is_tag_pos")
load("lipid_search/POS/lipid_table_pos")
load("lipid_search/POS/lipid_tag_pos")
load("lipid_search/POS/express_data_abs_pos1")
load("lipid_search/POS/express_data_abs_pos2")
load("lipid_search/POS/variable_info_abs_pos")

load("lipid_search/NEG/is_table_neg")
load("lipid_search/NEG/is_tag_neg")
load("lipid_search/NEG/lipid_table_neg")
load("lipid_search/NEG/lipid_tag_neg")
load("lipid_search/NEG/express_data_abs_neg1")
load("lipid_search/NEG/express_data_abs_neg2")
load("lipid_search/NEG/variable_info_abs_neg")


#####
sxtTools::setwd_project()
setwd(file.path("data", file_path, "Result"))

lipid_data <-
  readr::read_csv("lipid_data_um.csv")

rownames(express_data_abs_pos1) <-
  rownames(express_data_abs_pos2) <-
  variable_info_abs_pos$peak_name

variable_info_abs_pos <-
  variable_info_abs_pos %>%
  dplyr::filter(peak_name %in% lipid_data$peak_name)

express_data_abs_pos1 <-
  express_data_abs_pos1[match(variable_info_abs_pos$peak_name,
                              rownames(express_data_abs_pos1)),]

express_data_abs_pos2 <-
  express_data_abs_pos2[match(variable_info_abs_pos$peak_name,
                              rownames(express_data_abs_pos2)),]

rownames(lipid_table_pos) <- lipid_tag_pos$peak_name

lipid_tag_pos <-
  lipid_tag_pos %>%
  dplyr::filter(peak_name %in% lipid_data$peak_name)

lipid_table_pos <-
  lipid_table_pos[match(lipid_tag_pos$peak_name, rownames(lipid_table_pos)),]


rownames(express_data_abs_neg1) <-
  rownames(express_data_abs_neg2) <-
  variable_info_abs_neg$peak_name

variable_info_abs_neg <-
  variable_info_abs_neg %>%
  dplyr::filter(peak_name %in% lipid_data$peak_name)

express_data_abs_neg1 <-
  express_data_abs_neg1[match(variable_info_abs_neg$peak_name,
                              rownames(express_data_abs_neg1)),]

express_data_abs_neg2 <-
  express_data_abs_neg2[match(variable_info_abs_neg$peak_name,
                              rownames(express_data_abs_neg2)),]

rownames(lipid_table_neg) <- lipid_tag_neg$peak_name

lipid_tag_neg <-
  lipid_tag_neg %>%
  dplyr::filter(peak_name %in% lipid_data$peak_name)

lipid_table_neg <-
  lipid_table_neg[match(lipid_tag_neg$peak_name, rownames(lipid_table_neg)),]


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


sample_info <-
  data.frame(sample_id = c("Blank",
                           paste("test_sample", 1:4, sep = "")),
             class = c(rep("test", 5)))


dir.create(file.path("intensity_plot"),
           showWarnings = FALSE)


###output
###positive
for (i in names(match_item_pos)) {
  is_name <-
    match_item_pos[[i]]
  
  if (is.na(unique(is_name)[1])) {
    next()
  }
  
  dir.create(file.path("intensity_plot", i),
             showWarnings = FALSE)
  
  temp_idx <-
    which(lipid_tag_pos$Class == i)
  
  if (length(temp_idx) > 0) {
    for (j in temp_idx) {
      temp_is_name <-
        abs(lipid_tag_pos$rt[j] - is_tag_pos$rt[match(is_name, is_tag_pos$`Compound Name`)]) %>%
        which.min() %>%
        head(1) %>%
        `[`(is_name, .)
      
      temp_data1 <-
        is_table_pos[match(temp_is_name, is_tag_pos$`Compound Name`),] %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::left_join(sample_info, by = "sample_id")
      
      colnames(temp_data1)[2] <- "intensity"
      
      ###raw data
      temp_name <- lipid_tag_pos$name[j]
      temp_data2 <-
        lipid_table_pos[j, ] %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::left_join(sample_info, by = "sample_id")
      
      colnames(temp_data2)[2] <- "intensity"
      
      
      temp_data3 <-
        express_data_abs_pos2[match(temp_name, variable_info_abs_pos$name), ] %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::left_join(sample_info, by = "sample_id")
      
      colnames(temp_data3)[2] <- "intensity"
      
      temp_data <-
        rbind(
          data.frame(temp_data1, from = "is"),
          data.frame(temp_data2, from = "raw"),
          data.frame(temp_data3, from = "final")
        ) %>%
        dplyr::mutate(from = factor(from, levels = c("is", "raw", "final")))
      
      temp_plot <-
        temp_data %>%
        ggplot(aes(sample_id, intensity)) +
        geom_point(size = 4, shape = 21, aes(fill = class)) +
        ggrepel::geom_text_repel(aes(label = round(intensity, 4))) +
        ggsci::scale_fill_d3() +
        scale_y_continuous(expand = expansion(mult = c(0.5, 0.5))) +
        theme_bw() +
        labs(x = "", y = "Intensity") +
        theme() +
        facet_wrap(facets = vars(from),
                   ncol = 1,
                   scales = "free_y")
      
      temp_name <- paste(temp_name, ".pdf", sep = "")
      temp_name <-
        stringr::str_replace_all(temp_name, "\\/", "_")
      ggsave(
        temp_plot,
        filename = file.path("intensity_plot", i, temp_name),
        width = 7,
        height = 14
      )
    }
  }
}


###negative
for (i in names(match_item_neg)) {
  is_name <-
    match_item_neg[[i]]
  
  if (is.na(unique(is_name)[1])) {
    next()
  }
  
  dir.create(file.path("intensity_plot", i),
             showWarnings = FALSE)
  
  temp_idx <-
    which(lipid_tag_neg$Class == i)
  
  if (length(temp_idx) > 0) {
    for (j in temp_idx) {
      temp_is_name <-
        abs(lipid_tag_neg$rt[j] - is_tag_neg$rt[match(is_name, is_tag_neg$`Compound Name`)]) %>%
        which.min() %>%
        head(1) %>%
        `[`(is_name, .)
      
      temp_data1 <-
        is_table_neg[match(temp_is_name, is_tag_neg$`Compound Name`),] %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::left_join(sample_info, by = "sample_id")
      
      colnames(temp_data1)[2] <- "intensity"
      
      ###raw data
      temp_name <- lipid_tag_neg$name[j]
      temp_data2 <-
        lipid_table_neg[j, ] %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::left_join(sample_info, by = "sample_id")
      
      colnames(temp_data2)[2] <- "intensity"
      
      
      temp_data3 <-
        express_data_abs_neg2[match(temp_name, variable_info_abs_neg$name), ] %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "sample_id") %>%
        dplyr::left_join(sample_info, by = "sample_id")
      
      colnames(temp_data3)[2] <- "intensity"
      
      temp_data <-
        rbind(
          data.frame(temp_data1, from = "is"),
          data.frame(temp_data2, from = "raw"),
          data.frame(temp_data3, from = "final")
        ) %>%
        dplyr::mutate(from = factor(from, levels = c("is", "raw", "final")))
      
      temp_plot <-
        temp_data %>%
        ggplot(aes(sample_id, intensity)) +
        geom_point(size = 4, shape = 21, aes(fill = class)) +
        ggrepel::geom_text_repel(aes(label = round(intensity, 4))) +
        ggsci::scale_fill_d3() +
        scale_y_continuous(expand = expansion(mult = c(0.5, 0.5))) +
        theme_bw() +
        labs(x = "", y = "Intensity") +
        theme() +
        facet_wrap(facets = vars(from),
                   ncol = 1,
                   scales = "free_y")
      
      temp_name <- paste(temp_name, ".pdf", sep = "")
      temp_name <-
        stringr::str_replace_all(temp_name, "\\/", "_")
      ggsave(
        temp_plot,
        filename = file.path("intensity_plot", i, temp_name),
        width = 7,
        height = 14
      )
    }
  }
}








#####Class plot
dir.create("class_plot")
####each class
lipid_data_class <- readr::read_csv("lipid_data_class_um.csv")
for (i in 1:nrow(lipid_data_class)) {
  cat(i, " ")
  temp_data <-
    lipid_data_class[i, -1] %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample_id") %>%
    dplyr::left_join(sample_info, by = "sample_id")
  colnames(temp_data)[2] <- "intensity"
  
  temp_plot <-
    temp_data %>%
    ggplot(aes(sample_id, intensity)) +
    geom_point(size = 4, shape = 21, aes(fill = class)) +
    ggrepel::geom_text_repel(aes(label = round(intensity, 4))) +
    ggsci::scale_fill_d3() +
    scale_y_continuous(expand = expansion(mult = c(0.5, 0.5))) +
    theme_bw() +
    labs(x = "", y = "Intensity") +
    theme()
  
  ggsave(
    temp_plot,
    filename = file.path("class_plot",
                         paste(lipid_data_class$X1[i], ".pdf",
                               sep = "")),
    width = 7,
    height = 7
  )
  
}
