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
           fit.gaussian = FALSE,
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