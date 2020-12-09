##
no_source()

sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")
library(tidyverse)
setwd("data/lipid20201020/xcms/POS/")

library(metflow2)
library(xcms)
library(MSnbase)

metflow2::processData(
  path = ".",
  polarity = "positive",
  peakwidth = c(5, 60),
  threads = 5,
  output.tic = FALSE,
  output.bpc = FALSE,
  min.fraction = 0.5
)



sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")
library(tidyverse)
setwd("data/lipid20201020/xcms/NEG/")

library(metflow2)
library(xcms)
library(MSnbase)

metflow2::processData(
  path = ".",
  polarity = "negative",
  peakwidth = c(5, 60),
  threads = 5,
  output.tic = FALSE,
  output.bpc = FALSE,
  min.fraction = 0.5
)
