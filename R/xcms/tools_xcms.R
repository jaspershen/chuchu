mzCenter.wMean <- function(mz,intensity) {
  weighted.mean(mz, intensity)
}

mzCenter.mean <- function(mz,intensity) {
  mean(mz)
}

mzCenter.apex <- function(mz,intensity) {
  mz[which.max(intensity)]
}




estimateChromNoise <- function(x, trim = 0.05, minPts = 20) {
  gz <- which(x > 0)
  if (length(gz) < minPts)
    return(mean(x))
  
  mean(x[gz], trim = trim)
}


continuousPtsAboveThreshold <-
  function(y, threshold, num, istart = 1) {
    if (!is.double(y))
      y <- as.double(y)
    if (.C(
      "continuousPtsAboveThreshold",
      y,
      as.integer(istart - 1),
      length(y),
      threshold = as.double(threshold),
      num = as.integer(num),
      n = integer(1),
      PACKAGE = "xcms"
    )$n > 0)
      TRUE
    else
      FALSE
  }



getLocalNoiseEstimate <- function(d, td, ftd, noiserange, Nscantime, threshold, num) {
  
  if (length(d) < Nscantime) {
    
    ## noiserange[2] is full d-range
    drange <- which(td %in% ftd)
    n1 <- d[-drange] ## region outside the detected ROI (wide)
    n1.cp <- continuousPtsAboveThresholdIdx(n1, threshold=threshold,num=num) ## continousPtsAboveThreshold (probably peak) are subtracted from data for local noise estimation
    n1 <- n1[!n1.cp]
    if (length(n1) > 1)  {
      baseline1 <- mean(n1)
      sdnoise1 <- sd(n1)
    } else
      baseline1 <- sdnoise1 <- 1
    
    ## noiserange[1]
    d1 <- drange[1]
    d2 <- drange[length(drange)]
    nrange2 <- c(max(1,d1 - noiserange[1]) : d1, d2 : min(length(d),d2 + noiserange[1]))
    n2 <- d[nrange2] ## region outside the detected ROI (narrow)
    n2.cp <- continuousPtsAboveThresholdIdx(n2, threshold=threshold,num=num) ## continousPtsAboveThreshold (probably peak) are subtracted from data for local noise estimation
    n2 <- n2[!n2.cp]
    if (length(n2) > 1)  {
      baseline2 <- mean(n2)
      sdnoise2 <- sd(n2)
    } else
      baseline2 <- sdnoise2 <- 1
    
  } else {
    trimmed <- trimm(d,c(0.05,0.95))
    baseline1 <- baseline2 <- mean(trimmed)
    sdnoise1 <- sdnoise2 <- sd(trimmed)
  }
  
  c(min(baseline1,baseline2),min(sdnoise1,sdnoise2))
}