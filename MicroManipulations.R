# Working with full data sets sucks: make them smaller
# by grabbing ~300 rt scans and up to 80 m/z

library(MSnbase)
library(dplyr)
library(future.apply)

microfiles <- list.files("Z:/1_QEdata/Will/RectangulaRdata/microdata/", full.names = T)


# Read in small data ----

microblank <- readMSData(microfiles[1], mode = "inMemory", centroided. = F, msLevel. = 1)
micropool_1 <- readMSData(microfiles[2], mode = "inMemory", centroided. = F, msLevel. = 1)
micropool_2 <- readMSData(microfiles[3], mode = "inMemory", centroided. = F, msLevel. = 1)
micropool_3 <- readMSData(microfiles[4], mode = "inMemory", centroided. = F, msLevel. = 1)

par(mfcol=c(4,2))
par(mar=c(2,4,0,0.1))
mz_1 <- unlist(unname(mz(micropool_1)))
mz_2 <- unlist(unname(mz(micropool_2)))
mz_3 <- unlist(unname(mz(micropool_3)))
mz_raw <- c(mz_1, mz_2, mz_3)

massCheck <- function(mass){
  hist(mz_1[mz_1>mass&mz_1<(mass+0.001)],breaks = 100, main="")
  hist(mz_2[mz_2>mass&mz_2<(mass+0.001)],breaks = 100, main="")
  hist(mz_3[mz_3>mass&mz_3<(mass+0.001)],breaks = 100, main="")
  hist(mz_raw[mz_raw>mass&mz_raw<(mass+0.001)],breaks = 100, main="")
}
massCheck(60.044)
massCheck(117.078979+1.007267)
par(mfrow=c(1,1))



# Find medians ----
segment_size <- 7
lowest_mz <- floor(min(mz_raw))-floor(min(mz_raw))%%segment_size
highest_mz <- floor(max(mz_raw))+floor(max(mz_raw))%%segment_size
bin_minima <- seq(lowest_mz, highest_mz, segment_size)

getMeds <- function(bin_min, mz_raw){
  raw_i <- mz_raw[mz_raw>(bin_min)&mz_raw<(bin_min+segment_size)]
  span <- round(bin_min/segment_size)
  peak_mids <- list()
  for(j in seq(0, segment_size-0.1, 0.1)){
    raw_bin_j <- raw_i[raw_i>(bin_min+j)&raw_i<(bin_min+j+0.1)]
    if(length(raw_bin_j)>1){
      hist_bin_j <- hist(raw_bin_j, breaks = 10000, plot = F)
      peak_idxs <- MALDIquant:::.localMaxima(hist_bin_j$counts, span)
      peak_mids[[(j+0.1)*10]] <- hist_bin_j$mids[peak_idxs]
    }
  }
  return(peak_mids)
}

#plan(multiprocess, workers = min(availableCores()-1), 5)
#mids <- unlist(future_lapply(bin_minima, getMeds, mz_raw))
mids <- unlist(lapply(bin_minima, getMeds, mz_raw))
bin_halves <- mids[-length(mids)]+diff(mids)/2
bin_halves <- c(0, bin_halves, 1000)



# Check on medians ----
par(mfrow=c(1,2))
hist(mz_raw[mz_raw>60.044&mz_raw<60.045],breaks = 100, main="")
abline(v=mids, col=rgb(1,0,0,0.2))
abline(v=bin_halves, col=rgb(0,0,1,0.2))
hist(mz_raw[mz_raw>118.08&mz_raw<(118.09)],breaks = 1000, main="")
abline(v=mids, col=rgb(1,0,0,0.2))
abline(v=bin_halves, col=rgb(0,0,1,0.2))
par(mfrow=c(1,1))


will_plotXIC <- function (x, main = "", col = "grey", colramp = topo.colors, 
                          grid.color = "lightgrey", pch = 21, yl = NULL,
                          mn = NULL, ...) {
  x <- suppressWarnings(as(x, "data.frame"))
  bpi <- unlist(lapply(split(x$i, x$rt), max, na.rm = TRUE))
  brks <- lattice::do.breaks(range(x$i), nint = 256)
  par(mar = c(0, 4.5, 2, 1))
  plot(as.numeric(names(bpi)), bpi, xaxt = "n", col = col, 
       main = mn, bg = lattice:::level.colors(bpi, at = brks, col.regions = colramp), 
       xlab = "", pch = pch, ylab = "", las = 2, 
       ...)
  mtext(side = 4, line = 0, "Intensity", cex = par("cex.lab"))
  grid(col = grid.color)
  par(mar = c(3.5, 4.5, 0, 1))
  plot(x$rt, x$mz, main = "", pch = pch, col = col, xlab = "", 
       ylab = "", yaxt = "n", ylim = yl,
       bg = lattice:::level.colors(x$i, at = brks, col.regions = colramp), ...)
  axis(side = 2, las = 2)
  grid(col = grid.color)
  mtext(side = 1, line = 2.5, "Retention time", cex = par("cex.lab"))
  mtext(side = 4, line = 0, "m/z", cex = par("cex.lab"))
}
par(mfcol=c(2,4))
micropool_1 %>%
  filterMz(c(118.085, 118.089)) %>%
  will_plotXIC(col=NA, yl=c(118.085, 118.089))
abline(h=bin_halves, col="red")
micropool_2 %>%
  filterMz(c(118.085, 118.089)) %>%
  will_plotXIC(col=NA, yl=c(118.085, 118.089))
abline(h=bin_halves, col="red")
micropool_3 %>%
  filterMz(c(118.085, 118.089)) %>%
  will_plotXIC(col=NA, yl=c(118.085, 118.089))
abline(h=bin_halves, col="red")
microblank %>%
  filterMz(c(118.085, 118.089)) %>%
  will_plotXIC(col=NA, yl=c(118.085, 118.089))
abline(h=bin_halves, col="red")

par(mfcol=c(2,4))
micropool_1 %>%
  filterMz(c(60.044, 60.048)) %>%
  will_plotXIC(col=NA, yl=c(60.044, 60.048))
abline(h=bin_halves, col="red")
micropool_2 %>%
  filterMz(c(60.044, 60.048)) %>%
  will_plotXIC(col=NA, yl=c(60.044, 60.048))
abline(h=bin_halves, col="red")
micropool_3 %>%
  filterMz(c(60.044, 60.048)) %>%
  will_plotXIC(col=NA, yl=c(60.044, 60.048))
abline(h=bin_halves, col="red")
microblank %>%
  filterMz(c(60.044, 60.048)) %>%
  will_plotXIC(col=NA, yl=c(60.044, 60.048))
abline(h=bin_halves, col="red")



# Apply medians ----
middling <- function(scan){
  scan_mz <- scan@mz
  scan_cut <- cut(scan_mz, bin_halves)
  new_mz <- mids[scan_cut]
  scan@mz <- new_mz
  return(scan)
}

centered_data <- future_eapply(raw_data_mem@assayData, FUN = middling)
mod_data <- raw_data_mem
mod_data@assayData <- as.environment(centered_data)
