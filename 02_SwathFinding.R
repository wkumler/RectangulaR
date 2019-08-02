# Normalize mz swaths to the median value
# Read in MS data, find medians, normalize to medians, and write out MS data

# Setup things ----

library(dplyr)
library(future.apply)
library(MSnbase)

pos_data_path <- "Z:/1_QEdata/Will/RectangulaRdata/positive/"
if(!dir.exists(pos_data_path)){stop("Is the network drive mapped?")}

blank_file_name <- "190715_Blk_KM1906U14-Blk_C.mzML"
std_mix1_h2o_name <- "190715_Std_4uMStdsMix1InH2O_1.mzML"
std_mix2_h2o_name <- "190715_Std_4uMStdsMix2InH2O_1.mzML"
std_mix1_mat_name <- "190715_Std_4uMStdsMix1InMatrix_1.mzML"
std_mix2_mat_name <- "190715_Std_4uMStdsMix2InMatrix_1.mzML"
pool_1_name <- "190715_Poo_TruePooFK180310_Full1.mzML"
pool_2_name <- "190715_Poo_TruePooFK180310_Full2.mzML"
pool_3_name <- "190715_Poo_TruePooFK180310_Full3.mzML"



# Load data ----

pool_1_data <- readMSData(files = paste0(pos_data_path, pool_1_name), mode = "inMemory", 
                           centroided. = F, msLevel. = 1)
pool_2_data <- readMSData(files = paste0(pos_data_path, pool_2_name), mode = "inMemory", 
                          centroided. = F, msLevel. = 1)
pool_3_data <- readMSData(files = paste0(pos_data_path, pool_3_name), mode = "inMemory", 
                          centroided. = F, msLevel. = 1)
blank_data <- readMSData(files = paste0(pos_data_path, blank_file_name), mode = "inMemory", 
                         centroided. = F, msLevel. = 1)

par(mfrow=c(4,1))
mz_1 <- unlist(unname(mz(pool_1_data)))
hist(mz_1[mz_1>60.044&mz_1<60.045],breaks = 100, main="")
mz_2 <- unlist(unname(mz(pool_2_data)))
hist(mz_2[mz_2>60.044&mz_2<60.045],breaks = 100, main="")
mz_3 <- unlist(unname(mz(pool_3_data)))
hist(mz_3[mz_3>60.044&mz_3<60.045],breaks = 100, main="")
mz_raw <- c(mz_1, mz_2, mz_3)
hist(mz_raw[mz_raw>60.044&mz_raw<60.045],breaks = 100, main="")

par(mfrow=c(4,1))
hist(mz_1[mz_1>200.044&mz_1<200.045],breaks = 100, main="")
hist(mz_2[mz_2>200.044&mz_2<200.045],breaks = 100, main="")
hist(mz_3[mz_3>200.044&mz_3<200.045],breaks = 100, main="")
hist(mz_raw[mz_raw>200.044&mz_raw<200.045],breaks = 100, main="")

massCheck <- function(mass){
  par(mfrow=c(4,1))
  hist(mz_1[mz_1>mass&mz_1<(mass+0.001)],breaks = 100, main="")
  hist(mz_2[mz_2>mass&mz_2<(mass+0.001)],breaks = 100, main="")
  hist(mz_3[mz_3>mass&mz_3<(mass+0.001)],breaks = 100, main="")
  hist(mz_raw[mz_raw>mass&mz_raw<(mass+0.001)],breaks = 100, main="")
  par(mfrow=c(1,1))
}
massCheck(200.044)
massCheck(117.078979+1.007267)

# Find medians ----

span <- 7
lowest_mz <- floor(min(mz_raw))-floor(min(mz_raw))%%span
highest_mz <- floor(max(mz_raw))+floor(max(mz_raw))%%span

bin_minima <- seq(lowest_mz, highest_mz, span)
getMeds <- function(bin_min, mz_raw){
  span <- round(bin_min/span)
  raw_i <- mz_raw[mz_raw>(bin_min)&mz_raw<(bin_min+span)]
  peak_mids <- list()
  for(j in seq(0, span-0.1, 0.1)){
    raw_bin_j <- raw_i[raw_i>(bin_min+j)&raw_i<(bin_min+j+0.1)]
    if(length(raw_bin_j)>1){
      hist_bin_j <- hist(raw_bin_j, breaks = 10000, plot = F)
      peak_idxs <- MALDIquant:::.localMaxima(hist_bin_j$counts, span)
      peak_mids[[(j+0.1)*10]] <- hist_bin_j$mids[peak_idxs]
    }
  }
  return(peak_mids)
}

plan(multiprocess, workers = availableCores()-1)
mids <- unlist(future_lapply(bin_minima, getMeds, mz_raw))
bin_halves <- mids_clean[-length(mids_clean)]+diff(mids_clean)/2
bin_halves <- c(0, bin_halves, 1000)

# Check on median quality ----

par(mfcol=c(3,3))
par(mar=c(2,4,0,0.1))
palered <- rgb(1,0,0,0.2)

hist(mz_raw[mz_raw>60.044&mz_raw<60.045],breaks = 100, main="")
abline(v=mids, col=palered)
abline(v=bin_halves, col="blue")
hist(mz_raw[mz_raw>60.04&mz_raw<60.05],breaks = 1000, main="")
abline(v=mids, col=palered)
abline(v=bin_halves, col="blue")
hist(mz_raw[mz_raw>60.0&mz_raw<60.1],breaks = 10000, main="")
abline(v=mids, col=palered)
abline(v=bin_halves, col="blue")

hist(mz_raw[mz_raw>200.044&mz_raw<200.045],breaks = 100, main="")
abline(v=mids, col=palered)
abline(v=bin_halves, col="blue")
hist(mz_raw[mz_raw>200.04&mz_raw<200.05],breaks = 1000, main="")
abline(v=mids, col=palered)
abline(v=bin_halves, col="blue")
hist(mz_raw[mz_raw>200.0&mz_raw<200.1],breaks = 10000, main="")
abline(v=mids, col=palered)
abline(v=bin_halves, col="blue")

hist(mz_raw[mz_raw>810.101&mz_raw<810.102],breaks = 100, main="")
abline(v=mids, col=palered)
abline(v=bin_halves, col="blue")
hist(mz_raw[mz_raw>810.10&mz_raw<810.11],breaks = 1000, main="")
abline(v=mids, col=palered)
abline(v=bin_halves, col="blue")
hist(mz_raw[mz_raw>810.1&mz_raw<810.2],breaks = 10000, main="")
abline(v=mids, col=palered)
abline(v=bin_halves, col="blue")

par(mfcol=c(2,3))
pool_1_data %>%
  filterMz(c(200.04, 200.045)) %>%
  will_plotXIC(col=NA, yl=c(200.040, 200.045))
abline(h=bin_halves, col="red")
pool_2_data %>%
  filterMz(c(200.04, 200.045)) %>%
  will_plotXIC(col=NA, yl=c(200.040, 200.045))
abline(h=bin_halves, col="red")
pool_3_data %>%
  filterMz(c(200.04, 200.045)) %>%
  will_plotXIC(col=NA, yl=c(200.040, 200.045))
abline(h=bin_halves, col="red")
par(mfcol=c(2,1))
blank_data %>%
  filterMz(c(200.04, 200.045)) %>%
  will_plotXIC(col=NA, yl=c(200.040, 200.045))
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



# Check on data quality ----
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

par(mfcol=c(2,2))
raw_data_mem %>%
  filterRt(c(600, 750)) %>%
  filterMz(c(90.054, 90.058)) %>%
  will_plotXIC(col=NA, yl=c(90.054, 90.058), mn="Alanine")
mod_data %>%
  filterRt(c(600, 750)) %>%
  filterMz(c(90.054, 90.058)) %>%
  will_plotXIC(col=NA, yl=c(90.054, 90.058), mn="Alanine")

raw_data_mem %>%
  filterRt(c(650, 800)) %>%
  filterMz(c(308.088, 308.092)) %>%
  will_plotXIC(col=NA, yl=c(308.088, 308.092), mn="Glutathione")
mod_data %>%
  filterRt(c(650, 800)) %>%
  filterMz(c(308.088, 308.092)) %>%
  will_plotXIC(col=NA, yl=c(308.088, 308.092), mn="Glutathione")

par(mfrow=c(2,1))
mod_data %>%
  filterRt(c(700, 1000)) %>%
  filterMz(c(308.07, 308.12)) %>%
  will_plotXIC(col=NA, yl=c(308.07, 308.12), mn="Glutathione")
par(mfrow=c(1,1))


# Write out data ----

