# Normalize mz swaths to the median value
# Read in MS data, find medians, normalize to medians, and write out MS data

# Setup things ----

library(dplyr)
library(future.apply)
library(MSnbase)
library(beepr)

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
std_mix1_h2o_data <- readMSData(files = paste0(pos_data_path, std_mix1_h2o_name), mode = "inMemory", 
                                centroided. = F, msLevel. = 1)
std_mix1_mat_data <- readMSData(files = paste0(pos_data_path, std_mix1_mat_name), mode = "inMemory", 
                                centroided. = F, msLevel. = 1)
std_mix2_h2o_data <- readMSData(files = paste0(pos_data_path, std_mix2_h2o_name), mode = "inMemory", 
                                centroided. = F, msLevel. = 1)
std_mix2_mat_data <- readMSData(files = paste0(pos_data_path, std_mix2_mat_name), mode = "inMemory", 
                                centroided. = F, msLevel. = 1)

mz_1 <- unlist(unname(mz(pool_1_data)))
mz_2 <- unlist(unname(mz(pool_2_data)))
mz_3 <- unlist(unname(mz(pool_3_data)))
mz_raw <- c(mz_1, mz_2, mz_3)



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
bin_halves <- mids[-length(mids)]+diff(mids)/2
bin_halves <- c(0, bin_halves, 1000)



# Apply flattening and write out data ----
middling <- function(scan){
  scan_mz <- scan@mz
  scan_cut <- cut(scan_mz, bin_halves)
  new_mz <- mids[scan_cut]
  scan@mz <- new_mz
  return(scan)
}
for(i in ls(pattern = "data$")){
  centered_data <- future_eapply(i@assayData, FUN = middling)
  i@assayData <- as.environment(centered_data)
  writeMSData(i, file = paste0("Z:/1_QEdata/Will/RectangulaRdata/pos_norm/", i, ".mzML"))
}


