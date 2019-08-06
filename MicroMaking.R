# After raw conversion, make small versions of each file to play with

library(MSnbase)
library(dplyr)

files_det <- list.files("Z:/1_QEdata/Will/RectangulaRdata/positive/", full.names = T)

for(i in files_det){
  readMSData(i, mode = "inMemory", centroided. = F, msLevel. = 1) %>%
    filterRt(c(0, 300)) %>%
    filterMz(c(0, 120)) %>%
    writeMSData(file = gsub("positive", "microdata", i))
}

files_fin <- sapply(files_det, gsub, pattern="positive", replacement="microdata", USE.NAMES = F)
files_new <- list.files("Z:/1_QEdata/Will/RectangulaRdata/microdata/", full.names = T)
if(!identical(files_new, files_fin)){stop("Not all files converted")}

lapply(files_fin, readMSData, mode = "onDisk", centroided = F, msLevel. = 1)
