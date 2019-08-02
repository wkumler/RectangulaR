# Converting relevant MS data from .raw to .mzML (in profile mode)
# Note that these files are written to and from the Z: drive
# Should only need to be run once

# Setup things ----

library(dplyr)
library(future.apply)
library(beepr)



# Identify files to be converted ----

base_path <- "Z:/1_QEdata/LTC/DATA/HILIC/190718_DepthProfiles_FK180310/"
if(!dir.exists(base_path)){stop("Is the network drive mapped?")}

blank_1_raw_path <- paste0(base_path, "190715_Blk_KM1906U14-Blk_C.raw")
stds_h2o_1_raw_path <- paste0(base_path, "190715_Std_4uMStdsMix1InH2O_1.raw")
stds_h2o_2_raw_path <- paste0(base_path, "190715_Std_4uMStdsMix2InH2O_1.raw")
stds_mat_1_raw_path <- paste0(base_path, "190715_Std_4uMStdsMix1InMatrix_1.raw")
stds_mat_2_raw_path <- paste0(base_path, "190715_Std_4uMStdsMix2InMatrix_1.raw")
full_poo_1_raw_path <- paste0(base_path, "190715_Poo_TruePooFK180310_Full1.raw")
full_poo_2_raw_path <- paste0(base_path, "190715_Poo_TruePooFK180310_Full2.raw")
full_poo_3_raw_path <- paste0(base_path, "190715_Poo_TruePooFK180310_Full3.raw")

filepaths <- sapply(ls(pattern = "*._raw_path$"), get, USE.NAMES = F)
if(any(!sapply(filepaths, file.exists))){stop("Unable to find all files")}



# Generate system conversion commands ----
msconvertMaker <- function(filepath){
  pos <- "msconvert" %>%
    paste("--mzML") %>%
    paste('--filter "polarity positive"') %>%
    paste('-z', filepath) %>%
    paste("-o Z:/1_QEdata/Will/RectangulaRdata/positive")
  neg <- "msconvert" %>%
    paste("--mzML") %>%
    paste('--filter "polarity negative"') %>%
    paste('-z', filepath) %>%
    paste("-o Z:/1_QEdata/Will/RectangulaRdata/negative/")
  return(c(pos, neg))
}
mscommands <- msconvertMaker(filepaths)



# And convert them ----
plan(multiprocess, workers = availableCores()-1)
converted <- future_lapply(mscommands, system)



# Then check whether all files were converted ----
files_tbd <- c(
  paste0("Z:/1_QEdata/Will/RectangulaRdata/positive/",
                   sapply(basename(filepaths), gsub, pattern="raw", 
                          replacement="mzML", USE.NAMES = F)),
  paste0("Z:/1_QEdata/Will/RectangulaRdata/negative/",
         sapply(basename(filepaths), gsub, pattern="raw", 
                replacement="mzML", USE.NAMES = F))
)
files_det <- c(
  list.files("Z:/1_QEdata/Will/RectangulaRdata/positive/", full.names = T),
  list.files("Z:/1_QEdata/Will/RectangulaRdata/negative/", full.names = T)
)

if(any(unlist(converted))!=0){warning("Unable to convert all files???")}

if(length(files_det)!=length(files_tbd)){
  print(paste("Failed files:", setdiff(files_tbd, files_det)))
}
beep(2)
