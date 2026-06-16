library(dplyr)
library(tidyr)
library(ggplot2)
library(terra)
library(tidyterra)
library(maptiles)

grid <- rast("inputs/HELCOM_WP3_grid.tif")
out_repo <- "fish"
extra_dir <- paste0("../models/outputs_extra/")
out_dir <- paste0("../models/outputs/")

df <- read.csv(paste0("outputs/", out_repo, "/", out_repo, "_details.csv"))
spec_list <- df$scientific_name

############ Copy and export files/plots ############
#---------------------------------------------------#
# Binary TIFs
od <- paste0("outputs/", out_repo, "/binary_tif")
unlink(od, recursive = TRUE)
dir.create(od, recursive = TRUE)

for(i in 1:length(spec_list)){
  spec_full <- df$scientific_name[i]
  spec_short <- df$short_name[i] # shorten name if too long (issues writing files)
  
  fp <- paste0(extra_dir, "/", df$run_name[i],"/", spec_short)
  l <- list.files(fp, full.names = TRUE)
  l <- l[grepl("_TSSbin.tif", l)]
  r <- rast(l)
  r <- r*grid
  fp <- paste0(od, "/", spec_full, ".tif")
  
  writeRaster(r, fp, 
              datatype = "INT1U", NAflag = 255,
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=1", "ZLEVEL=9"),
              overwrite = TRUE)
}

#---------------------------------------------------#
# Full probability TIFs
od <- paste0("outputs/", out_repo, "/probability_tif")
unlink(od, recursive = TRUE)
dir.create(od, recursive = TRUE)

for(i in 1:length(spec_list)){
  spec_full <- df$scientific_name[i]
  spec_short <- df$short_name[i] # shorten name if too long (issues writing files)
  
  fp <- paste0(extra_dir, "/", df$run_name[i],"/", spec_short)
  l <- list.files(fp, full.names = TRUE)
  l <- l[grepl("_EMwmeanByTSS.tif", l)]
  r <- rast(l)
  r <- r*grid
  fp <- paste0(od, "/", spec_full, ".tif")
  r[r > 1000] <- 1000 # Fix in some cases probability > 1000
  
  r <- round(r)
  writeRaster(r, fp,
              datatype = "INT2U", NAflag = 65535,
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"),
              overwrite = TRUE)

  mn <- terra::minmax(r)[1]
  mx <- terra::minmax(r)[2]
  cat(spec_full, "probability: min =", mn, " - max =", mx, "\n")
}

#---------------------------------------------------#
# Binary confidence TIFs
od <- paste0("outputs/", out_repo, "/binary_confidence_tif")
unlink(od, recursive = TRUE)
dir.create(od, recursive = TRUE)

cat_names <- c("Unsuitable (high confidence)",
               "Unsuitable (low confidence)",
               "Suitable (low confidence)",
               "Suitable (high confidence)")

for(i in 1:length(spec_list)){
  spec_full <- df$scientific_name[i]
  spec_short <- df$short_name[i] # shorten name if too long (issues writing files)
  
  fp <- paste0(extra_dir, "/", df$run_name[i],"/", spec_short)
  l <- list.files(fp, full.names = TRUE)
  l <- l[grepl("_binary_confidence.tif", l)]
  r <- rast(l)
  r <- r*grid
  
  levels(r) <- data.frame(
    value = 1:4,
    suitability = cat_names
  )
  
  fp <- paste0(od, "/", spec_full, ".tif")
  writeRaster(r, fp, 
              datatype = "INT1U", NAflag = 255,
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=1", "ZLEVEL=9"),
              overwrite = TRUE)
}


#---------------------------------------------------#
# Confidence TIFs
od <- paste0("outputs/", out_repo, "/confidence_tif")
unlink(od, recursive = TRUE)
dir.create(od, recursive = TRUE)

for(i in 1:length(spec_list)){
  spec_full <- df$scientific_name[i]
  spec_short <- df$short_name[i] # shorten name if too long (issues writing files)
  
  fp <- paste0(extra_dir, "/", df$run_name[i],"/", spec_short)
  l <- list.files(fp, full.names = TRUE)
  l <- l[grepl("_confidence.tif", l)]
  l <- l[!grepl("binary_confidence.tif", l)]
  r <- rast(l)
  r <- r*grid
  fp <- paste0(od, "/", spec_full, ".tif")
  
  r <- round(r)
  writeRaster(r, fp,
              datatype = "INT2U", NAflag = 65535,
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"),
              overwrite = TRUE)
}

#---------------------------------------------------#
# Delta TIFs
od <- paste0("outputs/", out_repo, "/delta_tif")
unlink(od, recursive = TRUE)
dir.create(od, recursive = TRUE)

for(i in 1:length(spec_list)){
  spec_full <- df$scientific_name[i]
  spec_short <- df$short_name[i] # shorten name if too long (issues writing files)
  
  fp <- paste0(extra_dir, "/", df$run_name[i],"/", spec_short)
  l <- list.files(fp, full.names = TRUE)
  l <- l[grepl("_EMwmeanByTSS.tif", l)]
  r_wmean <- rast(l)
  
  l <- list.files(fp, full.names = TRUE)
  l <- l[grepl("_TSSbin.tif", l)]
  r_bin <- rast(l)
  
  em_delta <- r_wmean
  em_delta[r_bin == 0] <- 0 # Ensemble threshold binary mask
  em_delta_prob <- em_delta
  em_delta_prob <- em_delta_prob*grid
  em_delta[em_delta == 0] <- NA
  
  em_delta <- ((em_delta - minmax(em_delta)[1,]) / (minmax(em_delta)[2,] - minmax(em_delta)[1,])) * (1000 - 1) + 1
  em_delta <- round(em_delta)
  em_delta[is.na(em_delta)] <- 0
  em_delta <- em_delta*grid
  r <- em_delta
  
  fp <- paste0(od, "/", spec_full, ".tif")

  r <- round(r)
  writeRaster(r, fp,
              datatype = "INT2U", NAflag = 65535,
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"),
              overwrite = TRUE)
  
  mn <- terra::minmax(r)[1]
  mx <- terra::minmax(r)[2]
  cat(spec_full, "delta: min =", mn, " - max =", mx, "\n")
  
}

#---------------------------------------------------#
# Strict niche TIFs
od <- paste0("outputs/", out_repo, "/binary_strict_tif")
unlink(od, recursive = TRUE)
dir.create(od, recursive = TRUE)

cat_names <- c("Likely unsuitable",
               "Possibly suitable",
               "Likely suitable",
               "Very likely suitable")

for(i in 1:length(spec_list)){
  spec_full <- df$scientific_name[i]
  spec_short <- df$short_name[i] # shorten name if too long (issues writing files)
  
  fp <- paste0(extra_dir, "/", df$run_name[i],"/", spec_short)
  l <- list.files(fp, full.names = TRUE)
  l <- l[grepl("_EMwmeanByTSS.tif", l)]
  r_prob <- rast(l)
  
  l <- list.files(fp, full.names = TRUE)
  l <- l[grepl("_TSSbin.tif", l)]
  r_bin <- rast(l)
  
  l <- list.files(fp, full.names = TRUE)
  l <- l[grepl("_confidence.tif", l)]
  l <- l[!grepl("binary_confidence.tif", l)]
  r_conf <- rast(l)
  
  # Determine mid-point of presence probabilities
  r_delta <- r_bin*r_prob
  v <- terra::values(r_delta, na.rm = TRUE)
  v <- v[v > 0]
  q <- (min(v) + max(v))/2
  
  delta_thresh <- ifel(r_delta > q, 1, 0)
  conf_thresh <- ifel(r_conf > 50, 1, 0)
  
  # score: 0 = both low, 1 = one low, 2 = both high
  score <- conf_thresh + delta_thresh
  
  # 1 = Likely unsuitable (absent in SDM)
  # 2 = Possibly suitable (present, both low)
  # 3 = Likely suitable (present, one low)
  # 4 = Very likely suitable (present, both high)
  r_niche <- ifel(r_bin == 0, 1, score + 2)
  r_niche <- r_niche*grid
  
  levels(r_niche) <- data.frame(
    value = 1:4,
    suitability = cat_names
  )
  
  fp <- paste0(od, "/", spec_full, ".tif")
  r <- r_niche
  
  writeRaster(r, fp, 
              datatype = "INT1U", NAflag = 255,
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=1", "ZLEVEL=9"),
              overwrite = TRUE)
  
}

#---------------------------------------------------#
# Binary envelope TIFs
od <- paste0("outputs/", out_repo, "/binary_envelope_tif")
unlink(od, recursive = TRUE)
dir.create(od, recursive = TRUE)

df_mask <- read.csv("inputs/WP3_species_mask.csv")
hard <- rast("inputs/hard_binary_10perc_250m.tif")

for(i in 1:length(spec_list)){
  spec_full <- df$scientific_name[i]
  spec_short <- df$short_name[i] # shorten name if too long (issues writing files)
  
  fp <- paste0(extra_dir, "/", df$run_name[i],"/", spec_short)
  l <- list.files(fp, full.names = TRUE)
  l <- l[grepl("_TSSbin.tif", l)]
  r <- rast(l)
  
  # Mask by hard substrate?
  if(df_mask$hard_substrate_mask[df_mask$scientific_name == spec_full] == 1){r <- r*hard}
  
  r <- r*grid
  fp <- paste0(od, "/", spec_full, ".tif")
  
  writeRaster(r, fp, 
              datatype = "INT1U", NAflag = 255,
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=1", "ZLEVEL=9"),
              overwrite = TRUE)
}
