library(dplyr)
library(tidyr)
library(ggplot2)
library(terra)
library(tidyterra)
library(maptiles)
#---------------------------------------------------#
# Probability PDFs with land

spec_group_all <- c("fish", "macrophytes", "invertebrates")

for(j in 1:length(spec_group_all)){
  
  spec_group <- spec_group_all[j]
  
  grid <- rast("inputs/HELCOM_WP3_grid.tif")
  land <- rast("inputs/WP3_land_mask.tif")
  sea <- (land-1)*-1
  sea[land == 1] <- NA
  sea[grid == 1] <- NA
  
  # Build a land-only data frame once per species group (land doesn't change)
  land_df <- as.data.frame(land, xy = TRUE, na.rm = TRUE)
  names(land_df)[3] <- "land"
  land_df <- land_df[land_df$land == 1, ]
  
  sea_df <- as.data.frame(sea, xy = TRUE, na.rm = TRUE)
  names(sea_df)[3] <- "sea"
  sea_df <- sea_df[sea_df$sea == 1, ]
  
  source_dir <- paste0("outputs/", spec_group, "/probability_tif")      # adjust as needed
  out_dir    <- paste0("outputs/", spec_group, "/_pngs/probability_png")
  
  l  <- list.files(source_dir, full.names = TRUE, pattern = ".tif")
  l  <- l[!grepl(".tif.aux.xml", l)]
  ls <- basename(l)
  n  <- tools::file_path_sans_ext(ls)
  
  #unlink(out_dir, recursive = TRUE)
  dir.create(out_dir, recursive = TRUE)
  
  for(i in 1:length(n)){
    
    fp <- paste0(source_dir, "/", ls[i])
    dest_fp <- paste0(out_dir, "/", n[i], ".png")
    
    r <- rast(fp)
    #r <- r/10
    # Mask out land so it doesn't contaminate the colour ramp
    r[land == 1] <- NA
    
    dfr <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
    names(dfr)[3] <- "value"
    
    temp_plot <- ggplot() +
      # Land underneath
      geom_raster(data = land_df, aes(x = x, y = y), fill = "grey80") +
      geom_raster(data = sea_df, aes(x = x, y = y), fill = "grey90") +
      # Continuous data on top
      geom_raster(data = dfr, aes(x = x, y = y, fill = value)) +
      coord_equal() +
      scale_fill_distiller(name = "Probability",
                           na.value = "transparent",
                           palette = "Spectral",
                           direction = -1) +
      theme(panel.background = element_rect(fill = "white"), 
            plot.background = element_rect(fill = "white", color = NA),
            legend.background = element_rect(fill = "white", color = NA),
            legend.key = element_rect(fill = "white", color = NA),
            panel.grid = element_blank(),
            axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank())
    
    png(filename = dest_fp,
        width = 9.5,
        height = 8.5,
        units = "in",
        res = 600)
    print(temp_plot)
    dev.off()
  }
}