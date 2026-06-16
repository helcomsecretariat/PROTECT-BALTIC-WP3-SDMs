library(dplyr)
library(tidyr)
library(ggplot2)
library(terra)
library(tidyterra)
library(maptiles)
#---------------------------------------------------#
# Binary strict niche PNGs with land

grid <- rast("inputs/HELCOM_WP3_grid.tif")
land <- rast("inputs/WP3_land_mask.tif")
sea <- (land*-1)+1
sea[grid == 1] <- 0

spec_group_all <- c("fish", "macrophytes", "invertebrates")

for(j in 1:length(spec_group_all)){
  
  spec_group <- spec_group_all[j]
  
  source_dir <- paste0("outputs/", spec_group, "/binary_strict_tif")
  out_dir <- paste0("outputs/", spec_group, "/_pngs/binary_strict_png")
  
  l <- list.files(source_dir, full.names = TRUE, pattern = ".tif")
  l <- l[!grepl(".tif.aux.xml", l)]
  ls <- basename(l)
  n <- tools::file_path_sans_ext(ls)
  
  #unlink(out_dir, recursive = TRUE)
  dir.create(out_dir, recursive = TRUE)
  
  cat_names <- c("Likely unsuitable",
                 "Possibly suitable",
                 "Likely suitable",
                 "Very likely suitable")
  
  for(i in 1:length(n)){
    fp <- paste0(source_dir, "/", ls[i])
    dest_fp <- paste0(out_dir, "/", n[i], ".png")
    
    r <- rast(fp)
    r[land == 1] <- 0
    r[sea == 1] <- -1
    
    # Export plot
    levels(r) <- NULL
    dfr <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
    labs <- c("Sea", "Land", cat_names)
    cols <- c("grey90", 'grey80', 'white','#ffffbf','#feb24c', '#d7301f')
    names(cols) <- labs
    names(dfr)[3] <- "class"
    
    dfr$category <- factor(dfr$class, levels = -1:4, labels = labs)
    
    temp_plot <- ggplot(dfr, aes(x, y, fill = category)) +
      geom_raster(key_glyph = "polygon") +
      coord_equal() +
      scale_fill_manual(
        values = cols,
        name = "",
        breaks = labs[3:6]  # exclude "Land" from the legend
      ) +
      guides(fill = guide_legend(override.aes = list(color = "black", linewidth = 0.4))) +
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