# Packages
library(terra)
library(dplyr)
library(tidyr)
library(FNN)
library(corrplot)
library(usdm)
library(ggplot2)
predictor_dir <- "predictors_all/"

############# Invertebrates #############
pred_select <- c(
  "aspect", "bpi_4km_20km", "bpi_500m_5km", "chloro_a_mean",
  "depth", "depth_attenuated_exposure", "substrate_hard_binary",
  "substrate_mixed_binary", "nitrate_surface_mean", "oxygen_bottom_mean",
  "salinity_bottom_mean", "salinity_surface_mean", "sea_ice_fraction_mean", 
  "secchi_depth_mean", "slope", "temperature_bottom_mean", "rugosity", "wave_significant_height_mean" 
)

l <- list.files("../observations/invertebrates", full.names = TRUE)
ls <- list.files("../observations/invertebrates")
n <- tools::file_path_sans_ext(ls)

for(i in 1:length(l)){
  df <- readRDS(l[i])
  df <- df %>% select(-c(quantity))
  
  if(i == 1){obs <- df}
  if(i > 1){obs <- rbind(obs, df)}
  #obs <- unique(obs) # Only take unique combinations (possibly same grid cells)?
  
  cat("Iteration", i, "complete - ")
  
}

obs <- obs[complete.cases(obs), ]
obs <- obs[, pred_select]

# Correlation plot
m <- cor(obs)
corrplot(m, method = "number", order = "AOE", diag = FALSE)

# VIF test
vifs <- vif(obs)
#v <- vifstep(obs, th = 10, size = nrow(obs))
v <- vifstep(obs, th = 10, size = nrow(obs), keep = c("chloro_a_mean", "salinity_bottom_mean"))
v

pred_select <- v
pred_select <- pred_select@results$Variables

new_predictor_dir <- "predictors_invertebrates"
file.remove(list.files(new_predictor_dir, full.names = TRUE))
dir.create(new_predictor_dir)

l <- list.files(predictor_dir, pattern = ".tif$", full.names = TRUE)
ls <- list.files(predictor_dir, pattern = ".tif$", full.names = FALSE)
ln <- tools::file_path_sans_ext(ls)

for(i in 1:length(pred_select)){
  p <- pred_select[i]
  f <- l[which(ln %in% p)]
  file.copy(f, paste0(new_predictor_dir, "/", p, ".tif"))
}

############# Macrophytes #############
pred_select <- c(
  "bpi_4km_20km", "bpi_500m_5km", "depth", "substrate_soft", "substrate_sand",
  "substrate_hard", "substrate_coarse", "slope", "depth_attenuated_exposure",
  "secchi_depth_mean", "nitrate_surface_mean", "ph_surface_mean", "salinity_surface_mean", 
  "temperature_surface_mean"
)

l <- list.files("../observations/macrophytes", full.names = TRUE)
ls <- list.files("../observations/macrophytes")
n <- tools::file_path_sans_ext(ls)

for(i in 1:length(l)){
  df <- readRDS(l[i])
  df <- df %>% select(-c(quantity, country))
  
  if(i == 1){obs <- df}
  if(i > 1){obs <- rbind(obs, df)}
  #obs <- unique(obs) # Only take unique combinations (possibly same grid cells)?
  
  cat("Iteration", i, "complete - ")
  
}

obs <- obs[complete.cases(obs), ]
obs <- obs[, pred_select]

# Correlation plot
m <- cor(obs)
corrplot(m, method = "number", order = "AOE", diag = FALSE)

# VIF test
vifs <- vif(obs)
v <- vifstep(obs, th = 10, size = nrow(obs))
v

pred_select <- v
pred_select <- pred_select@results$Variables

new_predictor_dir <- "predictors_macrophytes"
file.remove(list.files(new_predictor_dir, full.names = TRUE))
dir.create(new_predictor_dir)

l <- list.files(predictor_dir, pattern = ".tif$", full.names = TRUE)
ls <- list.files(predictor_dir, pattern = ".tif$", full.names = FALSE)
ln <- tools::file_path_sans_ext(ls)

for(i in 1:length(pred_select)){
  p <- pred_select[i]
  f <- l[which(ln %in% p)]
  file.copy(f, paste0(new_predictor_dir, "/", p, ".tif"))
}

############# Fish #############
pred_select <- c(
  "bpi_500m_5km", "depth", "substrate_soft", "substrate_sand",
  "substrate_hard", "slope", "depth_attenuated_exposure",
  "secchi_depth_mean", "oxygen_bottom_mean", "primary_prod_surface_mean",
  "salinity_bottom_mean", "salinity_surface_mean", "sea_ice_fraction_mean",
  "temperature_bottom_mean", "wave_significant_height_mean"
)

l <- list.files("../observations/fish", full.names = TRUE)
ls <- list.files("../observations/fish")
n <- tools::file_path_sans_ext(ls)

for(i in 1:length(l)){
  df <- readRDS(l[i])
  df <- df %>% select(-c(quantity, method))
  
  if(i == 1){obs <- df}
  if(i > 1){obs <- rbind(obs, df)}
  #obs <- unique(obs) # Only take unique combinations (possibly same grid cells)?
  
  cat("Iteration", i, "complete - ")
}

obs <- obs[complete.cases(obs), ]
obs <- obs[, pred_select]
obs <- as.data.frame(obs)

# Correlation plot
m <- cor(obs)
corrplot(m, method = "number", order = "AOE", diag = FALSE)

# VIF test
vifs <- vif(obs)
v <- vifstep(obs, th = 10, size = nrow(obs))
v

pred_select <- v
pred_select <- pred_select@results$Variables

new_predictor_dir <- "predictors_fish"
file.remove(list.files(new_predictor_dir, full.names = TRUE))
dir.create(new_predictor_dir)

l <- list.files(predictor_dir, pattern = ".tif$", full.names = TRUE)
ls <- list.files(predictor_dir, pattern = ".tif$", full.names = FALSE)
ln <- tools::file_path_sans_ext(ls)

for(i in 1:length(pred_select)){
  p <- pred_select[i]
  f <- l[which(ln %in% p)]
  file.copy(f, paste0(new_predictor_dir, "/", p, ".tif"))
}

###### Final predictor list
l <- list.files("predictors_fish")
l <- tools::file_path_sans_ext(l)
df <- data.frame("species_group" = "fish", "predictor" = l)
dff <- df

l <- list.files("predictors_macrophytes")
l <- tools::file_path_sans_ext(l)
df <- data.frame("species_group" = "macrophytes", "predictor" = l)
dff <- rbind(dff, df)

l <- list.files("predictors_invertebrates")
l <- tools::file_path_sans_ext(l)
df <- data.frame("species_group" = "invertebrates", "predictor" = l)
dff <- rbind(dff, df)

write.csv(dff, "predictor_list.csv", row.names = FALSE)
