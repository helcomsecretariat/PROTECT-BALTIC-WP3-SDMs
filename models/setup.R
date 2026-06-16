########## Setup ##########
library(dplyr)
library(tidyr)
library(terra)

########## Invertebrates ##########
l <- list.files("../observations/invertebrates", full.names = TRUE)
ls <- list.files("../observations/invertebrates")
n <- tools::file_path_sans_ext(ls)

for(i in 1:length(l)){
  df <- readRDS(l[i])
  df$scientific_name <- n[i]
  
  df <- df %>% select(scientific_name, quantity)
  
  df <- df %>%
    group_by(scientific_name) %>%
    summarise(quantity = sum(quantity)) %>%
    ungroup()
  
  if(i == 1){spec_list <- df}
  if(i > 1){spec_list <- rbind(spec_list, df)}
  
}

# Order species list by quantity
spec_list <- spec_list %>% arrange(desc(quantity))
spec_list$run_initiated <- 0
spec_list$run_completed <- 0
spec_list$run_time <- 0
spec_list$models <- NA
spec_list$num_presence <- 0
spec_list$num_absence <- 0
spec_list$run_name <- NA
spec_list$observation_methods <- NA
spec_list$method_predict_to <- NA
spec_list$mean_TSS_ensemble <- NA
spec_list$mean_TSS_all_submodels <- NA
spec_list$metric_select_and_weight <- NA
spec_list$metric_binary <- NA
write.csv(spec_list, "spec_list_invertebrates.csv", row.names = FALSE)

########## Fish ##########
l <- list.files("../observations/fish", full.names = TRUE)
ls <- list.files("../observations/fish")
n <- tools::file_path_sans_ext(ls)

for(i in 1:length(l)){
  df <- readRDS(l[i])
  df$scientific_name <- n[i]
  
  df <- df %>% select(scientific_name, quantity)
  
  df <- df %>%
    group_by(scientific_name) %>%
    summarise(quantity = sum(quantity)) %>%
    ungroup()
  
  if(i == 1){spec_list <- df}
  if(i > 1){spec_list <- rbind(spec_list, df)}
  
}

# Specify use of method co-variate in models
spec_list$use_method_cov <- TRUE # Default TRUE as most species performed best with co-variate
spec_change <- c("Coregonus.maraena", "Gobius.niger", "Liparis.liparis",
                 "Lycodes.vahlii", "Spinachia.spinachia")

spec_list$use_method_cov[spec_list$scientific_name %in% spec_change] <- FALSE


# Order species list by quantity
spec_list <- spec_list %>% arrange(desc(quantity))
spec_list$run_initiated <- 0
spec_list$run_completed <- 0
spec_list$run_time <- 0
spec_list$models <- NA
spec_list$num_presence <- 0
spec_list$num_absence <- 0
spec_list$run_name <- NA
spec_list$observation_methods <- NA
spec_list$method_predict_to <- NA
spec_list$mean_TSS_ensemble <- NA
spec_list$mean_TSS_all_submodels <- NA
spec_list$metric_select_and_weight <- NA
spec_list$metric_binary <- NA
write.csv(spec_list, "spec_list_fish.csv", row.names = FALSE)

########## Macrophytes ##########
l <- list.files("../observations/macrophytes", full.names = TRUE)
ls <- list.files("../observations/macrophytes")
n <- tools::file_path_sans_ext(ls)

for(i in 1:length(l)){
  df <- readRDS(l[i])
  df$scientific_name <- n[i]
  
  df <- df %>% select(scientific_name, quantity)
  
  df <- df %>%
    group_by(scientific_name) %>%
    summarise(quantity = sum(quantity)) %>%
    ungroup()
  
  if(i == 1){spec_list <- df}
  if(i > 1){spec_list <- rbind(spec_list, df)}
  
}

# Specify use of method co-variate in models
spec_list$use_method_cov <- FALSE # Default FALSE
spec_change <- c("Acrosiphonia.arcta", "Ceramium.virgatum", "Chaetomorpha.linum",
                 "Elodea.nuttallii", "Eudesme.virescens", "Fucus.vesiculosus.radicans",
                 "Halosiphon.tomentosus", "Ranunculus.peltatus.peltatus", "Ruppia.maritima",
                 "Spermothamnion.repens", "Spongomorpha.aeruginosa", "Stictyosiphon.tortilis",
                 "Tolypella.nidifica", "Zannichellia.major", "Zannichellia.palustris",
                 "Agarophyton.vermiculophyllum", "Cladophora.glomerata", "Spongonema.tomentosum")

spec_list$use_method_cov[spec_list$scientific_name %in% spec_change] <- TRUE

# Order species list by quantity
spec_list <- spec_list %>% arrange(desc(quantity))
spec_list$run_initiated <- 0
spec_list$run_completed <- 0
spec_list$run_time <- 0
spec_list$models <- NA
spec_list$num_presence <- 0
spec_list$num_absence <- 0
spec_list$run_name <- NA
spec_list$observation_methods <- NA
spec_list$method_predict_to <- NA
spec_list$mean_TSS_ensemble <- NA
spec_list$mean_TSS_all_submodels <- NA
spec_list$metric_select_and_weight <- NA
spec_list$metric_binary <- NA
write.csv(spec_list, "spec_list_macrophytes.csv", row.names = FALSE)
