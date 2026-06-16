library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(terra)
library(worrms)
library(dplyr)
library(purrr)
library(tibble)
library(rgbif)
#source("fetch_species_ids.R")

#grid <- rast("inputs/HELCOM_WP3_grid.tif")

############ Filter species according to evaluation ############
#--------------------------------------------------------------#
############ Invertebrates ############
mod_id <- "run3_invertebrates" # Unique ID for model runs
out_repo <- "invertebrates"
extra_dir <- paste0("../models/outputs_extra/", mod_id)
out_dir <- paste0("../models/outputs/", mod_id)
dir.create(paste0("outputs/", out_repo), recursive = TRUE)

spec_list <- list.files(extra_dir)
spec_list <- spec_list[spec_list != "_run_finished" & spec_list != "_run_started"]

run_eval <- read.csv(paste0(extra_dir, "/_run_finished/_", mod_id, "_details.csv"))
#run_eval <- filter(run_eval, run_initiated == 1)
run_eval$scientific_name <- gsub(" ", ".", run_eval$scientific_name)
run_eval$short_name <- substr(run_eval$scientific_name, 1, 25)
run_eval <- run_eval[order(run_eval$scientific_name), ]
rownames(run_eval) <- NULL

# For invertebrates, specify which species should use the fish model
run_eval <- filter(run_eval, scientific_name != "Carcinus.maenas")
run_eval$output_repo <- out_repo

# Add flags
run_eval <- run_eval %>% 
  add_column(flag = "green", .after = "scientific_name") %>%
  add_column(comment = NA, .after = "flag") %>%
  add_column(use_method_cov = FALSE, .after = "quantity") 

run_eval$flag[run_eval$num_presence < 150] <- "yellow"
run_eval$comment[run_eval$num_presence < 150] <- "Low number of positive cases (presences) in model training data."

# Fetch evaluations
run_eval$kept_algos <- NA
run_eval$n_kept_algos <- NA

for(i in 1:nrow(run_eval)){
  fp <- paste0(extra_dir, "/", run_eval$short_name[i], "/", run_eval$short_name[i], "_validation_scores.csv")
  if(!file.exists(fp)){next}
  df <- read.csv(fp)
  
  # Check which models were kept in the final ensemble
  df <- df[df$metric.eval == "TSS",]
  df <- df[df$kept_in_ensemble == TRUE,] # only models kept in the ensemble?
  
  k <- unique(df$algo)
  run_eval$kept_algos[i] <- paste0(k, collapse = ",")
  run_eval$n_kept_algos[i] <- length(k)
}


# Remove extremely poor models entirely
run_keep <- run_eval %>%
  filter(mean_TSS_ensemble >= 0.5) %>%
  filter(mean_TSS_all_submodels >= 0.4) %>%
  filter(n_kept_algos >= 2) %>%
  filter(num_presence >= 0)

run_filter <- run_eval[which(!run_eval$scientific_name %in% run_keep$scientific_name),]

if(nrow(run_filter) > 0){
  run_filter <- run_eval[which(!run_eval$scientific_name %in% run_keep$scientific_name),]
  run_filter$flag <- "red"
  run_filter$comment <- "Species model excluded because it did not pass the evaluation criteria."
}

run_keep <- run_keep[order(run_keep$scientific_name),]
#spec_list <- run_keep$scientific_name

############ Fetch common names and AphiaIDs ############
# Note, in some cases, AphiaID will have to be added manually,
# e.g, if the map/model represents two species
#run_keep <- fetch_species_ids(run_keep)

readr::write_excel_csv(run_filter, paste0("outputs/", out_repo, "/", out_repo, "_exclude.csv"))
readr::write_excel_csv(run_keep, paste0("outputs/", out_repo, "/",out_repo, "_details.csv"))

#--------------------------------------------------------------#
############ Fish ############
mod_id <- "run3_fish" # Unique ID for model runs
out_repo <- "fish"
extra_dir <- paste0("../models/outputs_extra/", mod_id)
out_dir <- paste0("../models/outputs/", mod_id)
dir.create(paste0("outputs/", out_repo), recursive = TRUE)

spec_list <- list.files(extra_dir)
spec_list <- spec_list[spec_list != "_run_finished" & spec_list != "_run_started"]

run_eval <- read.csv(paste0(extra_dir, "/_run_finished/_", mod_id, "_details.csv"))
#run_eval <- filter(run_eval, run_initiated == 1)
run_eval$scientific_name <- gsub(" ", ".", run_eval$scientific_name)
run_eval$short_name <- substr(run_eval$scientific_name, 1, 25)
run_eval <- run_eval[order(run_eval$scientific_name), ]
rownames(run_eval) <- NULL

# Specify which species are invertebrates
run_eval$output_repo <- out_repo
run_eval$output_repo[run_eval$scientific_name == "Cancer.pagurus"] <- "invertebrates"
run_eval$output_repo[run_eval$scientific_name == "Homarus.gammarus"] <- "invertebrates"
run_eval$output_repo[run_eval$scientific_name == "Hyas.araneus"] <- "invertebrates"
run_eval$output_repo[run_eval$scientific_name == "Loligo.forbesii"] <- "invertebrates"
run_eval$output_repo[run_eval$scientific_name == "Macropodia.rostrata"] <- "invertebrates"
run_eval$output_repo[run_eval$scientific_name == "Nephrops.norvegicus"] <- "invertebrates"
run_eval$output_repo[run_eval$scientific_name == "Pagurus.bernhardus"] <- "invertebrates"
run_eval$output_repo[run_eval$scientific_name == "Crangon.crangon"] <- "invertebrates"
run_eval$output_repo[run_eval$scientific_name == "Alloteuthis.subulata"] <- "invertebrates"
run_eval$output_repo[run_eval$scientific_name == "Carcinus.maenas"] <- "invertebrates"

# For fish, specify which species should use the invertebrate model
run_eval <- filter(run_eval, scientific_name != "Pagurus.bernhardus")
run_eval <- filter(run_eval, scientific_name != "Crangon.crangon")

# Add flags
run_eval <- run_eval %>% 
  add_column(flag = "green", .after = "scientific_name") %>%
  add_column(comment = NA, .after = "flag")

run_eval$flag[run_eval$num_presence < 150] <- "yellow"
run_eval$comment[run_eval$num_presence < 150] <- "Low number of positive cases (presences) in model training data."

# Fetch evaluations
run_eval$kept_algos <- NA
run_eval$n_kept_algos <- NA

for(i in 1:nrow(run_eval)){
  fp <- paste0(extra_dir, "/", run_eval$short_name[i], "/", run_eval$short_name[i], "_validation_scores.csv")
  if(!file.exists(fp)){next}
  df <- read.csv(fp)
  
  # Check which models were kept in the final ensemble
  df <- df[df$metric.eval == "TSS",]
  df <- df[df$kept_in_ensemble == TRUE,] # only models kept in the ensemble?
  
  k <- unique(df$algo)
  run_eval$kept_algos[i] <- paste0(k, collapse = ",")
  run_eval$n_kept_algos[i] <- length(k)
}

# Remove extremely poor models entirely
run_keep <- run_eval %>%
  filter(mean_TSS_ensemble >= 0.5) %>%
  filter(mean_TSS_all_submodels >= 0.4) %>%
  filter(n_kept_algos >= 2) %>%
  filter(num_presence >= 0)

run_filter <- run_eval[which(!run_eval$scientific_name %in% run_keep$scientific_name),]

if(nrow(run_filter) > 0){
  run_filter <- run_eval[which(!run_eval$scientific_name %in% run_keep$scientific_name),]
  run_filter$flag <- "red"
  run_filter$comment <- "Species model excluded because it did not pass the evaluation criteria."
}

run_keep <- run_keep[order(run_keep$scientific_name),]
#spec_list <- run_keep$scientific_name

############ Fetch common names and AphiaIDs ############
# Note, in some cases, AphiaID will have to be added manually,
# e.g, if the map/model represents two species
#run_keep <- fetch_species_ids(run_keep)

readr::write_excel_csv(run_filter, paste0("outputs/", out_repo, "/", out_repo, "_exclude.csv"))
readr::write_excel_csv(run_keep, paste0("outputs/", out_repo, "/",out_repo, "_details.csv"))

#--------------------------------------------------------------#
############ Macrophytes ############
mod_id <- "run3_macrophytes" # Unique ID for model runs
out_repo <- "macrophytes"
extra_dir <- paste0("../models/outputs_extra/", mod_id)
out_dir <- paste0("../models/outputs/", mod_id)
dir.create(paste0("outputs/", out_repo), recursive = TRUE)

spec_list <- list.files(extra_dir)
spec_list <- spec_list[spec_list != "_run_finished" & spec_list != "_run_started"]

run_eval <- read.csv(paste0(extra_dir, "/_run_finished/_", mod_id, "_details.csv"))
#run_eval <- filter(run_eval, run_initiated == 1)
run_eval$scientific_name <- gsub(" ", ".", run_eval$scientific_name)
run_eval$short_name <- substr(run_eval$scientific_name, 1, 25)
run_eval <- run_eval[order(run_eval$scientific_name), ]
rownames(run_eval) <- NULL

run_eval$output_repo <- out_repo

# Add flags
run_eval <- run_eval %>% 
  add_column(flag = "green", .after = "scientific_name") %>%
  add_column(comment = NA, .after = "flag")

run_eval$flag[run_eval$num_presence < 150] <- "yellow"
run_eval$comment[run_eval$num_presence < 150] <- "Low number of positive cases (presences) in model training data."

# Fetch evaluations
run_eval$kept_algos <- NA
run_eval$n_kept_algos <- NA

for(i in 1:nrow(run_eval)){
  fp <- paste0(extra_dir, "/", run_eval$short_name[i], "/", run_eval$short_name[i], "_validation_scores.csv")
  if(!file.exists(fp)){next}
  df <- read.csv(fp)
  
  # Check which models were kept in the final ensemble
  df <- df[df$metric.eval == "TSS",]
  df <- df[df$kept_in_ensemble == TRUE,] # only models kept in the ensemble?
  
  k <- unique(df$algo)
  run_eval$kept_algos[i] <- paste0(k, collapse = ",")
  run_eval$n_kept_algos[i] <- length(k)
}


# Remove extremely poor models entirely
run_keep <- run_eval %>%
  filter(mean_TSS_ensemble >= 0.5) %>%
  filter(mean_TSS_all_submodels >= 0.4) %>%
  filter(n_kept_algos >= 2) %>%
  filter(num_presence >= 0)

run_filter <- run_eval[which(!run_eval$scientific_name %in% run_keep$scientific_name),]

if(nrow(run_filter) > 0){
  run_filter <- run_eval[which(!run_eval$scientific_name %in% run_keep$scientific_name),]
  run_filter$flag <- "red"
  run_filter$comment <- "Species model excluded because it did not pass the evaluation criteria."
}

run_keep <- run_keep[order(run_keep$scientific_name),]
#spec_list <- run_keep$scientific_name

############ Fetch common names and AphiaIDs ############
# Note, in some cases, AphiaID will have to be added manually,
# e.g, if the map/model represents two species
#run_keep <- fetch_species_ids(run_keep)

readr::write_excel_csv(run_filter, paste0("outputs/", out_repo, "/", out_repo, "_exclude.csv"))
readr::write_excel_csv(run_keep, paste0("outputs/", out_repo, "/",out_repo, "_details.csv"))

#--------------------------------------------------------------#
############ Fix invertebrate and fish lists ############
df_invert <- read.csv("outputs/invertebrates/invertebrates_details.csv")
df_fish <- read.csv("outputs/fish/fish_details.csv")

spec1 <- select(df_invert, scientific_name, run_name, output_repo)
spec2 <- select(df_fish, scientific_name, run_name,output_repo)

dups <- left_join(spec1, spec2, by = "scientific_name")

df_invert <- rbind(df_invert, df_fish)
df_invert <- filter(df_invert, output_repo == "invertebrates")
df_fish <- filter(df_fish, output_repo == "fish")

df_invert <- df_invert[order(df_invert$scientific_name), ]
rownames(df_invert) <- NULL

df_fish <- df_fish[order(df_fish$scientific_name), ]
rownames(df_fish) <- NULL

#write.csv(df_invert, "outputs/invertebrates/invertebrates_details.csv", row.names = FALSE, fileEncoding = "UTF-8-BOM")
#write.csv(df_fish, "outputs/fish/fish_details.csv", row.names = FALSE, fileEncoding = "UTF-8-BOM")

readr::write_excel_csv(df_invert, "outputs/invertebrates/invertebrates_details.csv")
readr::write_excel_csv(df_fish, "outputs/fish/fish_details.csv")

#--------------------------------------------------------------#
############ Add internal flags ############
# Fish
df <- read.csv("outputs/fish/fish_details.csv")
df_int <- read.csv("inputs/fish_internal_flag.csv")

df1 <- df %>% select(scientific_name, flag, comment)
df2 <- df_int

dfj <- left_join(df1, df2, by = "scientific_name")

dfj <- dfj |> unite("comment", c("internal_comment", "comment"), 
                    remove = FALSE, na.rm = TRUE, sep = " ")

flag_levels <- c("red", "yellow", "green")

f1 <- factor(dfj$flag, levels = flag_levels, ordered = TRUE)
f2 <- factor(dfj$internal_flag, levels = flag_levels, ordered = TRUE)

dfj$combined_flag <- as.character(pmin(f1, f2))

df$flag <- dfj$combined_flag
df$comment <- dfj$comment

df$comment <- stringr::str_trim(df$comment)
n <- gsub("\\.", " ", df$scientific_name)
df <- df %>%
  add_column(full_name = n, .after = "scientific_name")

#write.csv(df, "outputs/fish/fish_details.csv", row.names = FALSE, fileEncoding = "UTF-8-BOM")
readr::write_excel_csv(df, "outputs/fish/fish_details.csv")

# Invertebrates
df <- read.csv("outputs/invertebrates/invertebrates_details.csv")
df_int <- read.csv("inputs/invertebrates_internal_flag.csv")
df_int$scientific_name <- gsub(" ", ".", df_int$scientific_name)
#df_int$internal_flag[df_int$internal_flag == ""] <- "green"

df1 <- df %>% select(scientific_name, flag, comment)
df2 <- df_int

dfj <- left_join(df1, df2, by = "scientific_name")

dfj <- dfj |> unite("comment", c("internal_comment", "comment"), 
                    remove = FALSE, na.rm = TRUE, sep = " ")

dfj$internal_flag[dfj$internal_flag == ""] <- "green"
dfj$internal_flag[is.na(dfj$internal_flag)] <- "green"

flag_levels <- c("red", "yellow", "green")

f1 <- factor(dfj$flag, levels = flag_levels, ordered = TRUE)
f2 <- factor(dfj$internal_flag, levels = flag_levels, ordered = TRUE)

dfj$combined_flag <- as.character(pmin(f1, f2))

df$flag <- dfj$combined_flag
df$comment <- dfj$comment

df$comment <- stringr::str_trim(df$comment)
n <- gsub("\\.", " ", df$scientific_name)
df <- df %>%
  add_column(full_name = n, .after = "scientific_name")

#write.csv(df, "outputs/invertebrates/invertebrates_details.csv", row.names = FALSE, fileEncoding = "UTF-8-BOM")
readr::write_excel_csv(df, "outputs/invertebrates/invertebrates_details.csv")

# Macrophytes
df <- read.csv("outputs/macrophytes/macrophytes_details.csv")
df_int <- read.csv("inputs/macrophytes_internal_flag.csv")
df_int$scientific_name <- gsub(" ", ".", df_int$scientific_name)
#df_int$internal_flag[df_int$internal_flag == ""] <- "green"

df1 <- df %>% select(scientific_name, flag, comment)
df2 <- df_int

dfj <- left_join(df1, df2, by = "scientific_name")

dfj <- dfj |> unite("comment", c("internal_comment", "comment"), 
                    remove = FALSE, na.rm = TRUE, sep = " ")

dfj$internal_flag[dfj$internal_flag == ""] <- "green"
dfj$internal_flag[is.na(dfj$internal_flag)] <- "green"

flag_levels <- c("red", "yellow", "green")

f1 <- factor(dfj$flag, levels = flag_levels, ordered = TRUE)
f2 <- factor(dfj$internal_flag, levels = flag_levels, ordered = TRUE)

dfj$combined_flag <- as.character(pmin(f1, f2))

df$flag <- dfj$combined_flag
df$comment <- dfj$comment

df$comment <- stringr::str_trim(df$comment)
n <- gsub("\\.", " ", df$scientific_name)
df <- df %>%
  add_column(full_name = n, .after = "scientific_name")

# Fix some of the macrophyte full names
df$full_name[df$scientific_name == "Coccotylus.Phyllophora"] <- "Coccotylus truncatus/Phyllophora pseudoceranoides"
df$full_name[df$scientific_name == "Fucus.vesiculosus.radicans"] <- "Fucus vesiculosus/Fucus radicans"
df$full_name[df$scientific_name == "Ectocarpus.Pylaiella"] <- "Ectocarpus siliculosus/Pylaiella littoralis"

#write.csv(df, "outputs/macrophytes/macrophytes_details.csv", row.names = FALSE, fileEncoding = "UTF-8-BOM")
readr::write_excel_csv(df, "outputs/macrophytes/macrophytes_details.csv")

#--------------------------------------------------------------#
############ Final adjustments ############
dff <- read.csv("outputs/fish/fish_details.csv")
dfm <- read.csv("outputs/macrophytes/macrophytes_details.csv")
dfi <- read.csv("outputs/invertebrates/invertebrates_details.csv")

## Remove ID and name columns (no longer needed)
#dff <- dff %>% select(-c(common_name, AphiaID, GBIF_speciesKey, quantity, run_initiated, run_completed))
#dff <- dff %>% add_column(species_group = "fish", .after = "full_name")
#
#dfm <- dfm %>% select(-c(common_name, AphiaID, GBIF_speciesKey, quantity, run_initiated, run_completed))
#dfm <- dfm %>% add_column(species_group = "macrophytes", .after = "full_name")
#
#dfi <- dfi %>% select(-c(common_name, AphiaID, GBIF_speciesKey, quantity, run_initiated, run_completed))
#dfi <- dfi %>% add_column(species_group = "invertebrates", .after = "full_name")

## Add species group column
dff <- dff %>% add_column(species_group = "fish", .after = "full_name")
dfm <- dfm %>% add_column(species_group = "macrophytes", .after = "full_name")
dfi <- dfi %>% add_column(species_group = "invertebrates", .after = "full_name")

## Add mobility
dff <- dff %>% add_column(mobility = "mobile", .after = "species_group")
dfm <- dfm %>% add_column(mobility = "non-mobile", .after = "species_group")
dfi <- dfi %>% add_column(mobility = "non-mobile", .after = "species_group")

# Specify mobile invertebrates
# Mobility of some invertebrates is somewhat debatable, but here defined as 
# non-mobile if sedentary enough to form an important part of the benthic community
dfi$mobility[dfi$run_name == "run3_fish"] <- "mobile"
#dfi$mobility[dfi$scientific_name == "Asterias.rubens"] <- "mobile"

## Export
readr::write_excel_csv(dff, "outputs/fish/fish_details.csv")
readr::write_excel_csv(dfm, "outputs/macrophytes/macrophytes_details.csv")
readr::write_excel_csv(dfi, "outputs/invertebrates/invertebrates_details.csv")

df <- rbind(dff, dfm, dfi)
readr::write_excel_csv(df, "outputs/WP3_species_list.csv")

