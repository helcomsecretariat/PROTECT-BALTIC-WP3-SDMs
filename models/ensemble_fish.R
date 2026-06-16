########## Set up reproduciple R package environment ##########
# Ensure correct working directory and package environment
# The "setup_environment.R" script must be run first to set up the package environment
options(show.error.messages = TRUE)
wd <- "C:/rdir/pb_sdms_r/models"
setwd(wd)
#library(renv)
#renv::load(wd)
#renv::deactivate()

########## Model and script settings ##########
mod_id <- "run3_fish" # Unique ID for model runs
run_name <- mod_id
predictor_dir <- "../predictors/predictors_fish"
obs_dir <- "../observations/fish"
spec_list_file <- "spec_list_fish.csv"
grid_loc <- "inputs/grid_marine_250m_v3.tif"
source("_functions/cv_partition.R")
source("_functions/calc_tss.R")

# Models to include in the ensemble
mod_methods <- c("RFd", "XGBOOST", "GAM") # Subset
n_rep <- 5 # Number of cross-validation repetitions
cv_perc <- 0.8 # Percentage of training/testing data
presence_thresh <- 0 # Minimum number of presences required to run the model
method_pred <- TRUE # Add method as a predictor variable/covariate?
method_pred_thresh <- 0 # Set to 0 to always include all methods

do_full_run <- FALSE # Whether to do a full run with all data (and only use CV reps for evaluation)
sink_log <- TRUE
mod_metric_eval <- c('TSS', 'ROC') # All evaluation metrics used
mod_metric_select <- 'TSS' # Metric for model selection and weighting
mod_metric_thresh <- 0.4 # Threshold for model inclusion in the ensemble (check threshold is suitable for selection metric)
mod_metric_bin <- 'TSS' # Metric for determining the binary presence/absence threshold

# Custom model settings
# RF custom settings
rf_mtry <- 4 # Approximation: round(sqrt(nlyr(preds)))
rf_ntree <- 1500 # Suggest 1500
rf_nodesize <- 2 # Suggest 2

# XGBOOST custom settings
xg_nrounds <- 1000 # Suggest 1000
xg_early_stopping_rounds <- 20 # Suggest 20
xg_max_depth <- 4 # Suggest 4
xg_eta <- 0.01 # Suggest 0.01

# GAM custom settings
gam_method = "REML" # Suggest "REML"
gam_select = TRUE # Suggest TRUE
gam_gamma = 1.4 # Suggest 1.4

########## Setup ##########
# Packages
library(terra)
library(biomod2)
library(caret)
library(dplyr)
library(tidyr)
library(FNN)
library(vroom)
library(viridis)
library(ggplot2)
library(patchwork)

pdf(NULL)
start_time <- Sys.time()

# Create directories
if(!dir.exists("outputs")){dir.create("outputs")} # Main outputs from biomod
if(!dir.exists("outputs_extra")){dir.create("outputs_extra")} # Desired outputs
start_dir <- paste0("outputs_extra/", mod_id, "/_run_started")
finish_dir <- paste0("outputs_extra/", mod_id, "/_run_finished")
if(!dir.exists(start_dir)){dir.create(start_dir, recursive = TRUE)}
if(!dir.exists(finish_dir)){dir.create(finish_dir, recursive = TRUE)}

# Import species list
spec_list <- read.csv(spec_list_file)
l <- list.files(start_dir)
spec_list$run_initiated[spec_list$scientific_name %in% l] <- 1
spec_miss <- spec_list %>%
  filter(run_initiated == 0) %>%
  filter(quantity > presence_thresh)

spec_full <- first(spec_miss$scientific_name)
spec_full_ <- gsub(" ", ".", spec_full)
spec <- spec_full
if(nchar(spec) > 25){spec <- substr(spec, 1, 25)} # shorten name if too long (issues writing files)
spec_ <- gsub(" ", ".", spec)
write.csv(NA, paste0(start_dir, "/", spec_full)) # Record that run was initiated in _run_started folder
spec_list$run_initiated[spec_list$scientific_name == spec_full] <- 1

# Biomod file/folder names and creation of output folders
resp_n <- gsub(" ", "_", spec)
proj_n <- mod_id
out_dir <- paste0("outputs/", mod_id)
if(!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}
extra_dir <- paste0("outputs_extra/", mod_id, "/", spec_)
if(!dir.exists(extra_dir)){dir.create(extra_dir, recursive = TRUE)}

# Save console outputs to file
if(sink_log) {
  log_file <- paste0(extra_dir, "/", spec_, "_log.txt")
  log_file <- file(log_file, "at")
  sink(log_file)
  sink(log_file, type = "message")
}

########## Observation data ##########
obs <- readRDS(paste0(obs_dir, "/", spec_full_, ".rds"))

# Convert to binary
obs$quantity[obs$quantity > 0] <- 1

# Remove spaces from method
obs$method <- gsub(" ", "_", obs$method)

########## Predictors ##########
grid <- rast(grid_loc)
l <- list.files(predictor_dir, pattern = ".tif$", full.names = TRUE)
ls <- list.files(predictor_dir, pattern = ".tif$", full.names = FALSE)
preds <- rast(l)
names(preds) <- tools::file_path_sans_ext(ls)

pred_mask <- preds*0
pred_mask <- app(pred_mask, fun = "sum")
grid <- pred_mask + 1 # Overwrite grid as all cells with predictor values

# Take random grid coordinates (due to restricted coordinates in source data)
grid_4326 <- project(grid, "EPSG:4326")
grid_not_na <- values(grid_4326)
grid_not_na <- which(!is.na(grid_not_na))
grid_not_na <- sample(grid_not_na, nrow(obs), replace = FALSE)

grid_xy <- xyFromCell(grid_4326, cell = grid_not_na)

obs$longitude <- grid_xy[,1]
obs$latitude <- grid_xy[,2]

xy <- vect(obs[, c("longitude", "latitude")], geom=c("longitude", "latitude"), crs = "EPSG:4326")
xy <- project(xy, grid)

# Select columns with selected predictors
df_preds <- obs[, colnames(obs) %in% names(preds)]

########## Add method co-variate to predictor stack ##########
# Get summary statistics about the methods
method_tab <- obs %>% group_by(method) %>%
  summarise(num_samples = n(),
            num_presence = sum(quantity),
            perc_presence = mean(quantity)*100)

# Identify method with highest prevalence
method_most <- method_tab$method[which.max(method_tab$perc_presence)]
obs$method <- as.factor(obs$method)
obs$method <- droplevels(obs$method)
method_levels <- levels(obs$method)

method_pred <- spec_list$use_method_cov[spec_list$scientific_name == spec_full_]  # Only use method as predictor if specified in the settings
if(length(method_levels) == 1){method_pred <- FALSE} # Only use method as a predictor if more than 1 method in the data

if(method_pred){
  
  dummies <- model.matrix(~ method, data = obs)
  dummies <- as.data.frame(dummies)
  dummies <- dummies[,-1, drop = FALSE]
  
  df_preds <- cbind(df_preds, dummies)
  
  
  for(i in 2:length(method_levels)){
    r <- grid
    r[] <- 0 # Defaults all method layers to zero (which represents the intercept)
    if(method_most == method_levels[i]){r[] <- 1}
    
    if(i == 2){dummy_preds <- r}
    if(i > 2){dummy_preds <- c(dummy_preds, r)}
  }
  
  names(dummy_preds) <- colnames(dummies)
  varnames(dummy_preds) <- colnames(dummies)
  
  preds <- c(preds, dummy_preds) # Add dummy predictors to predictor stack
  
}

########## Finalize observations ##########
# Remove NA points
obs_na <- complete.cases(df_preds)
obs <- obs[obs_na, ]
xy <- xy[obs_na, ]
df_preds <- df_preds[obs_na, ]
xy_pres <- xy[which(obs$quantity == 1)]
xy_abs <- xy[which(obs$quantity == 0)]

# Presence/absence map export
#clrs <- ifelse(obs$quantity == 1, "red", "black")
#temp_plot <- plot(xy, col = clrs, pch = 16, cex = 1.5)
pdf(file = paste0(extra_dir, "/", spec_, "_point_map.pdf"),
    width = 8,
    height = 8)
plot(grid, col = "lightgrey")
plot(xy_abs, col = "black", pch = 16, cex = 0.3, add = TRUE)
plot(xy_pres, col = "red", pch = 16, cex = 0.3, add = TRUE)
dev.off()

########## Model setup ##########
## Prepare model
dat <- BIOMOD_FormatingData(resp.var = obs$quantity,
                            expl.var = df_preds,
                            resp.name = resp_n,
                            dir.name = out_dir,
                            resp.xy = terra::crds(xy))

dat

# Set training/testing data
# Set seed based on species number
spec_id <-  which(spec_list$scientific_name == spec_full)
set.seed(123456 + spec_id)
calib <- cv_partition(x = obs$quantity, train_perc = cv_perc, n = n_rep)
calib_exp <- cbind(calib, obs$quantity) # Export CV sets if we need to double-check prevalence
write.csv(calib_exp, paste0(extra_dir, "/", spec_, "_CV_permutations.csv"), row.names = FALSE)

mod_cv <- bm_CrossValidation(bm.format = dat,
                             strategy = 'user.defined',
                             user.table = calib)

# Set custom settings
user.RFd <- vector("list", n_rep)
names(user.RFd) <- paste0("_allData_RUN", 1:n_rep)
user.RFd <- lapply(user.RFd, function(x) list(mtry = rf_mtry,
                                              ntree = rf_ntree,
                                              nodesize = rf_nodesize))

user.XGBOOST <- vector("list", n_rep)
names(user.XGBOOST) <- paste0("_allData_RUN", 1:n_rep)
user.XGBOOST <- lapply(user.XGBOOST, function(x) list(nrounds = xg_nrounds,
                                                      early_stopping_rounds = xg_early_stopping_rounds,
                                                      params = list(max_depth = xg_max_depth, eta = xg_eta)))

user.GAM <- vector("list", n_rep)
names(user.GAM) <- paste0("_allData_RUN", 1:n_rep)
user.GAM <- lapply(user.GAM, function(x) list(method = gam_method,
                                              select = gam_select,
                                              gamma = gam_gamma))

user_val <- list(RFd.binary.randomForest.randomForest = user.RFd,
                 XGBOOST.binary.xgboost.xgboost = user.XGBOOST,
                 GAM.binary.mgcv.gam = user.GAM)

user_opt <- bm_ModelingOptions(data.type = 'binary',
                               models = mod_methods,
                               strategy = "user.defined",
                               user.val = user_val,
                               user.base = "bigboss",
                               bm.format = dat,
                               calib.lines = mod_cv)

########## Single models ##########
mod <- BIOMOD_Modeling(bm.format = dat,
                       models = mod_methods,
                       modeling.id = mod_id,
                       CV.strategy = 'user.defined',
                       CV.user.table = mod_cv,
                       OPT.strategy = 'user.defined',
                       OPT.user = user_opt,
                       CV.do.full.models = do_full_run,
                       var.import = 1,
                       metric.eval = mod_metric_eval,
                       do.progress = FALSE)

# Model evaluation plot export
bm_plot <- bm_PlotEvalMean(mod, dataset = "validation", do.plot = FALSE)
temp_plot <- bm_plot[["plot"]]
pdf(file = paste0(extra_dir, "/", spec_, "_evaluation.pdf"),
    width = 8,
    height = 8)
print(temp_plot)
dev.off()

########## Ensemble model ##########
# Check model overfitting
mod_eval <- get_evaluations(mod)
mod_eval <- filter(mod_eval, is.na(validation) == FALSE)
mod_eval$overfit <- mod_eval$calibration - mod_eval$validation
mod_select <- mod_eval$full.name

# Run with 'all' for single ensemble
emod <- BIOMOD_EnsembleModeling(bm.mod = mod,
                                models.chosen = mod_select,
                                em.by = 'all',
                                em.algo = c('EMwmean', 'EMca'), # Weighted mean and committee average
                                metric.eval = mod_metric_eval,
                                metric.select = mod_metric_select,
                                metric.select.thresh = mod_metric_thresh, # Cutoff for inclusion/removal
                                metric.select.dataset = 'validation',
                                var.import = 2,
                                EMci.alpha = 0.05,
                                #EMwmean.decay = 'proportional' # Select one of the EMwmean.decay options
                                #EMwmean.decay = 1.5
                                EMwmean.decay = function(x){exp(5*x)})

# Save validation scores and kept models
mod_eval$kept_in_ensemble <- mod_eval$full.name %in% get_kept_models(emod)
write.csv(mod_eval, paste0(extra_dir, "/", spec_, "_validation_scores.csv"), row.names = FALSE)

# Response curves plot export (single models)
bm_plot <- bm_PlotResponseCurves(bm.out = mod, 
                                 models.chosen = get_built_models(mod),
                                 fixed.var = 'median',
                                 do.progress = FALSE,
                                 do.plot = FALSE)
temp_plot <- bm_plot[["plot"]]
pdf(file = paste0(extra_dir, "/", spec_, "_response_curves_single.pdf"),
    width = nlyr(preds)*1.4,
    height = nlyr(preds)*0.9)
print(temp_plot)
dev.off()

# Response curves plot export (ensemble)
bm_plot <- bm_PlotResponseCurves(bm.out = emod, 
                                 models.chosen = get_built_models(emod)[1],
                                 fixed.var = 'median',
                                 do.progress = FALSE,
                                 do.plot = FALSE)
temp_plot <- bm_plot[["plot"]]
pdf(file = paste0(extra_dir, "/", spec_, "_response_curves_ensemble.pdf"),
    width = nlyr(preds)*1.4,
    height = nlyr(preds)*0.9)
print(temp_plot)
dev.off()

# Variable importance plot export
bm_plot <- bm_PlotVarImpBoxplot(mod,
                                group.by = c("algo", "expl.var", "expl.var"),
                                do.plot = TRUE)

temp_plot <- bm_plot[["plot"]]
pdf(file = paste0(extra_dir, "/", spec_, "_var_importance.pdf"),
    width = nlyr(preds)*1.4,
    height = nlyr(preds)*0.9)
print(temp_plot)
dev.off()

# Save variable importance scores
var_imp <- get_variables_importance(emod)
var_imp <- var_imp %>%
  filter(algo == "EMwmean") %>%
  group_by(expl.var) %>%
  summarise(var.imp = mean(var.imp))

write.csv(var_imp, paste0(extra_dir, "/", spec_, "_var_importance.csv"), row.names = FALSE)

########## Full ensemble evaluation ##########
emod_preds <- get_predictions(emod)
emod_preds <- filter(emod_preds, algo == "EMwmean")

emod_thresh <- get_evaluations(emod)
emod_thresh <- emod_thresh$cutoff[emod_thresh$algo == "EMwmean" & emod_thresh$metric.eval == "TSS"]

for(i in 1:n_rep){
  cv <- !calib[,i] # testing data
  cv_preds <- emod_preds$pred[cv]
  cv_preds <- ifelse(cv_preds >= emod_thresh, 1, 0)
  
  temp <- calc_tss(obs = obs$quantity[cv], pred = cv_preds)
  temp <- temp$TSS
  
  if(i == 1){emod_TSS <- temp}
  if(i > 1){emod_TSS <- c(emod_TSS, temp)}
}

emod_TSS <- mean(emod_TSS)
emod_TSS_sub <- mod_eval[mod_eval$metric.eval == mod_metric_select, ]
emod_TSS_sub <- mean(emod_TSS_sub$validation)

########## Predictions (individual models) ##########
mod_proj <- BIOMOD_Projection(bm.mod = mod,
                              proj.name = proj_n,
                              new.env = preds,
                              models.chosen = mod_select,
                              metric.binary = mod_metric_bin,
                              build.clamping.mask = FALSE,
                              keep.in.memory = FALSE,
                              do.stack = FALSE, # Check to see if this help reduce memory usage
                              omit.na = FALSE
)

# Individual projections plot export
temp_plot <- plot(mod_proj, do.plot = FALSE)
pdf(file = paste0(extra_dir, "/", spec_, "_individual_projections.pdf"),
    width = n_rep*10,
    height = length(mod_methods)*5)
print(temp_plot)
dev.off()

########## Predictions (ensemble) ##########
emod_proj <- BIOMOD_EnsembleForecasting(bm.em = emod,
                                        bm.proj = mod_proj,
                                        models.chosen = 'all',
                                        metric.binary = mod_metric_bin,
                                        keep.in.memory = FALSE,
                                        do.stack = FALSE,
)

# Binary and probability TIF save to outputs_extra
em_list <- terra::unwrap(emod_proj@proj.out@link)
em_wmean <- rast(em_list[grepl("EMwmean", em_list) & !grepl("bin.tif", em_list)])
em_wmean_bin <- rast(em_list[grepl("EMwmean", em_list) & grepl(paste0(mod_metric_bin, "bin.tif"), em_list)]) # Use ensemble binary threshold
em_wmean[is.na(grid)] <- NA
em_wmean_bin[is.na(grid)] <- NA
writeRaster(em_wmean, paste0(extra_dir, "/", spec_, "_ensemble_EMwmeanBy", mod_metric_select, ".tif"), overwrite = TRUE)
writeRaster(em_wmean_bin, paste0(extra_dir, "/", spec_, "_ensemble_", mod_metric_bin, "bin.tif"), overwrite = TRUE)

# Binary and probability plot
# Ensemble projections plot export
temp_plot <- plot(emod_proj, do.plot = FALSE)
pdf(file = paste0(extra_dir, "/", spec_, "_ensemble_projections.pdf"),
    width = 15,
    height = 10)
print(temp_plot)
dev.off()

# Ensemble probability and binary plot export
df1 <- as.data.frame(em_wmean, xy = TRUE)
colnames(df1) <- c("x", "y", "z")
p1 <- ggplot(df1, aes(x, y, fill = z)) +
  geom_raster() +
  scale_fill_viridis() +
  coord_equal() +
  theme(panel.background = element_rect(fill = "grey60"), panel.grid = element_blank(),
        axis.title = element_blank(), axis.text = element_blank()) +
  ggtitle("Ensemble weighted mean")

df2 <- as.data.frame(em_wmean_bin, xy = TRUE)
colnames(df2) <- c("x", "y", "z")
p2 <- ggplot(df2, aes(x, y, fill = z)) +
  geom_raster() +
  scale_fill_viridis() +
  coord_equal() +
  theme(panel.background = element_rect(fill = "grey60"), panel.grid = element_blank(),
        axis.title = element_blank(), axis.text = element_blank()) +
  ggtitle("Ensemble binary")

p <- p1 | p2
pdf(file = paste0(extra_dir, "/", spec_, "_binary.pdf"), width = 18, height = 10)
print(p)
dev.off()

########## Confidence/uncertainty maps ##########
em_list <- emod_proj@proj.out@link

em_ca <- rast(em_list[grepl("EMca", em_list) & !grepl("bin.tif", em_list)])
em_ca <- em_ca * grid
em_ca <- ((abs(em_ca - 500)) / 500) * 100 # Confidence maps with 100 giving the highest confidence

# Export plots
# Uncertainty plot
df <- as.data.frame(em_ca, xy = TRUE)
colnames(df) <- c("x", "y", "z")
temp_plot <- ggplot(df, aes(x = x, y = y, fill = z)) +
  geom_raster() +
  scale_fill_viridis(name = "Value", option = "D", limits = c(0,100)) +
  coord_equal() +
  theme(panel.background = element_rect(fill = "grey60"), panel.grid = element_blank())

pdf(file = paste0(extra_dir, "/", spec_, "_confidence.pdf"),
    width = 15,
    height = 10)
print(temp_plot)
dev.off()

# Export uncertainty raster
writeRaster(em_ca, paste0(extra_dir, "/", spec_, "_confidence.tif"), overwrite = TRUE)

########## Binary predictions combined with confidence ##########
em_list <- terra::unwrap(emod_proj@proj.out@link)
em_wmean_bin <- rast(em_list[grepl("EMwmean", em_list) & grepl(paste0(mod_metric_bin, "bin.tif"), em_list)]) # Use ensemble binary threshold
em_bin <- em_wmean_bin*grid

em_list <- emod_proj@proj.out@link
em_ca <- rast(em_list[grepl("EMca", em_list) & !grepl("bin.tif", em_list)])
em_ca <- em_ca * grid
em_ca <- ((abs(em_ca - 500)) / 500) * 100 # Confidence maps with 100 giving the highest confidence

conf <- (em_ca >= 50) + 1
pres <- em_bin
pres[pres != 1] <- NA
pres <- (conf * pres) + 2

conf <- (em_ca >= 50) + 1
conf <- -(conf-3)
abs <- em_bin == 0
abs[abs != 1] <- NA
abs <- (conf * abs)

pres[is.na(pres)] <- 0
abs[is.na(abs)] <- 0

em_bin_conf <- pres + abs
em_bin_conf[is.na(grid)] <- NA

# Export plot
labs <- c("Absent (high confidence)", "Absent (low confidence)", "Present (low confidence)", "Present (high confidence)")
cols <- c('white','#ffffb3','#f46d43', '#a50026')

df <- as.data.frame(em_bin_conf, xy = TRUE, na.rm = TRUE)
names(df)[3] <- "class"

df$category <- factor(df$class, levels = 1:4, labels = labs)

temp_plot <- ggplot(df, aes(x, y, fill = category)) +
  geom_raster() +
  coord_equal() +
  scale_fill_manual(values = cols, name = "") +
  theme(panel.background = element_rect(fill = "grey80"), panel.grid = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())+
  ggtitle(spec_full)

pdf(file = paste0(extra_dir, "/", spec_, "_binary_confidence.pdf"),
    width = 9.5,
    height = 9)
print(temp_plot)
dev.off()

# Export raster
writeRaster(em_bin_conf, paste0(extra_dir, "/", spec_, "_binary_confidence.tif"), overwrite = TRUE)

########## Export run details ##########
# Record run time and export run details
end_time <- Sys.time()
run_time <- as.numeric(difftime(end_time, start_time, units = "hours"))

#spec_list <- read.csv("spec_list.csv")
spec_list$run_completed[spec_list$scientific_name == spec_full] <- 1
spec_list$run_time[spec_list$scientific_name == spec_full] <- run_time
spec_list$models[spec_list$scientific_name == spec_full] <- paste(mod_methods, collapse = ", ")
spec_list$num_presence[spec_list$scientific_name == spec_full] <- nrow(obs[obs$quantity == 1, ])
spec_list$num_absence[spec_list$scientific_name == spec_full] <- nrow(obs[obs$quantity == 0, ])
spec_list$run_name[spec_list$scientific_name == spec_full] <- run_name
spec_list$metric_select_and_weight <- mod_metric_select
spec_list$metric_binary <- mod_metric_bin
spec_list$mean_TSS_ensemble <- emod_TSS
spec_list$mean_TSS_all_submodels <- emod_TSS_sub

if(method_pred){
  spec_list$observation_methods <- paste(method_levels, collapse = ",")
  spec_list$method_predict_to <- method_most
}

write.csv(spec_list[spec_list$scientific_name == spec_full,], paste0(finish_dir, "/", spec_full_, ".csv"), row.names = FALSE)

# Save full run log
l <- list.files(start_dir)
ls <- list.files(start_dir, full.names = TRUE)
lf <- list.files(finish_dir, full.names = TRUE)
lf <- lf[!grepl("_details.csv", lf, fixed = TRUE)]
spec_finish <- vroom(lf)
spec_list <- read.csv(spec_list_file)
spec_list$run_initiated[spec_list$scientific_name %in% l] <- 1
m1 <- match(spec_list$scientific_name, spec_finish$scientific_name)
m2 <- !is.na(m1)
spec_list[m2, ] <- spec_finish[m1[m2], ]
write.csv(spec_list, paste0(finish_dir, "/_", mod_id, "_details.csv"), row.names = FALSE)

if(sink_log){
  sink(type = "message")
  sink()
  close(log_file)
}
