library(vroom)

mod_id <- "run1_macrophytes" # Unique ID for model runs
spec_list <- read.csv("spec_list_macrophytes.csv")

start_dir <- paste0("outputs_extra/", mod_id, "/_run_started")
finish_dir <- paste0("outputs_extra/", mod_id, "/_run_finished")

l <- list.files(start_dir)
ls <- list.files(start_dir, full.names = TRUE)
lf <- list.files(finish_dir, full.names = TRUE)
lf <- lf[!grepl("_details.csv", lf, fixed = TRUE)]
spec_finish <- vroom(lf)

spec_list$run_initiated[spec_list$scientific_name %in% l] <- 1
m1 <- match(spec_list$scientific_name, spec_finish$scientific_name)
m2 <- !is.na(m1)
spec_list[m2, ] <- spec_finish[m1[m2], ]


#######
# Remove species files for re-run

mod_id <- "run2_macrophytes" # Unique ID for model runs
extra_dir <- paste0("../models/outputs_extra/", mod_id)
out_dir <- paste0("../models/outputs/", mod_id)
start_dir <- paste0("outputs_extra/", mod_id, "/_run_started")
finish_dir <- paste0("outputs_extra/", mod_id, "/_run_finished")

spec_list <- read.csv(paste0("outputs_extra/", mod_id, "/_run_finished/_", mod_id, "_details.csv"))

spec_rem <- species <- c(
  "Acrosiphonia.arcta",
  "Ceramium.virgatum",
  "Chaetomorpha.linum",
  "Elodea.nuttallii",
  "Eudesme.virescens",
  "Fucus.vesiculosus.radicans",
  "Halosiphon.tomentosus",
  "Ranunculus.peltatus.peltatus",
  "Ruppia.maritima",
  "Spermothamnion.repens",
  "Spongomorpha.aeruginosa",
  "Stictyosiphon.tortilis",
  "Tolypella.nidifica",
  "Zannichellia.major",
  "Zannichellia.palustris"
)

spec_rem_short <- spec_rem
# shorten name if too long (issues writing files)
spec_rem_short <- ifelse(
    nchar(spec_rem_short) > 25,
    substr(spec_rem_short, 1, 25),
    spec_rem_short) 

for(i in 1:length(spec_rem)){
  spec_full <- spec_rem[i]
  spec_short <- spec_rem_short[i]
  
  f <- paste0(out_dir, "/", spec_short)
  unlink(f, recursive = TRUE)
  
  f <- paste0(extra_dir, "/", spec_short)
  unlink(f, recursive = TRUE)
  
  f <- paste0(extra_dir, "/_run_started/", spec_full)
  unlink(f, recursive = TRUE)
  
  f <- paste0(extra_dir, "/_run_finished/", spec_full, ".csv")
  unlink(f, recursive = TRUE)
  
}


