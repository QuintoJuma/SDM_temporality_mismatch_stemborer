# ================================================================
# Title:     SDM Workflow: Ensemble modeling + niche similarity
# Author:    Quinto Juma Meltus
# Date:      2025-09-08
# Purpose:   Reproducible script for species distribution modeling (SDM)
#            using the 'sdm' package (ensemble of RF, GLM, BRT),
#            performing ensemble prediction, thresholding, niche
#            similarity, and validation.
# Notes:     - Replace placeholder paths below with your local paths.
#            - This script uses 'terra' for raster I/O and 'sf' for
#              vector I/O (more modern than 'sp' / 'rgdal').
#            - If you need to run MaxEnt, set 'maxent_path' to the
#              location of maxent.jar on your machine.
# ================================================================

# -------------------------
# 0. Setup: parameters
# -------------------------
# (Replace these generic paths with your real file locations)
path_predictors_dir <- "path/to/predictors_tif_dir"       # e.g., "data/predictors/"
path_shapefiles_dir  <- "path/to/shapefiles_dir"          # e.g., "data/shapefiles/"
path_output_dir      <- "path/to/output_dir"              # e.g., "results/"
maxent_path          <- "path/to/maxent.jar"              # e.g., "tools/maxent/maxent.jar"
path_validation_csv  <- "path/to/validation_data.csv"     # optional validation points CSV

# Create output dir if missing
if (!dir.exists(path_output_dir)) dir.create(path_output_dir, recursive = TRUE)

# Reproducibility
set.seed(42)

# -------------------------
# 1. Package management
# -------------------------
required_pkgs <- c(
  "terra", "sf", "sdm", "dismo", "rJava", "ggplot2", "dplyr"
)

# Install missing packages (CRAN). For publication scripts we check then stop with instruction.
install_if_missing <- function(pkgs){
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if(length(missing)){
    message("The following packages are missing: ", paste(missing, collapse = ", "))
    message("Please install them before running this script, e.g.:")
    message("install.packages(c('", paste(missing, collapse = "', '"), "'))")
    stop("Missing packages. Aborting.")
  }
}
install_if_missing(required_pkgs)

# Load libraries (quietly)
library(terra)   # raster I/O and processing
library(sf)      # vector I/O and transforms
library(sdm)     # species distribution modelling
library(dismo)   # utilities (and MaxEnt integration)
library(rJava)   # for MaxEnt via dismo/sdm if needed
library(ggplot2)
library(dplyr)

# Print session info snippet for reproducibility (can be extended)
message("R version: ", R.version.string)
message("Package versions (selected): sdm=", packageVersion("sdm"),
        ", terra=", packageVersion("terra"), ", sf=", packageVersion("sf"))

# -------------------------
# 2. Java / MaxEnt initialization (optional)
# -------------------------
# If you use MaxEnt, provide maxent_path; otherwise skip. Do NOT include system-specific paths in the script.
if(file.exists(maxent_path)){
  # If Java home needs to be changed, set JAVA_HOME before rJava loads.
  # Sys.setenv(JAVA_HOME = "path/to/java")  # uncomment & set if needed
  # Initialize rJava safely
  tryCatch({
    .jinit()
    java_ver <- .jcall("java/lang/System", "S", "getProperty", "java.version")
    message("rJava initialized. Java version: ", java_ver)
  }, error = function(e){
    warning("Unable to initialize rJava. MaxEnt-related models may fail.\n", e$message)
  })
  message("MaxEnt JAR found at: ", maxent_path)
} else {
  message("MaxEnt JAR not found at the configured path. If you will use MaxEnt, set 'maxent_path' correctly.")
}

# -------------------------
# 3. Load predictors (rasters)
# -------------------------
# We expect many .tif files in path_predictors_dir (e.g., bioclim layers, topography, ...)
tif_files <- list.files(path = path_predictors_dir, pattern = "\\.tif$", full.names = TRUE)
if(length(tif_files) == 0) stop("No .tif predictor files found in path_predictors_dir (", path_predictors_dir, ")")

message("Loading predictors (this may take a moment)...")
predictors <- rast(tif_files)  # terra::rast reads a raster stack/collection
names(predictors) <- make.names(names(predictors), unique = TRUE) # sanitize names
plot(predictors[[1]], main = "Example predictor (first layer)")

# -------------------------
# 4. Load species occurrence data (shapefiles or CSV)
# -------------------------
# The script expects two species shapefiles in path_shapefiles_dir:
#  - Busseola fusca
#  - Chilo partellus
# If you have CSVs, you can read them with read.csv() and convert to sf.

shp_files <- list.files(path_shapefiles_dir, pattern = "\\.shp$", full.names = TRUE)
if(length(shp_files) == 0) stop("No shapefiles found in path_shapefiles_dir. Provide presence points as shapefiles or adapt the script to use CSVs.")

# Helper: read shapefile by name (generic)
read_named_shp <- function(base_name){
  f <- grep(base_name, shp_files, value = TRUE, ignore.case = TRUE)
  if(length(f) == 0) stop("Shapefile containing '", base_name, "' not found in ", path_shapefiles_dir)
  st_read(f[1], quiet = TRUE)
}

# Read species data: adapt names if your files differ
busseola_sf <- read_named_shp("busseola_fusca")   # Busseola fusca occurrences
chilo_sf     <- read_named_shp("chilo_partellus")   # Chilo partellus occurrences

# Ensure geometries are points
if(!all(st_geometry_type(busseola_sf) %in% c("POINT", "MULTIPOINT"))) stop("Busseola shapefile must contain point geometries.")
if(!all(st_geometry_type(chilo_sf) %in% c("POINT", "MULTIPOINT"))) stop("Chilo shapefile must contain point geometries.")

# Reproject species data to the predictor CRS if necessary
pred_crs <- crs(predictors)
busseola_sf <- st_transform(busseola_sf, crs = pred_crs)
chilo_sf    <- st_transform(chilo_sf, crs = pred_crs)

# -------------------------
# 5. Prepare occurrence tables for sdm::sdmData
# -------------------------
# sdm::sdmData expects a column named 'Occ' with 1 = presence, 0 = absence (or presence-only with background)
# If your shapefiles contain presence records only, create Occ = 1.
prep_occ_table <- function(sf_obj, occ_col = "Occ"){
  df <- sf_obj %>% st_coordinates() %>% as.data.frame()
  # Preserve any attribute columns if present
  attrs <- st_drop_geometry(sf_obj)
  combined <- cbind(attrs, df)
  # Ensure Occ column exists: if not present, create it as 1 (presence)
  if(!("Occ" %in% names(combined))){
    combined$Occ <- 1
  } else {
    # convert to numeric 0/1 if factor/character
    combined$Occ <- as.numeric(as.character(combined$Occ))
  }
  # Standardize coordinate names for sdmData formulas (x, y or long, lat)
  names(combined)[(ncol(combined)-1):ncol(combined)] <- c("x", "y")
  combined
}

busseola_df <- prep_occ_table(busseola_sf)
chilo_df     <- prep_occ_table(chilo_sf)

# Quick checks
message("Number of Busseola presence records: ", nrow(busseola_df))
message("Number of Chilo presence records: ", nrow(chilo_df))

# -------------------------
# 6. Create sdmData objects
# -------------------------
# For presence-only data we supply background samples via bg = list(...)
# The coordinate names in the formula must match the column names created above.
sdm_data1 <- sdmData(
  formula = Occ ~ . + coords(x + y), 
  train   = busseola_df, 
  predictors = predictors,
  bg = list(n = 1000, method = 'gRandom', remove = TRUE)  # adjust n as needed
)

sdm_data2 <- sdmData(
  formula = Occ ~ . + coords(x + y), 
  train   = chilo_df, 
  predictors = predictors,
  bg = list(n = 800, method = 'gRandom', remove = TRUE)   # adjust n as needed
)

message("sdmData objects created for both species.")

# -------------------------
# 7. Model training (ensemble)
# -------------------------
# Choose modeling methods and replication settings. Adjust as needed for computation time.
methods_vec <- c("rf", "glm", "brt")   # random forest, generalized linear model, boosted regression tree

# Train Busseola models
message("Training Busseola models (this may take time)...")
t0 <- Sys.time()
model1 <- sdm(Occ ~ ., data = sdm_data1, methods = methods_vec,
              replication = "cv", test.percent = 30, cv.folds = 5)  # use 5-fold CV (adjust)
t1 <- Sys.time()
message("Busseola training time: ", round(difftime(t1,t0,units="mins"), 2), " minutes")

# Train Chilo models
message("Training Chilo models (this may take time)...")
t0 <- Sys.time()
model2 <- sdm(Occ ~ ., data = sdm_data2, methods = methods_vec,
              replication = "cv", test.percent = 30, cv.folds = 5)
t1 <- Sys.time()
message("Chilo training time: ", round(difftime(t1,t0,units="mins"), 2), " minutes")

# Save model objects (optional)
saveRDS(model1, file = file.path(path_output_dir, "model_busseola.rds"))
saveRDS(model2, file = file.path(path_output_dir, "model_chilo.rds"))

# -------------------------
# 8. Ensemble predictions
# -------------------------
# Use model ensembles weighted by AUC (or other metric). Save outputs to output_dir.
message("Creating ensemble predictions...")

ensemble_outfile1 <- file.path(path_output_dir, "ensemble_busseola.tif")
ensemble1 <- ensemble(model1, newdata = predictors, filename = ensemble_outfile1,
                      setting = list(method = "weighted", stat = "AUC"))
message("Busseola ensemble saved to: ", ensemble_outfile1)

ensemble_outfile2 <- file.path(path_output_dir, "ensemble_chilo.tif")
ensemble2 <- ensemble(model2, newdata = predictors, filename = ensemble_outfile2,
                      setting = list(method = "weighted", stat = "AUC"))
message("Chilo ensemble saved to: ", ensemble_outfile2)

# Quick visualization (first ensemble layer if multiple)
plot(ensemble1, main = "Busseola ensemble suitability")
plot(ensemble2, main = "Chilo ensemble suitability")

# -------------------------
# 9. Model evaluation and variable importance
# -------------------------
eval1 <- getEvaluation(model1)
eval2 <- getEvaluation(model2)
# Basic printout
message("Busseola model evaluation (summary):"); print(summary(eval1))
message("Chilo model evaluation (summary):"); print(summary(eval2))

# Variable importance
vi1 <- getVarImp(model1)
vi2 <- getVarImp(model2)
plot(vi1, main = "Variable importance - Busseola")
plot(vi2, main = "Variable importance - Chilo")

# -------------------------
# 10. Thresholding: presence-absence classification
# -------------------------
# Using pa() from sdm to classify ensemble predictions by optimized threshold
message("Converting ensemble suitability to presence-absence (optimized threshold).")
pa_busseola <- pa(ensemble1, model1, id = "ensemble", opt = 2) # opt=2 => max(sens+spec)
pa_chilo    <- pa(ensemble2, model2, id = "ensemble", opt = 2)

# Save PA rasters
pa_file1 <- file.path(path_output_dir, "presence_absence_busseola.tif")
writeRaster(pa_busseola, filename = pa_file1, overwrite = TRUE)
message("Busseola presence-absence raster saved to: ", pa_file1)

pa_file2 <- file.path(path_output_dir, "presence_absence_chilo.tif")
writeRaster(pa_chilo, filename = pa_file2, overwrite = TRUE)
message("Chilo presence-absence raster saved to: ", pa_file2)

# -------------------------
# 11. Niche similarity (geographic/suitability space)
# -------------------------
# sdm::nicheSimilarity can compare two suitability rasters; ensure same extent/resolution/CRS.
# If you want statistical metrics, you can compute I-statistic, Schoener's D, etc. sdm::nicheSimilarity accepts objects from sdm/ensemble output.
message("Computing niche similarity metrics between ensembles...")
ns <- tryCatch({
  niche_sim <- nicheSimilarity(ensemble1, ensemble2, stat = c("Imod", "Icor", "R"))
  niche_sim
}, error = function(e){
  warning("nicheSimilarity failed: ", e$message)
  NULL
})
if(!is.null(ns)) print(ns)

# Export a simple raster showing one similarity metric for visualization.
# Here we create a uniform raster with the Imod value as an example visualization.
if(!is.null(ns) && !is.null(ns$Imod)){
  sim_raster <- ensemble1  # copy structure
  values(sim_raster) <- as.numeric(ns$Imod[1]) # use the first Imod result
  sim_file <- file.path(path_output_dir, "niche_similarity_Imod.tif")
  writeRaster(sim_raster, sim_file, overwrite = TRUE)
  message("Niche similarity raster exported to: ", sim_file)
  plot(sim_raster, main = "Niche similarity (Imod)")
}

# -------------------------
# 12. Optional: PCA-based environmental niche analysis
# -------------------------
# If desired, compute PCA on predictor values sampled from the study area.
message("Performing PCA on predictor stack (sample-based).")
# Sample points from predictor rasters (use a manageable sample size)
samp_cells <- sample(seq_len(ncell(predictors)), size = min(5000, ncell(predictors)))
pred_vals <- terra::values(predictors)[samp_cells, , drop = FALSE]
pred_vals <- na.omit(as.data.frame(pred_vals))
pca_res <- prcomp(pred_vals, scale. = TRUE, center = TRUE)
summary(pca_res)
# Save/load if needed
saveRDS(pca_res, file = file.path(path_output_dir, "predictors_pca.rds"))

# -------------------------
# 13. Validation using independent data (optional)
# -------------------------
# If you have a CSV with validation points and a column 'actual_presence' (0/1),
# we can extract predicted suitability and compute basic metrics.
if(file.exists(path_validation_csv)){
  message("Loading validation dataset from: ", path_validation_csv)
  val_df <- read.csv(path_validation_csv)
  # Expect columns named lon, lat, actual_presence
  if(!all(c("lon", "lat", "actual_presence") %in% names(val_df))){
    warning("Validation CSV must contain columns: lon, lat, actual_presence. Skipping validation.")
  } else {
    # Convert to terra/spatVector and extract suitability from similarity raster (or ensemble)
    val_pts <- vect(val_df, geom = c("lon", "lat"), crs = pred_crs)
    # Choose the raster to extract from; here we use sim_raster if present, else ensemble1
    raster_to_use <- if(exists("sim_raster")) sim_raster else ensemble1
    extracted <- terra::extract(raster_to_use, val_pts)
    # Combine and compute metrics
    val_df$predicted_value <- extracted[,2]  # depending on extract output, adjust index
    # If predicted_value is NA, report
    val_df <- val_df %>% filter(!is.na(predicted_value))
    threshold <- median(val_df$predicted_value, na.rm = TRUE)  # example threshold; change as needed
    val_df$pred_presence <- ifelse(val_df$predicted_value >= threshold, 1, 0)
    cm <- table(val_df$actual_presence, val_df$pred_presence)
    accuracy <- sum(diag(cm))/sum(cm)
    precision <- ifelse(sum(cm[,2])==0, NA, cm[2,2]/sum(cm[,2]))
    recall <- ifelse(sum(cm[2,])==0, NA, cm[2,2]/sum(cm[2,]))
    f1 <- ifelse(is.na(precision) | is.na(recall) | (precision+recall)==0, NA, 2*(precision*recall)/(precision+recall))
    message("Validation metrics (using chosen threshold):")
    message(sprintf("Accuracy = %.3f; Precision = %.3f; Recall = %.3f; F1 = %.3f", accuracy, precision, recall, f1))
    # Save validation results
    write.csv(val_df, file = file.path(path_output_dir, "validation_results.csv"), row.names = FALSE)
    message("Validation results saved.")
  }
} else {
  message("No validation CSV provided. Skipping validation step.")
}

# -------------------------
# 14. Housekeeping and session information
# -------------------------
message("Script completed. Key outputs (in path_output_dir):")
list.files(path_output_dir, full.names = TRUE)

# Save sessionInfo for reproducibility
writeLines(capture.output(sessionInfo()), con = file.path(path_output_dir, "sessionInfo.txt"))

# End of script
# ================================================================
