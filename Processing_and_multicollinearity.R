# ===========================================================
# Title: Environmental Raster Preparation, Extraction, and Correlation Analysis
# Author: Quinto Juma Meltus
# Affiliations:
#   1. International Centre of Insect Physiology and Ecology (icipe),
#      P.O. Box 30772 00100, Nairobi, Kenya
#   2. School of Agricultural, Earth, and Environmental Sciences,
#      University of KwaZulu-Natal, Pietermaritzburg 3209, South Africa
#
# Description:
# This script prepares environmental predictor rasters for species
# distribution modeling (SDM) by:
#   - Reprojecting/resampling rasters to a common grid
#   - Cropping/masking to a study area mask
#   - Extracting raster values to point shapefiles (e.g., pest occurrences)
#   - Performing correlation analyses (Pearson, Spearman), VIF, and clustering
#   - Generating correlation plots and dendrograms
#
# The script is generic (no hard-coded local working directories).
# Edit the "USER PARAMETERS" section and run.

# ---------------------------
# 0. Install / load packages
# ---------------------------
required_pkgs <- c(
  "terra", "sf", "raster", "rgdal", "ggplot2", "corrplot",
  "usdm", "car", "reshape2", "ggdendro", "dplyr"
)

install_if_missing <- function(pkgs){
  missing <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
  if(length(missing)) install.packages(missing, dependencies = TRUE)
}
install_if_missing(required_pkgs)

# load libraries
library(terra)     # raster/vector operations (preferred)
library(sf)        # vector handling
library(raster)    # legacy raster functions (some convenience)
library(ggplot2)
library(corrplot)
library(usdm)      # VIF (vifstep, vifcor)
library(car)       # vif (for models)
library(reshape2)
library(ggdendro)
library(dplyr)

# ---------------------------
# USER PARAMETERS (edit here)
# ---------------------------
params <- list(
  env_dir = "path/to/environmental_rasters",        # folder with .tif rasters
  mask_raster = "path/to/mask_or_template.tif",     # template for CRS/resolution/extent
  pest_shapefiles = list(                            # list of shapefile paths (pests)
    busseola = "path/to/busseola_fusca_points.shp",
    chilo    = "path/to/chilo_partellus_points.shp"
  ),
  output_dir = "path/to/output_folder",             # where to save outputs
  categorical_pattern = "landcover|lulc|class|categor", # regex for categorical rasters
  resample_continuous = "bilinear",  # options: "bilinear", "ngb" (nearest neighbor)
  resample_categorical = "ngb",      # typically "ngb"
  overwrite = TRUE,
  save_plots = TRUE,
  correlation_methods = c("pearson", "spearman"),  # methods to compute
  vif_threshold = 10                               # threshold to flag high VIF
)

# Create output directory if missing
if(!dir.exists(params$output_dir)) dir.create(params$output_dir, recursive = TRUE)

# ---------------------------
# Helper functions
# ---------------------------

# list rasters in folder
list_rasters <- function(dir){
  files <- list.files(dir, pattern = "\\.tif$|\\.tiff$", full.names = TRUE, recursive = TRUE)
  if(length(files) == 0) stop("No raster files found in env_dir")
  files
}

# Detect categorical rasters by filename pattern
is_categorical_raster <- function(filename, pattern){
  grepl(pattern, basename(filename), ignore.case = TRUE)
}

# load mask/template raster and ensure it's a terra raster
load_template <- function(path){
  if(is.null(path) || !file.exists(path)) stop("Mask/template raster not found.")
  r <- terra::rast(path)
  r
}

# Reproject & resample a single raster to match template
align_raster_to_template <- function(raster_path, template, method = "bilinear", is_cat = FALSE, outpath = NULL, overwrite = FALSE){
  r <- terra::rast(raster_path)
  # if CRS differ, project to template
  if(!terra::crs(r) == terra::crs(template)){
    r <- terra::project(r, terra::crs(template))
  }
  # resample / align
  # use terra::resample with method
  resampled <- terra::resample(r, template, method = method)
  # mask / crop to template extent
  resampled <- terra::mask(resampled, template)
  # write if outpath provided
  if(!is.null(outpath)){
    terra::writeRaster(resampled, filename = outpath, overwrite = overwrite)
  }
  return(resampled)
}

# Extract raster values at point shapefile(s) and return dataframe
extract_to_points <- function(rast_stack, points_shp_path, id_field = NULL){
  pts_sf <- st_read(points_shp_path, quiet = TRUE)
  # transform points to raster crs
  pts_sf <- st_transform(pts_sf, crs(rast_stack))
  # terra extract prefers SpatVector
  pts_sv <- terra::vect(pts_sf)
  vals <- terra::extract(rast_stack, pts_sv)
  # combine with attributes
  out <- cbind(st_drop_geometry(pts_sf), vals)
  return(out)
}

# Correlation matrix and plots
run_correlation_analysis <- function(df, methods = c("pearson", "spearman"), outprefix, save = TRUE){
  numeric_df <- df[, sapply(df, is.numeric), drop = FALSE]
  if(ncol(numeric_df) < 2) stop("Need at least 2 numeric predictors for correlation.")
  results <- list()
  for(m in methods){
    corr_mat <- cor(numeric_df, use = "pairwise.complete.obs", method = m)
    results[[m]] <- corr_mat
    if(save){
      png(filename = file.path(params$output_dir, paste0(outprefix, "_corr_", m, ".png")), width = 1200, height = 1000, res = 150)
      corrplot::corrplot(corr_mat, method = "color", type = "upper", order = "hclust",
                         tl.cex = 0.8, tl.col = "black", addCoef.col = "black", number.cex = 0.6,
                         title = paste("Correlation (", m, ")", sep = ""))
      dev.off()
    }
  }
  return(results)
}

# VIF using usdm::vifstep (works on raster stack or dataframe)
run_vif_analysis <- function(df_or_raster, threshold = 10, outprefix, save = TRUE){
  # usdm expects either a data.frame or raster stack (RasterStack). Convert if terra SpatRaster
  if(inherits(df_or_raster, "SpatRaster")) {
    # convert to raster::stack
    rs <- raster::stack(as.list(df_or_raster))
    vs <- usdm::vifstep(rs, th = threshold)
    vif_table <- vs@results
  } else if(is.data.frame(df_or_raster)){
    vs <- usdm::vifstep(df_or_raster, th = threshold)
    vif_table <- vs@results
  } else {
    stop("df_or_raster must be a data.frame or SpatRaster")
  }
  if(save){
    write.csv(vif_table, file = file.path(params$output_dir, paste0(outprefix, "_vif_results.csv")), row.names = FALSE)
  }
  return(list(vif = vs, table = vif_table))
}

# Dendrogram based on correlation distance
plot_dendrogram <- function(corr_mat, outpng, title = "Dendrogram"){
  dist_mat <- as.dist(1 - abs(corr_mat))
  hc <- hclust(dist_mat, method = "average")
  # ggplot style
  png(outpng, width = 1200, height = 800, res = 150)
  plot(hc, main = title, xlab = "", sub = "")
  dev.off()
  return(hc)
}

# ---------------------------
# MAIN PROCESS
# ---------------------------
message("Starting processing...")

# 1. List rasters and template
rasters <- list_rasters(params$env_dir)
template <- load_template(params$mask_raster)

# 2. Prepare an aligned raster stack and save each aligned raster
aligned_list <- list()
for(rpath in rasters){
  cat("Processing:", basename(rpath), "...\n")
  is_cat <- is_categorical_raster(rpath, params$categorical_pattern)
  method <- ifelse(is_cat, params$resample_categorical, params$resample_continuous)
  outname <- paste0(tools::file_path_sans_ext(basename(rpath)), "_aligned.tif")
  outpath <- file.path(params$output_dir, outname)
  aligned <- align_raster_to_template(rpath, template, method = method, is_cat = is_cat, outpath = outpath, overwrite = params$overwrite)
  aligned_list[[basename(rpath)]] <- aligned
}

# combine into a SpatRaster (terra)
aligned_stack <- terra::rast(aligned_list)
names(aligned_stack) <- make.names(names(aligned_stack), unique = TRUE)

# Save combined stack (optional)
stack_outfile <- file.path(params$output_dir, "aligned_env_stack.tif")
terra::writeRaster(aligned_stack, filename = stack_outfile, overwrite = params$overwrite)

# 3. Extract raster values to pest points and save CSVs
extraction_results <- list()
for(pest_name in names(params$pest_shapefiles)){
  shp_path <- params$pest_shapefiles[[pest_name]]
  if(!file.exists(shp_path)){
    warning("Shapefile not found for: ", pest_name, " (", shp_path, "). Skipping.")
    next
  }
  cat("Extracting to points for:", pest_name, "...\n")
  df <- extract_to_points(aligned_stack, shp_path)
  csv_out <- file.path(params$output_dir, paste0("extracted_values_", pest_name, ".csv"))
  write.csv(df, csv_out, row.names = FALSE)
  extraction_results[[pest_name]] <- df
}

# 4. Correlation analysis and plots for each pest (and combined if desired)
for(pest_name in names(extraction_results)){
  df <- extraction_results[[pest_name]]
  # drop non-predictor columns (try to keep only numeric columns)
  numeric_df <- df %>% select(where(is.numeric))
  # remove ID columns or coordinates if present (assume first 1-3 non-predictor columns)
  if(ncol(numeric_df) < 2){
    warning("Not enough numeric columns for correlation for ", pest_name)
    next
  }
  # run correlation
  corr_results <- run_correlation_analysis(numeric_df, methods = params$correlation_methods,
                                           outprefix = paste0("corr_", pest_name), save = params$save_plots)
  # VIF
  vif_res <- run_vif_analysis(numeric_df, threshold = params$vif_threshold,
                              outprefix = paste0("vif_", pest_name), save = TRUE)
  # dendrogram from Pearson (if available)
  if("pearson" %in% names(corr_results)){
    dend_out <- file.path(params$output_dir, paste0("dendrogram_", pest_name, ".png"))
    plot_dendrogram(corr_results$pearson, outpng = dend_out, title = paste("Dendrogram (Pearson) -", pest_name))
  } else {
    # fallback to Spearman if pearson not present
    dend_out <- file.path(params$output_dir, paste0("dendrogram_", pest_name, ".png"))
    plot_dendrogram(corr_results[[1]], outpng = dend_out, title = paste("Dendrogram -", pest_name))
  }
}

# 5. Combined analysis across both pests (if both exist)
if(all(c("busseola","chilo") %in% names(extraction_results))){
  combined_df <- bind_rows(extraction_results$busseola, extraction_results$chilo, .id = "pest_source")
  numeric_combined <- combined_df %>% select(where(is.numeric))
  # run correlations and VIF
  run_correlation_analysis(numeric_combined, methods = params$correlation_methods, outprefix = "corr_combined", save = params$save_plots)
  run_vif_analysis(numeric_combined, threshold = params$vif_threshold, outprefix = "vif_combined", save = TRUE)
  # dendrogram
  corr_mat_comb <- cor(numeric_combined, use = "pairwise.complete.obs", method = "pearson")
  plot_dendrogram(corr_mat_comb, outpng = file.path(params$output_dir, "dendrogram_combined.png"), title = "Dendrogram (Combined)")
}

message("All done. Outputs are in: ", normalizePath(params$output_dir))

# Optionally print short summary
list.files(params$output_dir, pattern = "\\.(csv|png|tif)$", full.names = TRUE)



