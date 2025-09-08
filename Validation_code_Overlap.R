################################################################################
# Title: Validation of Species Distribution Models Across Scenarios
# Author: Quinto Juma Meltus
# Purpose:
#   This script extracts model prediction values from raster layers at species
#   occurrence points, evaluates prediction performance across multiple
#   scenarios, and summarizes metrics such as precision, recall, F1 score, and
#   accuracy. Results are saved to CSV and optionally visualized.
#
# Notes:
#   - Replace "path/to/..." with your actual file paths.
#   - Requires R (≥ 4.2.0) and packages: terra, sp, tidyr, ggplot2
################################################################################

# --- Load Required Packages ---------------------------------------------------
library(terra)     # Raster handling
library(sp)        # Spatial points (legacy, used for compatibility)
library(tidyr)     # Data reshaping
library(ggplot2)   # Visualization

# --- Define Input and Output Paths --------------------------------------------
raster_files <- list(
  scenario_1 = "path/to/niche_sc1.tif",
  scenario_2 = "path/to/niche_sc2.tif",
  scenario_3 = "path/to/niche_sc3.tif"
)

csv_file    <- "path/to/validation_data.csv"       # Input species validation data
output_csv  <- "path/to/validation_summary.csv"    # Output metrics summary

# --- Load Species Occurrence Data ---------------------------------------------
species_points <- read.csv(csv_file)

# Convert to SpatialPointsDataFrame (requires 'lon' and 'lat' columns)
coordinates(species_points) <- ~ lon + lat

# Assign WGS84 CRS (EPSG:4326)
proj4string(species_points) <- CRS("EPSG:4326")

# --- Load Raster Scenarios ----------------------------------------------------
raster_scenario_1 <- rast(raster_files$scenario_1)
raster_scenario_2 <- rast(raster_files$scenario_2)
raster_scenario_3 <- rast(raster_files$scenario_3)

# Align CRS if needed
if (crs(species_points) != crs(raster_scenario_1)) {
  species_points <- spTransform(species_points, crs(raster_scenario_1))
}

# --- Extract Raster Values at Occurrence Points -------------------------------
species_points$scenario_1 <- extract(raster_scenario_1, species_points)
species_points$scenario_2 <- extract(raster_scenario_2, species_points)
species_points$scenario_3 <- extract(raster_scenario_3, species_points)

# Save updated points with extracted values
write.csv(as.data.frame(species_points), 
          "path/to/validation_points_with_values.csv", 
          row.names = FALSE)

# --- Define Evaluation Metrics Function ---------------------------------------
calculate_metrics <- function(data, column_name) {
  # Define counts
  true_positive  <- sum(data[[column_name]] == 3, na.rm = TRUE)
  false_positive <- sum(data[[column_name]] %in% c(0, 1, 2), na.rm = TRUE)
  false_negative <- sum(data[[column_name]] %in% c(0, 1, 2) & data$actual == 3, na.rm = TRUE)
  true_negative  <- sum(data[[column_name]] == 0, na.rm = TRUE)
  
  # Calculate metrics
  precision <- true_positive / (true_positive + false_positive)
  recall    <- true_positive / (true_positive + false_negative)
  f1_score  <- 2 * ((precision * recall) / (precision + recall))
  accuracy  <- (true_positive + true_negative) / nrow(data)
  
  # Handle division by zero
  precision <- ifelse(is.nan(precision), 0, precision)
  recall    <- ifelse(is.nan(recall), 0, recall)
  f1_score  <- ifelse(is.nan(f1_score), 0, f1_score)
  
  # Return tidy data frame
  return(data.frame(
    Metric = c("True Positives", "False Positives", "False Negatives", "True Negatives", 
               "Precision", "Recall", "F1 Score", "Accuracy"),
    Value  = c(true_positive, false_positive, false_negative, true_negative, 
               precision, recall, f1_score, accuracy)
  ))
}

# --- Compute Metrics for Each Scenario ----------------------------------------
metrics_scenario_1 <- calculate_metrics(species_points, "scenario_1")
metrics_scenario_2 <- calculate_metrics(species_points, "scenario_2")
metrics_scenario_3 <- calculate_metrics(species_points, "scenario_3")

# Combine results
all_metrics <- rbind(
  cbind(Scenario = "Scenario 1", metrics_scenario_1),
  cbind(Scenario = "Scenario 2", metrics_scenario_2),
  cbind(Scenario = "Scenario 3", metrics_scenario_3)
)

# Save summary metrics
write.csv(all_metrics, output_csv, row.names = FALSE)

# --- Optional: Visualization --------------------------------------------------
# Convert metrics to long format
metrics_long <- gather(all_metrics, key = "Variable", value = "Value", -Scenario, -Metric)

# Plot comparison
ggplot(metrics_long, aes(x = Metric, y = Value, fill = Scenario)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Validation Metrics Across Scenarios",
       x = "Metric", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# --- Reproducibility Information ----------------------------------------------
sessionInfo()
