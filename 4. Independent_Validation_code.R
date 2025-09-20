############################################################
# Title: Raster-Based Validation of Species Distribution Models
# Author: Quinto Juma Meltus
# Date: YYYY-MM-DD
# Description:
#   This script validates predicted presence–absence maps 
#   of two stemborers (Busseola fusca and Chilo partellus) 
#   against independent validation points. It extracts 
#   raster predictions, compares them to observed data, 
#   and computes evaluation metrics (Accuracy, Recall, 
#   Precision, F1-score).
#
# Notes:
#   - Replace "path/to/..." with your own file locations.
#   - Validation data must include columns: 
#       Longitude, Latitude, overlap (binary observed value).
#
# Required R packages:
#   terra, readr, dplyr
#
############################################################

# ----------------------------
# Load required libraries
# ----------------------------
library(terra)   # modern raster handling
library(readr)   # for CSV input/output
library(dplyr)   # data manipulation

# ----------------------------
# Input file paths
# ----------------------------
raster1_path <- "path/to/busseola_fusca_Current_presence_absence.tif"
raster2_path <- "path/to/Chilo_partellus_Current_presence_absence.tif"
csv_path     <- "path/to/new_val_data.csv"

# ----------------------------
# Load data
# ----------------------------
r1 <- rast(raster1_path)  # Busseola fusca
r2 <- rast(raster2_path)  # Chilo partellus

points_df <- read_csv(csv_path)  # expects: Longitude, Latitude, overlap
coords <- points_df[, c("Longitude", "Latitude")]

# ----------------------------
# Extract raster predictions at point locations
# ----------------------------
vals_r1 <- terra::extract(r1, coords)[,2]  # second column = values
vals_r2 <- terra::extract(r2, coords)[,2]

# ----------------------------
# Combine with observed overlap
# ----------------------------
df <- points_df %>%
  mutate(
    pred_busseola = vals_r1,
    pred_chilo    = vals_r2,
    observed      = overlap
  ) %>%
  filter(!is.na(pred_busseola) & !is.na(pred_chilo) & !is.na(observed))

# ----------------------------
# Define evaluation metrics function
# ----------------------------
compute_metrics <- function(predicted, observed) {
  TP <- sum(predicted == 1 & observed == 1)
  FP <- sum(predicted == 1 & observed == 0)
  FN <- sum(predicted == 0 & observed == 1)
  TN <- sum(predicted == 0 & observed == 0)
  
  Accuracy  <- (TP + TN) / (TP + FP + FN + TN)
  Recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  Precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA
  F1        <- if (!is.na(Recall) && !is.na(Precision) && (Recall + Precision) > 0) {
    2 * (Precision * Recall) / (Precision + Recall)
  } else {
    NA
  }
  
  tibble(
    True_Positives   = TP,
    False_Positives  = FP,
    False_Negatives  = FN,
    True_Negatives   = TN,
    Accuracy         = round(Accuracy, 3),
    Recall           = round(Recall, 3),
    Precision        = round(Precision, 3),
    F1_Score         = round(F1, 3)
  )
}

# ----------------------------
# Compute metrics for each species
# ----------------------------
metrics_busseola <- compute_metrics(df$pred_busseola, df$observed) %>%
  mutate(Species = "Busseola fusca")

metrics_chilo <- compute_metrics(df$pred_chilo, df$observed) %>%
  mutate(Species = "Chilo partellus")

# ----------------------------
# Combine and export results
# ----------------------------
all_metrics <- bind_rows(metrics_busseola, metrics_chilo) %>%
  select(Species, everything())

# Save results to CSV
write_csv(all_metrics, "validation_metrics_by_species.csv")

# Print to console
print(all_metrics)

############################################################
# End of Script
############################################################
