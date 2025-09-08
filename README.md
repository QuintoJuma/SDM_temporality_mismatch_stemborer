🌱 SDM Temporality Mismatch in Stemborer Distribution
🧩 About the Project

This repository provides a fully reproducible workflow for analyzing Species Distribution Models (SDMs) of stemborers (Busseola fusca and Chilo partellus) in Kenya, focusing on temporal mismatches between current and future climatic scenarios.

Key highlights:

✅ Hypothesis testing for SDM performance

✅ Raster-based validation using observed occurrence points

✅ Model evaluation metrics: Accuracy, Precision, Recall, F1-score

🗂 Repository Structure
├── data/
│   ├── rasters/              # Current/future SDM projections
│   ├── shapefiles/           # Kenya administrative boundaries
│   └── validation/           # Observed species occurrence CSVs
├── scripts/
│   ├── 01_SDM_building.R    # Build SDMs and extract raster predictions
│   ├── 02_SDM_validation.R  # Compute validation metrics
│   ├── 03_Range_Shift_Analysis.R # Calculate county/national range shifts
│   └── 04_Visualization.R   # Generate plots and heatmaps
├── outputs/
│   ├── metrics/              # Validation metrics CSVs
│   ├── plots/                # PNGs and heatmaps
│   └── shift_summaries/      # County/national range shift tables
└── README.md

🔄 Workflow Diagram
flowchart TD
    A[Raw Data: Occurrence Points & Environmental Rasters] --> B[Species Distribution Models (SDMs)]
    B --> C[Raster Predictions for Each Scenario (BAU & PESS)]
    C --> D[Validation Metrics Extraction]
    D --> E[County & National Range Shift Analysis]
    E --> F[Plots: Heatmaps, Barplots, Combined Figures]
    F --> G[Summary Tables: Metrics & Area Changes]


This diagram shows the full pipeline from raw data to final outputs.

⚙️ Setup Instructions

