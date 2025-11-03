# 🌿 SDM Temporality Mismatch in Stemborer Distribution

### 🧬 About the Project

This repository presents a complete and reproducible workflow for assessing temporal mismatches in species distribution modeling (SDM) of *Busseola fusca* and *Chilo partellus*.  
The analysis investigates how temporal gaps between species occurrence records and climate data influence model accuracy and prediction reliability.

The workflow includes scripts for data preprocessing, full SDM modeling, hypothesis testing, and independent validation using current occurrence datasets.  
Evaluation of model performance is conducted using standard classification metrics such as True Positives, False Positives, False Negatives, True Negatives, Precision, Recall, F1-score, and Accuracy.  
All outputs; including raster layers, validation tables, and summary statistics, are structured for transparent and reproducible reporting.

---

## 📂 Repository Structure

| File | Description |
|------|--------------|
| **1. Environmental_Raster_Preparation_Extraction.R** | Prepares and extracts relevant environmental variables for modeling. |
| **2. Complete_SDM_modeling_workflow.R** | Implements the full SDM workflow using MaxEnt, including model calibration and projection. |
| **3. Hypothesis_testing.R** | Tests the effect of temporal mismatches and model assumptions statistically. |
| **4. Independent_Validation_code.R** | Performs independent validation using external occurrence data. |
| **5. Validation_code_Overlap.R** | Evaluates prediction overlaps and consistency across different datasets. |
| **Species Data/** | Directory containing cleaned and formatted occurrence and environmental data used in the analyses. |

---

## 🧩 Dependencies

Core R packages used in the analysis include:  
`raster`, `terra`, `dplyr`, `readr`, `sf`, `ggplot2`, `maxnet`, and `ENMeval`.

---

## 🚀 Getting Started

1. Clone this repository to your local machine.  
2. Open the R project and install all listed dependencies.  
3. Run the scripts in numerical order for a complete analysis workflow.  
4. Replace input paths with your own data where necessary.  
5. Review outputs in the designated folders for metrics and raster maps.

---

## 🧠 Applications

- Evaluating the reliability of SDMs with temporal mismatches  
- Understanding data-timing effects on ecological modeling outcomes  
- Enhancing reproducibility and validation of pest distribution models under changing climates  

---

## 🗺️ License

This project is distributed under the MIT License, allowing reuse, modification, and distribution with appropriate attribution.

---

## 📜 Citation

> Juma, Q. (2025). *SDM Temporality Mismatch in Stemborer Distribution*. GitHub repository.  

---

## 💬 Contact

For inquiries or collaboration, please contact:  
**Quinto Juma**  
Email: [jmeltus@icipe.com]  or [meltusquinto@gmail.com]
Affiliation: International Centre of Insect Physiology and Ecology (icipe)
