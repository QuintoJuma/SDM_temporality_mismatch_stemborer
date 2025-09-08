############################################################
# Title: Model Performance Comparison between Scenarios
# Author: Quinto Juma Meltus
# Date: 2025-09-01
# Description:
#   This script compares internal and independent (external) 
#   model performance metrics between two scenarios. It 
#   summarizes evaluation metrics, conducts paired t-tests,
#   interprets results, and generates comparative plots.
#
# Notes:
#   - Replace "path/to/..." with the actual location of your data files.
#   - Ensure the CSVs follow the expected structure (see below).
#
# Required Input Files:
#   1. scenario_internal_model_testing_metrics.csv
#        Columns: Scenario, Averaged_AUC, Averaged_TSS
#   2. independent_validation_testing.csv
#        Columns: Scenario, Precision_Accuracy, F1_Score
#
# Required R packages:
#   dplyr, ggplot2, tidyr
#
############################################################

# ----------------------------
# Load required libraries
# ----------------------------
library(dplyr)
library(ggplot2)
library(tidyr)

# ----------------------------
# Load data
# ----------------------------
# Replace file paths with your own
internal <- read.csv("path/to/scenario_internal_model_testing_metrics.csv")
independent <- read.csv("path/to/independent_validation_testing.csv")

# ----------------------------
# Internal model metrics
# ----------------------------
# Summarize mean AUC and TSS by scenario
internal_summary <- internal %>%
  group_by(Scenario) %>%
  summarise(
    mean_AUC = mean(Averaged_AUC, na.rm = TRUE),
    mean_TSS = mean(Averaged_TSS, na.rm = TRUE)
  )

print("Internal Model Metrics Summary:")
print(internal_summary)

# Hypothesis testing for internal metrics
auc_test_internal <- t.test(
  internal$Averaged_AUC[internal$Scenario == 1],
  internal$Averaged_AUC[internal$Scenario == 2],
  paired = TRUE,
  alternative = "greater" # H1: Scenario 1 > Scenario 2
)

tss_test_internal <- t.test(
  internal$Averaged_TSS[internal$Scenario == 1],
  internal$Averaged_TSS[internal$Scenario == 2],
  paired = TRUE,
  alternative = "greater"
)

print("Internal Model Metrics Hypothesis Tests:")
print(auc_test_internal)
print(tss_test_internal)

# ----------------------------
# Independent validation metrics
# ----------------------------
# Summarize mean precision accuracy and F1 score by scenario
independent_summary <- independent %>%
  group_by(Scenario) %>%
  summarise(
    mean_Precision_Accuracy = mean(Precision_Accuracy, na.rm = TRUE),
    mean_F1 = mean(F1_Score, na.rm = TRUE)
  )

print("Independent Validation Metrics Summary:")
print(independent_summary)

# Hypothesis testing for independent validation metrics
pa_test_independent <- t.test(
  independent$Precision_Accuracy[independent$Scenario == 1],
  independent$Precision_Accuracy[independent$Scenario == 2],
  paired = TRUE,
  alternative = "greater"
)

f1_test_independent <- t.test(
  independent$F1_Score[independent$Scenario == 1],
  independent$F1_Score[independent$Scenario == 2],
  paired = TRUE,
  alternative = "greater"
)

print("Independent Validation Hypothesis Tests:")
print(pa_test_independent)
print(f1_test_independent)

# ----------------------------
# Interpretation Section
# ----------------------------
cat("\n--- Interpretation ---\n")

if (auc_test_internal$p.value < 0.05 | tss_test_internal$p.value < 0.05) {
  cat("Internal metrics show a significant difference in at least one metric between Scenario 1 and Scenario 2.\n")
} else {
  cat("Internal metrics do not show a significant difference between Scenario 1 and Scenario 2.\n")
}

if (pa_test_independent$p.value < 0.05 | f1_test_independent$p.value < 0.05) {
  cat("Independent validation metrics indicate that Scenario 1 performs significantly better than Scenario 2.\n")
} else {
  cat("Independent validation metrics do not indicate a significant difference.\n")
}

cat("\nHypothesis:\n")
cat("H0: There is no significant difference in model performance between Scenario 1 and Scenario 2.\n")
cat("H1: Scenario 1 yields significantly better model performance than Scenario 2.\n")

# ----------------------------
# Visualization
# ----------------------------
# Prepare data for plotting
plot_data <- internal_summary %>%
  mutate(Source = "Internal") %>%
  rename(Metric1 = mean_AUC, Metric2 = mean_TSS) %>%
  bind_rows(
    independent_summary %>%
      mutate(Source = "External") %>%
      rename(Metric1 = mean_Precision_Accuracy, Metric2 = mean_F1)
  ) %>%
  tidyr::pivot_longer(
    cols = c(Metric1, Metric2),
    names_to = "Metric",
    values_to = "Value"
  )

# Rename metrics for clarity
plot_data$Metric <- factor(
  plot_data$Metric,
  levels = c("Metric1", "Metric2"),
  labels = c("AUC / Precision Accuracy", "TSS / F1 Score")
)

# Plot grouped bar chart
ggplot(plot_data, aes(x = Metric, y = Value, fill = factor(Scenario))) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  facet_wrap(~Source) +
  scale_fill_manual(
    values = c("#1b9e77", "#d95f02"),
    name = "Scenario",
    labels = c("Scenario 1", "Scenario 2")
  ) +
  labs(
    x = "Metric", y = "Mean Value",
    title = "Scenario Performance Comparison",
    subtitle = "Internal vs External Validation Metrics"
  ) +
  theme_minimal(base_size = 13)

############################################################
# End of Script
############################################################
