# SunflowerAssociationMappingFieldSaltStress

**Overview**
The study investigates the genetic tolerance of sunflower (Helianthus annuus) to high salinity soils. It involves genome-wide association studies using the Sunflower Association Mapping panel grown in naturally occurring saline soils. The analysis covers six phenotypes: days to bloom, height, leaf area, leaf mass, oil percentage, and yield.

**Contents**
R Scripts:
  1. BLUP_ANOVA_example_OilPercentage.R: Performs BLUP (Best Linear Unbiased Prediction) analysis for oil percentage.
  2. Kriging.R: Implements kriging to estimate soil salinity for individual plots.
  3. outlier_removal.R: Removes outliers from the phenotypic data.

Input CSV Files:
  1. Phenotypes_with_Salinity.csv: Contains phenotypic data with salinity measurements.
  2. oil_noOutlier.csv: Oil percentage data with outliers removed.
  3. krig_train_data.csv: Training data for kriging analysis.
  4. krig_test_data.csv: Test data for kriging analysis.

**Usage**
  1. Ensure you have R and the required packages installed (dplyr, ggplot2, scales, magrittr, gstat, sp, asreml, ASRgenomics, ASRgwas).
  2. Clone this repository:
       `clone https://github.com/your-username/sunflower-salinity-tolerance.git`
  3. Set your working directory to the cloned repository.
  4. Run the R scripts in the following order:
       - First, run outlier_removal.R to clean the data.
       - Then, run Kriging.R to estimate soil salinity for individual plots.
       - Finally, run BLUP_ANOVA_example_OilPercentage.R for the BLUP analysis of oil percentage.

**Key Analyses**
  1. Outlier Removal: The outlier_removal.R script removes outliers from various phenotypic traits using influence plots and studentized residuals.
  2. Kriging: The Kriging.R script performs kriging to estimate soil salinity for individual plots, using data from two years (2016 and 2017) and two soil depths (0-6 cm and 6-24 cm).
  3. BLUP Analysis: The BLUP_ANOVA_example_OilPercentage.R script demonstrates the BLUP analysis for oil percentage, including model fitting and heritability calculation.

**Results**
The scripts generate various outputs including:
  - Cleaned phenotypic data
  - Estimated soil salinity for individual plots
  - BLUP values for oil percentage
  - Variance components and heritability estimates
These results were used in the genome-wide association study to identify loci associated with salinity tolerance in sunflower.

**Citation**
If you use this code or data in your research, please cite:
[Insert paper citation here when available]
License
[Insert chosen license here]
Contact
For any queries regarding the code or data, please contact:
[Insert contact information here]
