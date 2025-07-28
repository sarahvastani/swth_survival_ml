### 01_format_capture_histories.R

* Start with Motus detections, filter detections depending on known issues with data (ex. poor quality stations, manual filtering)
* Produce capture histories based on detections over each of 300 days, then bin into 10 day capture occasions

### 02_estimate_phenotypes_migration.R

* Use Motus detections to estimate bearing & timing at a range of distances from release site & Pemberton (centre of hybrid zone)
* Only use first bearings for ML classification because later time points don't have enough data

### 03_get_survival_estimates.R

* Run CJS models on capture histories to get a binary estimate of survival

### 04_ml_workflow_survival_prediction.ipynb

* Execute machine learning pipeline to predict Swainson’s thrush migratory survival with Random Forest, including preprocessing, SMOTE balancing, hyperparameter tuning, and SHAP interpretation
* Generate code and outputs for Figures S1a, S1b, S3, and S4

### fig1_map_format.R

* Map Motus station detections for all birds
* Format Fig 1 - map + pheno network

### fig1_phenotype_network.R

* Make a correlation network & run bootstrapping to remove non-significant correlations
* Use this for Fig 1

### fig2_shap_etc.Rmd

* Make Fig 2
* Model evaluation fig - confusion matrix, SHAP values, feature importances

### figS1c_S2_tables_drawings.pdf

* Make Fig S1c and S2
* Visualize SMOTE resampling and KNN imputation using tables and drawings

### figS5_Workflow.pdf
* Visually summarize the machine learning pipeline to predict Swainson’s thrush migratory survival
