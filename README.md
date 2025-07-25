### 01_format_capture_histories.R

* Start with Motus detections, filter detections depending on known issues with data (ex. poor quality stations, manual filtering)
* Produce capture histories based on detections over each of 300 days, then bin into 10 day capture occasions

### 02_estimate_phenotypes_migration.R

* Use Motus detections to estimate bearing & timing at a range of distances from release site & Pemberton (centre of hybrid zone)
* Only use first bearings for ML classification because later time points don't have enough data

### 03_get_survival_estimates.R

* Run CJS models on capture histories to get a binary estimate of survival

### fig1_map_format.R

* Map Motus station detections for all birds
* Format Fig 1 - map + pheno network

### fig1_phenotype_network.R

* Make a correlation network & run bootstrapping to remove non-significant correlations
* Use this for Fig 1

### fig2_shap_etc.Rmd

* Make Fig 2
* Model evaluation fig - confusion matrix, SHAP values, feature importances