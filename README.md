
# Datasets

### Key points
* Meta-analysis codebook: https://github.com/christinabergmann/IDSPreference_ManyBabiesMeta/blob/master/data/codebook.xlsx

* Main-analysis dataset: `mb_ma_combined_prepped_0.125.csv`

* Dataset used for sensitivity analyses with more stringent inclusion criterion: `mb_ma_combined_prepped_0.75.csv`

* Dataset used for age-matched sensitivity analyses: `mb_ma_combined_prepped_0.125_age_matched`


### Fine print
* Fractional suffixes on dataset names (e.g., "0.125") correspond to inclusion criteria. The "0.75" dataset (i.e., with the more stringent subject inclusion criterion) has fewer estimates than the main dataset (equivalent to the "0.125" dataset) because some age groups were dropped completely if N<10. 




# Code

### Data prep code

To prep the main-analysis dataset and the various datasets used for sensitivity analyses, simply run `0_master_prep.R`. This file sets different combinations of global variables regarding the inclusion criteria and whether the dataset should be age-matched and then calls a sequence of other data-prep files (`1_get_tidy_MA_data.R`, `2_get_tidy_MB_data.R`, `3_merge_MA_MB_data.R`, and `4_prep_for_analysis.R`).  



### Analysis code

The analysis script conducts the main analyses as well as the sensitivity analyses. It uses only the datasets containing "combined_prepped" in the names. 


# Results files

* In `tables_to_prettify`, titles with uppercase model names represent main analyses
