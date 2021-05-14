
# Datasets

* `from_data_team` has the Dunst datasets that Maya received from the data-collection team. The corresponding MB datasets are on [Git]("https://raw.githubusercontent.com/manybabies/mb1-analysis-public/master/processed_data/03_data_diff_main.csv") and are automatically pulled in by the script `2_get_tidy_MB_data.R`.



### Key points
* Meta-analysis codebook: https://github.com/christinabergmann/IDSPreference_ManyBabiesMeta/blob/master/data/codebook.xlsx

* Main-analysis dataset: `mb_ma_combined_prepped_0.125.csv`

* Dataset used for sensitivity analyses with more stringent inclusion criterion: `mb_ma_combined_prepped_0.75.csv`

* Dataset used for age-matched sensitivity analyses: `mb_ma_combined_prepped_0.125_age_matched`


### Fine print
* Fractional suffixes on dataset names (e.g., "0.125") correspond to inclusion criteria. The "0.75" dataset (i.e., with the more stringent subject inclusion criterion) has fewer estimates than the main dataset (equivalent to the "0.125" dataset, corresponding to one of 8 trials per condition, the same as the published analyses) because some age groups were dropped completely if N<10. 



# Code

### Data prep code

To prep the main-analysis dataset and the various datasets used for sensitivity analyses, simply run `0_master_prep.R`. This file sets different combinations of global variables regarding the inclusion criteria and whether the dataset should be age-matched and then calls a sequence of other data-prep files (`1_get_tidy_MA_data.R`, `2_get_tidy_MB_data.R`, `3_merge_MA_MB_data.R`, and `4_prep_for_analysis.R`).  

The analyses with the original vs. corrected Dunst dataset are set up as a parallel directory structure in both the data-prep and analysis directories. The reproducible manuscript on Overleaf has a subdirectory called `corrected_dunst` for those results. 



### Analysis code

The analysis script conducts the main analyses as well as the sensitivity analyses. It uses only the datasets containing "combined_prepped" in the names. 


# Results files

* In `tables_to_prettify`, titles with uppercase model names represent main analyses
