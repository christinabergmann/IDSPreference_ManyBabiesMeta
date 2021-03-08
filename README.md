
# Data prep

Some parts of code need to be run multiple times to create different datasets for sensitivity analyses:

* To create the IPD age-matched replication dataset, we ran `2_get_tidy_MB_data.R`, `3_merge_MA_MB_data.R`, and `4_prep_for_analysis.R` sequentially with the parameter `age.matched = TRUE`.

* To create the replication dataset with more stringent inclusion criteria, we similarly ran `2_get_tidy_MB_data.R`, `3_merge_MA_MB_data.R`, and `4_prep_for_analysis.R` sequentially with the parameter `ic.dataset = TRUE`.

* After these data prep steps, the analysis script only needs to be run once. It conducts the main analyses as well as the sensitivity analyses. 

# Where is codebook?