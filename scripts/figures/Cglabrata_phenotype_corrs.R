## Purpose: Test relationships among phenotypes of interest 
## Author: Nancy Scott
## Email: scot0854@umn.edu
options(scipen = 999) 

source("../Candida_clinical_isolate_data/redcap_reports/MIC_data_summary.R")

# Load packages----
library(correlation)
library(patchwork)

# Get redcap data----
growth_curves <- '58045'

gc <- import_report(growth_curves) %>%
    filter(redcap_repeat_instrument != "NA", 
           !primary_id %in% c("AMS5122", "AMS5123", "AMS2401")) %>%
    select(primary_id, redcap_repeat_instance, 
           gc_date, drug_used, 
           gc_temp, gc_time, 
           k, r, t_gen, auc_l, auc_e) 

gc <- gc %>% 
    group_by(primary_id) %>% 
    filter(drug_used=="None", gc_date==max(gc_date), gc_time < 24.5) %>% 
    left_join(sample_info %>% select(primary_id,genus_species,series_id), 
              by = join_by(primary_id))

# Add MIC data ----
gc <- gc %>% 
    left_join(mic_info, by = join_by(primary_id, genus_species, series_id))

stationary_k <- mic_info %>% 
  group_by(primary_id) %>% 
  summarise(mean_k = mean(mean_no_drug_stationary_k),
            sd_k = sd(mean_no_drug_stationary_k))

# Subset to species of interest----
cglabrata_gc <- gc %>% 
  filter(genus_species=="C. glabrata")

# Scatterplots of potential correlations----
drugs <- as_labeller(c(fluconazole="Fluconazole", 
                       micafungin="Micafungin",
                       `amphotericin B` = "Amphotericin B"))


auc_stationary_k <- ggplot(cglabrata_gc, aes(y = auc_e, x=mean_no_drug_stationary_k)) +
  facet_wrap( ~ drug, labeller = drugs) +
  geom_point() + 
  theme_bw() +
  ylab("Area under the curve") +
  xlab("Stationary carrying capacy")

smg_stationary_k <- ggplot(cglabrata_gc, aes(y = mean_smg, x=mean_no_drug_stationary_k)) +
  facet_wrap( ~ drug, labeller = drugs) +
  geom_point() + 
  theme_bw() +
  ylab("Supra-MIC growth") +
  xlab("Stationary carrying capacy")

k_vs_stationary_k <- ggplot(cglabrata_gc, aes(y = k, x=mean_no_drug_stationary_k)) +
  facet_wrap( ~ drug, labeller = drugs) +
  geom_point() + 
  theme_bw() +
  ylab("Carrying capacity with shaking") +
  xlab("Stationary carrying capacy")

# Correlation testing----
cglabrata_drug_corrs <- cglabrata_gc %>% 
  group_by(drug) %>% 
  select(primary_id, drug, k, t_gen, auc_e, mean_smg, mean_no_drug_stationary_k) %>% 
  correlation(method = "pearson")
