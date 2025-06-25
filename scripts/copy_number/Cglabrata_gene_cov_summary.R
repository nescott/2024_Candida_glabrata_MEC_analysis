## Purpose: Add coverage to spreadsheet of del/dup ORFs
## Author: Nancy Scott
## Email: scot0854@umn.edu

library(tidyverse)
library(readxl)
library(writexl)

# Deletion/truncations
deletion_summary <- read.delim("ch4/Cglabrata_MEC_candidate_deletions.tsv")

del_cov <- read.delim("ch4/Cglabrata_MEC_candidate_deletion_coverage.tsv") 

del_merged <- deletion_summary %>% 
  left_join(del_cov, by = c("ID"="orf", "sample")) %>% 
  select(-c(score,phase,Name,parent_feature_type)) %>% 
  arrange(ID, sample)

write_xlsx(del_merged, "ch4/Cglabrata_MEC_deletion_cov.xlsx")

# Duplications
dup_summary <- read.delim("ch4/Cglabrata_MEC_candidate_duplications.tsv")
                    
dup_cov <- read.delim("ch4/Cglabrata_MEC_candidate_duplication_coverage.tsv")

dup_merged <- dup_summary %>% 
  left_join(dup_cov, by = c("ID"="orf", "sample")) %>% 
  select(-c(score,phase,Name,parent_feature_type)) %>% 
  arrange(ID, sample)

write_xlsx(dup_merged, "ch4/Cglabrata_MEC_duplication_cov.xlsx")