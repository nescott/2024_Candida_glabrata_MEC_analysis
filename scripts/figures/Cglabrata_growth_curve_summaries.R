## Purpose: Summary plots of growth curve metrics
## Author: Nancy Scott
## Email: scot0854@umn.edu
options(scipen = 999) 

library(patchwork)
library(ggbeeswarm)

# Variables----
save_path <- "~/umn/images/2023_growth_curves/"

growth_curves <- '58045'

species_colors <- c("#88CCEE", "#999933", "#CC6677", "#44AA99", "#117733", "#332288",
                    "#882255","#BBBBBB", "#AA4499", "#DDCC77", "black")

# Report import function----
import_report <- function(report_number) {
  url <- "https://redcap.ahc.umn.edu/redcap/api/"
  formData <- list("token"=Sys.getenv("redcap_api_key"),
                   content='report',
                   format='csv',
                   report_id=report_number,
                   csvDelimiter='',
                   rawOrLabel='label',
                   rawOrLabelHeaders='raw',
                   exportCheckboxLabel='true',
                   returnFormat='csv'
  )
  response <- httr::POST(url, body = formData, encode = "form")
  result <- httr::content(response, show_col_types = FALSE)
}

# Read in and join reports----
source("../Candida_clinical_isolate_data/redcap_reports/MIC_data_summary.R")

gc <- import_report(growth_curves) %>%
    filter(redcap_repeat_instrument != "NA", 
           !primary_id %in% c("AMS5122", "AMS5123")) %>%
    select(primary_id, redcap_repeat_instance, 
           gc_date, drug_used, 
           gc_temp, gc_time, 
           k, r, t_gen, auc_l, auc_e) 

gc <- gc %>% 
    group_by(primary_id) %>% 
    filter(drug_used=="None", gc_date==max(gc_date), gc_time < 24.1)  

gc <- gc %>% 
  left_join(mic_info, by=join_by(primary_id))

# for ordering species and colors for consistency----
species_count <- sample_info %>%
    group_by(genus_species) %>%
    summarize(species_count=n()) %>%
    arrange(desc(species_count))

species_colors <- species_colors %>% 
    set_names(species_count$genus_species)

gc$genus_species <- factor(gc$genus_species, levels = species_count$genus_species)

# Plot AUC per drug----
flc_auc <- gc %>% 
    filter(genus_species %in% c("C. glabrata"), drug=="fluconazole") %>% 
    ggplot(aes(x=mic50, y=auc_e)) + 
    geom_beeswarm(color = "#999933") +
    theme_bw() +
    scale_y_continuous(limits = c(0,25))+
    xlab("Fluconazole MIC50") +
    ylab("AUC, 24H")

mcf_auc <- gc %>% 
  filter(genus_species %in% c("C. glabrata"), drug=="micafungin") %>% 
  ggplot(aes(x=mic50, y=auc_e)) + 
  geom_beeswarm(color = "#657b1f") +
  theme_bw() +
  scale_y_continuous(limits = c(0,25))+
  xlab("Micafungin MIC50") +
  ylab(NULL)

amb_auc <- gc %>% 
  filter(genus_species %in% c("C. glabrata"), drug=="amphotericin B") %>% 
  ggplot(aes(x=mic50, y=auc_e)) + 
  geom_beeswarm(color = "#305d0b") +
  theme_bw() +
  scale_y_continuous(limits = c(0,25))+
  xlab("Amphotericin B MIC90") +
  ylab(NULL)

# Merge and save----
combined_auc <- flc_auc + mcf_auc + amb_auc

ggsave("images/Cglabrata/MIC_GC_SMG/2023_AUC_by_drug_mic.png", 
       combined_auc,
       width = 11, 
       height = 4, 
       units = "in",
       device = png, 
       dpi=300, 
       bg="white")

