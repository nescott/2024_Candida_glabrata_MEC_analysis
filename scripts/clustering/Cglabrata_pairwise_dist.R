## ---------------------------
## Purpose: Calculate and plot pairwise differing sites or alleles
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
#options(scipen = 999)
## ---------------------------
## load packages
library(vcfR)
library(poppr)
library(tidyverse)
library(readxl)
library(rstatix)
library(scales)
library(patchwork)

source("ch2/redcap_MIC_summary.R")

species <- "C. glabrata"
in_variant_file <- "~/umn/data/variants/Cglabrata/Cglabrata_MEC_bwa_filtered_annotated.vcf.gz"
#in_variant_file <- "~/umn/data/variants/Cglabrata/Cglabrata_series_unfixed_snps.vcf.gz"
#series_ids <- scan("~/umn/data/metadata/Cglabrata_series_samples.txt", what = character())

in_patient_data <- "~/umn/data/metadata/2024_Cglabrata_sorted_patient_info.xlsx"
ordered_patient_data <-  "~/umn/data/metadata/Cglabrata_MEC_raxml_midpoint_tips.csv"
save_dir <- "~/umn/images/Cglabrata/"

species_ploidy <- 1

# Clade colors
color_file <- read.table("ch3/data/clade_colors.txt")
clade_colors <- color_file[,2]
names(clade_colors) <- color_file[,1]

## create genlight object from vcf file (can be .gz) and separate sample population file
vcf <- read.vcfR(in_variant_file)

pop.data <- read_excel(in_patient_data)
pop.data[pop.data=="NA"] <- NA

pt_order <- read.csv(ordered_patient_data, header = TRUE)
pt_order <- pt_order %>%
  left_join(pop.data %>% filter(study=="MEC") %>% select(sample, mec_pt_code, mec_isolate_code, ST))

#pt_order <- pt_order %>% filter(sample %in% series_ids)

all(colnames(vcf@gt)[-1] %in% pop.data$sample[1:98])

gl_species <- vcfR2genlight(vcf)
ploidy(gl_species) <- species_ploidy

# Get distance matrix
# setting differences_only to TRUE matches plink 1.9 --genome full output (IBS0 + IBS1 per row)
# differences_only FALSE counts differences in alleles, not just genotypes
intra_species_false <- bitwise.dist(gl_species, percent = FALSE, mat = TRUE, differences_only = FALSE)
intra_species_true <- bitwise.dist(gl_species, percent = FALSE, mat = TRUE, differences_only = TRUE)

# Pivot longer, add names
intra_species_false <- intra_species_false[pt_order$sample, pt_order$sample]

snp_dist <- intra_species_false %>%
  pull_upper_triangle(diagonal = TRUE) %>% 
  as_tibble() %>%
  rename(Var1 = rowname) %>% 
  mutate(Var1 = rownames(intra_species_false))
 
snp_dist <- snp_dist %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>% 
  filter(value !="") %>% 
  mutate(value = as.integer(value))

snp_dist <- snp_dist %>% 
  left_join(pt_order %>% select(sample, mec_pt_code, ST), by=c("Var1" = "sample")) %>%
  rename(pt_code_1 = mec_pt_code)

snp_dist <- snp_dist %>% 
  left_join(pt_order %>% select(sample, mec_pt_code, ST), by=c("Var2" = "sample")) %>%
  rename(pt_code_2 = mec_pt_code)

# Get patient info for adding days
sample_info <- sample_info %>% 
  mutate(primary_id = case_when(!is.na(secondary_id) ~ secondary_id, .default = primary_id))

# Within-strain vals and summaries
intra_strain_dist <- snp_dist %>% 
  filter(pt_code_1==pt_code_2, Var1!=Var2) %>% 
  distinct(value,pt_code_1,pt_code_2, .keep_all = TRUE) %>% 
  left_join(sample_info %>% select(primary_id, relative_days), by=c("Var1" = "primary_id")) %>% 
  rename(relative_days1 = relative_days) %>% 
  left_join(sample_info %>% select(primary_id, relative_days), by=c("Var2" = "primary_id")) %>% 
  rename(relative_days2 = relative_days)
  
intra_strain_dist <- intra_strain_dist %>% 
  mutate(days_apart = abs(relative_days1 - relative_days2)) %>% 
  arrange(pt_code_1, value, days_apart)

intra_strain_summary <- intra_strain_dist %>%
  group_by(pt_code_1) %>%
  summarize(median_snp = round(median(value), digits = 0), max_snp = max(value), min_snp = min(value))

# Series dist estimate from manual curation
seriesbb <- read.csv("~/umn/data/variants/Cglabrata/Cglabrata_MEC_BB_SNPS.csv", header = T, row.names = 2)
seriesbb <- seriesbb[,2:9]
bbtest <- data.table::transpose(seriesbb)
colnames(bbtest) <- rownames(seriesbb)
rownames(bbest) <- colnames(seriesbb)
dist(as.matrix(bbtest), method = "manhattan")

# Isolate collection ranges mean don't bother reporting this
snp_time_corr <- correlation::correlation(intra_strain_dist, select = c("value", "days_apart"))
ggplot(intra_strain_dist, aes(y = value, x = days_apart, colour = pt_code_1)) + 
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = paletteer_d("ggthemes::Tableau_20"))

# Between strain vals and summary
inter_strain_dist <- snp_dist %>% 
  filter(pt_code_1!=pt_code_2) %>% 
  arrange(value)

intra_st_dist <- inter_strain_dist %>% 
  filter(ST.x == ST.y) %>% 
  arrange(value)

between_st_dist <- inter_strain_dist %>% 
  filter(ST.x != ST.y) %>% 
  arrange(value)

between_st_summary <- between_st_dist %>% 
  group_by(ST.x,ST.y) %>% 
  summarise(mean_snp_count = mean(value),
            max_snp_count = max(value),
            min_snp_count = min(value))

