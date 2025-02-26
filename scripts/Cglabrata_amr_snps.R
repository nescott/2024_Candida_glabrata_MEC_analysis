## ---------------------------
## Purpose: Find SNVs private to drug-resistant isolates
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
source("~/umn/thesis/ch2/redcap_MIC_summary.R", echo=TRUE)

cglab_flc <- resistant_isolates %>% filter(genus_species=="C. glabrata", eucast_breakpoint_fluconazole=="R")
cglab_flc_mutations <- read_delim("~/umn/data/variants/Cglabrata/Cglabrata_MEC_FLC_candidates.txt", 
                                  col_types = "cdccccccccccccccccccccccccccccccccccc")
cglab_flc_mutations <- cglab_flc_mutations %>%   
  pivot_longer(cols= -c(CHROM:GENE,type:aa_change), names_to = "sample", values_to = "genotype") %>%
  filter(sample !="MEC208")
flc_candidates <- cglab_flc_mutations %>% group_by(ORF,POS) %>% filter(genotype=="1")

pdr1_check <- flc_candidates %>% 
  filter(GENE=="PDR1") 

flc_mutator_check <- flc_candidates %>% filter(GENE %in% c("MLH1", "MSH2"))

cglab_mcf <- resistant_isolates %>% filter(genus_species=="C. glabrata", eucast_breakpoint_micafungin=="R")
cglab_mcf_mutations <- read_delim("~/umn/data/variants/Cglabrata/Cglabrata_MEC_MCF_candidates.txt", col_types = "cdcccccccccccccccccccc")
cglab_mcf_mutations <- cglab_mcf_mutations %>% pivot_longer(cols= -c(CHROM:GENE,type:aa_change), names_to = "sample", values_to = "genotype")
mcf_candidates <- cglab_mcf_mutations %>% filter(genotype=="1")

fks_check <- mcf_candidates %>% filter(GENE %in% c("FKS1", "FKS2"))
mcf_mutator_check <- mcf_candidates %>% filter(GENE %in% c("MLH1","MSH2"))

