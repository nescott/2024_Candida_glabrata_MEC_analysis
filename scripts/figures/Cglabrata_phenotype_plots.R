## ---------------------------
## Purpose: Summary plots of growth curve metrics
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999) 

source("~/umn/2024_Candida_clinical_isolate_phenotyping/redcap_reports/MIC_data_summary.R")
library(patchwork)
library(ggbeeswarm)

drugs <- as_labeller(c(fluconazole="Fluconazole", 
                       micafungin="Micafungin",
                       `amphotericin B` = "Amphotericin B"))

# MIC sub-plots
flc <- ggplot(mic_info %>% filter(genus_species == "C. glabrata", drug=="fluconazole"), aes(x=mic50))+
    geom_bar(fill = "#999933", just = 1) +
    scale_x_discrete(limits = c("0.5", "1", "2", "4", "8", "16", "32", ">32")) +
    scale_y_continuous(limits = c(0,100)) +
    theme_bw() +
    theme(axis.text.x = element_text(hjust = 1.4, size = 8)) +
    theme(axis.title = element_text(size = 13)) +
    geom_vline(data = filter(mic_info, genus_species=="C. glabrata"),
             aes(xintercept = "16"), linetype = 2) + 
    xlab("\nFluconazole MIC50") +
    ylab(NULL)

mcf <- ggplot(mic_info %>% filter(genus_species == "C. glabrata", drug=="micafungin"), aes(x=mic50))+
    geom_bar(fill = "#657b1f", just = 1) +
    scale_x_discrete(limits = c("0.016", "0.032", "0.064", "0.125", "0.256", "0.5", "1", ">1"),
                     labels = c ("0.01", "0.03", "0.06", "0.125", "0.25 ", "0.5 ", "1  ", ">1 ")) +
    scale_y_continuous(limits = c(0,100)) +
    theme_bw() +
    theme(axis.text.y = element_blank(), axis.ticks = element_blank(),) +
    theme(axis.text.x = element_text(hjust = 1, size = 8)) +
    theme(axis.title = element_text(size = 13)) +
    geom_vline(data = filter(mic_info, genus_species=="C. glabrata"),
             aes(xintercept = "0.032"), linetype = 2) +
    xlab("\nMicafungin MIC50") +
    ylab(NULL)

amb <- ggplot(mic_info %>% filter(genus_species == "C. glabrata", drug=="amphotericin B"), aes(x=mic50))+
    geom_bar(fill = "#305d0b", just = 1) +
    scale_x_discrete(limits = c("0.016", "0.032", "0.064", "0.125", "0.256", "0.5", "1", ">1"),
                     labels = c ("0.01", "0.03", "0.06", "0.125", "0.25 ", "0.5 ", "1  ", ">1  ")) +
    scale_y_continuous(limits = c(0,100)) +
    theme_bw() +
    theme(axis.text.y = element_blank(), axis.ticks = element_blank(),) +
    theme(axis.text.x = element_text(hjust = 1, size = 8)) +
    theme(axis.title = element_text(size = 13)) +
  geom_vline(data=filter(mic_info, genus_species=="C. glabrata"),
             aes(xintercept = "1"), linetype =2) +
    xlab("\nAmphotericin B MIC90") +
    ylab(NULL)

flc + mcf + amb
ggsave("images/Cglabrata/MIC_GC_SMG/Cglabrata_MEC_MICs.tiff", bg="white", width = 8, height = 2.93, units = "in", device = tiff, dpi=300)

# SMG subplots
flc_smg <- ggplot(mic_info%>% filter(genus_species == "C. glabrata", drug == "fluconazole"), 
                  aes(x=genus_species, y=mean_smg)) + 
    geom_beeswarm(size=1, cex = 2, color = "#999933") +
    theme_bw() +
    scale_y_continuous(limits = c(0,0.8)) +
    theme(axis.title = element_text(size = 13)) +
    ylab("Supra-MIC growth") +
    xlab("Fluconazole") +
    theme(axis.text.x = element_blank())

mcf_smg <- ggplot(mic_info%>% filter(genus_species == "C. glabrata", drug=="micafungin"), 
                  aes(x=genus_species, y=mean_smg)) + 
  geom_beeswarm(size=1, cex = 2, color = "#657b1f") +
  theme_bw() +
  scale_y_continuous(limits = c(0,0.8)) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(),) +
  theme(axis.title = element_text(size = 13)) +
  ylab(NULL) +
  xlab("Micafungin") +
  theme(axis.text.x = element_blank())

amb_smg <- ggplot(mic_info%>% filter(genus_species == "C. glabrata", drug=="amphotericin B"), 
                  aes(x=genus_species, y=mean_smg)) + 
  geom_beeswarm(size=1, cex = 2, color = "#305d0b") +
  theme_bw() +
  scale_y_continuous(limits = c(0,0.8)) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(),) +
  theme(axis.title = element_text(size = 13)) +
  ylab(NULL) +
  xlab("Amphotericin B") +
  theme(axis.text.x = element_blank())

flc_smg + mcf_smg + amb_smg

ggsave("Cglabrata_MEC_SMGs.tiff", bg="white", width = 11, height = 4, units = "in", device = tiff, dpi=300)
