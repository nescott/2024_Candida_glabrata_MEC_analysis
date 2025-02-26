## ---------------------------
## Purpose: Replot genome-view of LOH and CNV from saved excel file (see "genome_vis.R)
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# Input file variables
genome_df_file <- "~/umn/data/genome_plots/Cglabrata/2024-02-15_MEC335.xlsx"   #args[1]
sample_id <- "54-8"  #args[2] # or "YourID"
feature_file <- "~/umn/Candida_genome_visualization/ref_genome_files/Cglabrata_CGD_s05m03r02_features.txt"  #args[3] # or "path/to/features.txt"
label_file <- "~/umn/Candida_genome_visualization/ref_genome_files/Cglabrata_CGD_s05m03r02_chr_labels.txt"  #args[4] # or "path/to/chr_labels.txt"

# Load packages
library(readxl)
library(tidyverse)
library(ggplot2)

# Set variables
window <- 5000 # size of window used for rolling mean and snp density
ploidy <- 1

y_axis_labels <- c(1,2,3)  # manual y-axis labels, adjust as needed
inter_chr_spacing <- 100000 # size of space between chrs

save_dir <-  "" 
ref <- "CBS138" 

cnv_color <- "dodgerblue4"

feature_colors <- c("white", "grey26", "deepskyblue")
feature_shapes <- c(24,21,22)

ploidy_multiplier <- 3  # this multiplied by ploidy sets the max-y scale

chrom_outline_color <- "gray15"  # color of chromosome outlines

chrom_line_width <- 0.2  # line width of chromosome outlines

# X-axis labels overwrite input scaffold names in final plot
chr_ids <- scan(label_file, what = character())
chr_ids <- LETTERS[1:13]
# Dataframe of joined copy number, snps, and plotting positions per window
genome_depth <- read_xlsx(genome_df_file, sheet=1)

# Small dataframes for chrom. outlines and features
chroms <- genome_depth %>%
  group_by(index) %>%
  summarise(xmin=min(plot_pos), xmax=max(plot_pos), ymin=0, ymax=Inf)

features <- read_tsv(feature_file, show_col_types = FALSE)

features <- features %>%
     group_by(chr, index=consecutive_id(chr)) %>%
     left_join(chroms, by=join_by(index))

features <- features %>%
     mutate(plot_start = start + xmin, plot_end = end + xmin)

# Tick marks to center chromosome ID label
ticks <- tapply(genome_depth$plot_pos, genome_depth$index, quantile, probs =
                0.5, na.remove = TRUE)

# Plot linear genome
p <- ggplot(genome_depth) +
  geom_segment(aes(x = plot_pos,
                   y = ifelse(copy_number <= ploidy*ploidy_multiplier, copy_number, Inf),
                   xend = plot_pos, yend = ploidy), alpha = 0.9, color = cnv_color) +
  geom_rect(data=chroms, aes(group=index, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            linewidth = chrom_line_width, fill = NA,
            colour = chrom_outline_color, linejoin = "round", inherit.aes = FALSE) +
  geom_point(data = features, size = 1,
               aes(group=index, x=plot_start, y=ymin, shape = Feature, fill = Feature),
               position = position_nudge(y=0.07)) +
    scale_fill_manual(values = feature_colors, guide = "none") +
    scale_shape_manual(values = feature_shapes, guide = "none") +
  ylab(sample_id) +
  scale_x_continuous(name = NULL, expand = c(0, 0), breaks = ticks, labels=chr_ids, position = "top") +
  scale_y_continuous(limits = c(0, ploidy*ploidy_multiplier), breaks = y_axis_labels) +
  theme_classic() +
  theme(
        axis.ticks = element_line(color = NA),
        axis.line = element_blank(),
        axis.title.y = element_text(size = 12, color = "black", 
                                    margin = margin(r=5, l=0),
                                    angle = 0,
                                    vjust= 0.6),
        axis.text.y = element_text(size = 10, color = "black"),
        #legend.position = "bottom"
        axis.text.x = element_text(size=10, color ="black")
        )

# Save plot
ggsave(sprintf("%s%s_%s_%s_%sbp.pdf",
               save_dir,
               Sys.Date(),
               "CCNC",
               ref,
               window),
       p2/p,
       width = 6,
       height = 2,
       units = "in",
       #device = png,
       #dpi = 300,
       bg = "white")

################################################################################

small_cnv <- genome_depth %>% 
  filter(index==1)

small_chroms <- small_cnv %>%
  group_by(index) %>%
  summarise(xmin=min(plot_pos), xmax=max(plot_pos), ymin=0, ymax=Inf)

small_features <- features <- read_tsv(feature_file, show_col_types = FALSE)

small_features <- small_features %>%
  group_by(chr, index=consecutive_id(chr)) %>%
  filter(index %in% c(1)) %>% 
  left_join(small_chroms, by=join_by(index)) %>%
  mutate(plot_start = start + xmin, plot_end = end + xmin)

small_ticks <- tapply(small_cnv$plot_pos, small_cnv$index, quantile, probs =
                        0.5, na.remove = TRUE)
mec335 <- 
 ggplot(small_cnv) +
  geom_segment(aes(x = plot_pos,y = ifelse(copy_number <= ploidy*ploidy_multiplier, copy_number, Inf),
                   xend = plot_pos, yend = ploidy), alpha = 0.9, color = cnv_color,
               show.legend = FALSE) +
  geom_point(data = small_features, size = 1.3,
             aes(group=index, x=plot_start, y=ymin, shape = Feature, fill = Feature),
             position = position_nudge(y=0.1)) +
  scale_fill_manual(values = c("white", "grey26", "deepskyblue"), guide = "none") + 
  scale_shape_manual(values = c(24,21,22), guide = "none") +
  #ylab(sample_id) +
  geom_rect(data=small_chroms, aes(group=index, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            linewidth = chrom_line_width, fill = NA,
            colour = chrom_outline_color, linejoin = "round", inherit.aes = FALSE) +
  scale_x_continuous( name = NULL, expand = c(0, 0), breaks = small_ticks, labels=paste(sample_id, "Chr A"), position = "top") +
  scale_y_continuous(limits = c(0, ploidy*ploidy_multiplier), breaks = y_axis_labels) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(color = "grey90"),
        axis.ticks = element_line(color = NA),
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        #axis.title.y = element_text(size = 12, color = "black", 
         #                           margin = margin(r=5, l=0),
          #                          angle = 0,
           #                         vjust= 0.6),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size=11, color ="black"))


case54 <- (mec328 + mec329)/
  (mec330 + mec331)/
  (mec332 + mec333)/
  (mec334 + mec335)

ggsave(sprintf("%s%s_%s_%s_%sbp.pdf",
               save_dir,
               Sys.Date(),
               "case54",
               ref,
               window),
       case54,
       width = 4,
       height = 4,
       units = "in",
       #device = png,
       #dpi = 300,
       bg = "white")
