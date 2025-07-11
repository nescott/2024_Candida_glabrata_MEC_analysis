## Purpose: Plot genome-wide view of rolling mean read depth
## Email: scot0854@umn.edu
## Notes: Adapted from https://github.com/berman-lab/ymap
## and https://github.com/stajichlab/C_lusitaniae_popseq.
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# Input file variables
read_depth_file <-  args[1]
sample_id <- args[2]
feature_file <- args[3]
label_file <- args[4]
mito <- "mito_C_glabrata_CBS138"

# Load packages----
library(RcppRoll)
library(tidyverse)
library(ggplot2)
library(writexl)
library(readxl)

# Variables----
window <- 5000 # size of window for rolling mean and snp density
ploidy <- 1

y_axis_labels <- c(1,2,3)
inter_chr_spacing <- 150000

image_dir <- "~/umn/images/Cglabrata/2024-02-15_genome_plots/"
ref <- "cgd_s05m03r02"
data_dir <-"~/umn/data/genome_plots/Cglabrata/"

chr_ids <- scan(label_file, what = character())

snp_low <- "white"  # 2-color gradient scale
snp_high <- "black"  # 2-color gradient scale
cnv_color <- "dodgerblue4"  # copy number color

ploidy_multiplier <- 3  # this multiplied by ploidy sets the max-y scale

chrom_outline_color <- "gray15"  # color of chromosome outlines

chrom_line_width <- 0.2  # line width of chromosome outlines

# Base R doesn't have a mode calculation
Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

# Relative copy number calcs from samtools depth input----
genome_raw <- read.table(read_depth_file, header = TRUE)
genome_raw <- genome_raw %>%
  filter(chr != mito)
genome_raw$rolling_mean <- roll_mean(genome_raw$depth, window)[seq_len(length(genome_raw$chr))]

raw_genome_median <- median(genome_raw$depth) # includes all chromosomes, may need correcting

chr_median <- genome_raw %>%  # checking each chromosome for outliers relative to genome
  group_by(chr) %>%
  summarise(chr_mode = Modes(depth), chr_med = median(depth))  # can manually compare mode and median if questioning median

subset_chr_median <- chr_median %>%  # modest filtering to avoid aneuploidy skew of "normal" genome depth
  filter(chr_med <= raw_genome_median *1.15 & chr_med >= raw_genome_median * 0.85)

genome_median <- median(subset_chr_median$chr_med)  # filtered median used to calculate relative depth

# Reshape dataframe for plotting----
genome_window <- genome_raw %>%
  group_by(chr, index=consecutive_id(chr)) %>%
  reframe(position=unique((pos %/% window)*window+1))

genome_window <- genome_window %>%
  group_by(index) %>%
  arrange(index) %>%
  mutate(chr_length = max(position))

chrs <- as.vector(unique(genome_window$chr_length))
chr_plot <- c()
for(i in 1:length(chrs)){chr_plot[i] <- sum(chrs[1:i-1])}

# Calculate relative depth and copy number----
genome_depth <- genome_raw %>%
  group_by(chr, index=consecutive_id(chr)) %>%
  filter(pos %in% genome_window$position) %>%
  mutate(relative_depth = rolling_mean/genome_median) %>%
  mutate(copy_number= relative_depth * ploidy) %>%
  mutate(chr_sums=chr_plot[index]) %>%  # for proper x-axis plotting
  mutate(plot_pos=ifelse(index==1, pos, (pos+chr_sums+(inter_chr_spacing*(index-1))))) # for proper x-axis plotting

# Small dataframes for chrom. outlines and features----
chroms <- genome_depth %>%
  group_by(index) %>%
  summarise(xmin=min(plot_pos), xmax=max(plot_pos), ymin=0, ymax=Inf)

features <- read_tsv(feature_file, show_col_types = FALSE)

features <- features %>%
     group_by(chr, index=consecutive_id(chr)) %>%
     left_join(chroms, by=join_by(index))

features <- features %>%
     mutate(plot_start = start + xmin, plot_end = end + xmin)

# Tick marks to center chromosome ID label----
ticks <- tapply(genome_depth$plot_pos, genome_depth$index, quantile, probs =
                0.5, na.remove = TRUE)

# Plot linear genome----
p <- ggplot(genome_depth) +
  geom_segment(aes(x = plot_pos,
                   y = ifelse(copy_number <= ploidy*ploidy_multiplier, copy_number, Inf),
                   xend = plot_pos, yend = ploidy), alpha = 0.9, color = cnv_color) +
  geom_rect(data=chroms, aes(group=index, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            linewidth = chrom_line_width, fill = NA,
            colour = chrom_outline_color, linejoin = "round", inherit.aes = FALSE) +
  geom_point(data = features, size = 2,
               aes(group=index, x=plot_start, y=ymin, shape = Feature, fill = Feature),
               position = position_nudge(y=0.07)) +
  scale_fill_manual(values = c("white", "grey26", "deepskyblue")) +
  scale_shape_manual(values = c(24,21,22)) +
  ylab(sample_id) +
  scale_x_continuous(name = NULL, expand = c(0, 0), breaks = ticks, labels=chr_ids) +
  scale_y_continuous(limits = c(0, ploidy*ploidy_multiplier), breaks = y_axis_labels) +
  theme_classic() +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.ticks = element_line(color = NA),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=12))

ggsave(sprintf("%s%s_%s_%s_%sbp.png", image_dir, Sys.Date(), sample_id, ref, window),
       p, width = 18, height = 1.7, units = "in", device = png, dpi = 300, bg = "white")

# Save dataframes as excel----
outfiles <- list(plotting_data=genome_depth,
                 read_depth_summary=chr_median,
                 raw_genome_median=as.data.frame(raw_genome_median),
                 corrected_genome_median=as.data.frame(genome_median))

write_xlsx(outfiles, path = sprintf("%s%s_%s.xlsx", data_dir, Sys.Date(), sample_id))
