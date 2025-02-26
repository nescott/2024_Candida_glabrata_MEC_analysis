#' Draw Linear Synteny Plots
#'
#' This function draws linear synteny plots.
#'
#' It requires:
#'
#' 1. The desired output file name;
#'
#' 2. Tab separated file of all chromosome, scaffold, or contig lengths and the species identifier,
#' in order from first target species in the alignment files followed by the first reference species in the alignment files
#' -- top of file -- to the last target species and reference species in the alignment files -- end of file.
#' in this format:
#' chromosome ID, chromosome length, species identifier
#'
#' 3. files containing the syntenic blocks - one file per alignment, in order from first target/reference
#' (most recent species pairwise alignment in ancestral reconstruction data) alignment file
#' to last target/reference (ancestor pairwise alignment in ancestral reconstruction data) alignment file
#' following this format:
#' reference chromosome, reference start position, reference end position, target chromosome,
#' target start position, target end position, orient, reference species identifier, target species identifier
#'
#' Please separate files by tab and ensure any species identifiers used between length and alignment files are matching (same identifiers and caseing)
#'
#'
#' There are optional parameters for some customization of this function:
#'
#' 1. The format for saving the image i.e. png or pdf can be altered by inputting: `fileformat = "pdf"` (the default value is "png")
#'
#' 2. The colour of the synteny bands can be altered by inputting a concatenated string of chromosome IDs with assigned colour values which can be found with R colour Pallette
#' e.g. `colours = c("1" = "red", "2" = "blue", "3" = "green","4" = "orange", "5" = "purple","X" = "grey")` if no colours are assigned default values will be used but colours MUST be assigned to all chromosomes
#'
#' 3. The width of the image created can be changed by inputting: `w = 13` (default)
#'
#' 4. The height of the image created can be changed by inputting: `h = 5` (default)
#'
#' 5. The opacity of the ribbons can be changed by inputting: `opacity = .5` (default)
#'
#' 6. The directory where the image file should be saved, as default the image is saved to temporary directory, change by inputting: `directory = "path/to/directory"`
#'
#' The function works using the chromosome length file to order the Y axis and provide chromosome lengths to draw chromosome ideograms and the alignment files provides coordinates to draw the alignment bands between ideograms
#'
#' Example: `draw.linear("outputname", "example_lengths.txt", "example_alignment_1.txt", "example_alignment_2.txt", "example_alignment_3.txt", directory = "path/to/directory", fileformat = "pdf")`
#'
#' @title Linear synteny plot
#' @param output output file name
#' @param sizefile Chromosome Size file
#' @param ... synteny files (any number of alignment files can be entered)
#' @param directory string containing file path to chosen directory to save image file
#' @param fileformat output file format specified using the format `fileformat = "pdf"` (the default is "png")
#' @param colours concatenated string of chromosome IDs and assigned colours if desired using the format `colours = c("1" = "red", "2" = "blue", "3" = "green", "X" = "grey")` if the no colours are assigned default values will be used
#' @param w width of output image using the format `w = 13` (default)
#' @param h height of output image using the format `h = 5` (default)
#' @param opacity opacity of syntenic bands using the format `opacity = .5` (default)
#' @return An image file showing the linear comparison drawings
#' @examples
#'
#' # Create objects containing file paths to external dataset
#' # (see vignette to follow examples with personal data)
#'
#' length.file <- system.file("extdata", "example_lengths.txt", package = "syntenyPlotteR")
#' file1 <- system.file("extdata", "example_alignment_1.txt", package = "syntenyPlotteR")
#' file2 <- system.file("extdata", "example_alignment_2.txt", package = "syntenyPlotteR")
#' file3 <- system.file("extdata", "example_alignment_3.txt", package = "syntenyPlotteR")
#'
#' # -----------------------------------------------------------------------------------
#'
#' # Run draw.linear function
#' # To run example and save file to working directory
#' # add directory parameter and set working directory
#' # To run example with personal data see vignette
#'
#' draw.linear("outputName", length.file, file1, file2, file3, fileformat = "pdf")
#' @export
#'

library(tidyverse)
draw.linear <- function(output, sizefile, ..., directory = NULL, fileformat = "png", colours = colours.default, u = "in",w = 6, h = 5, opacity = .6) {
  
  if (is.null(directory)) {
    directory <- tempdir()
  }
  
  synteny.data.reframing <- function(data, tar.y, ref.y, compiled.size) {
    synteny <- data.frame()
    for (i in c(1:nrow(data))) {
      reference <- data[i, "ref.species"]
      target <- data[i, "tar.species"]
      tar_chr <- data[i, "tarchr"]
      ref_chr <- data[i, "refchr"]
      dir <- data[i, "dir"]
      tar_sizes <- compiled.size[compiled.size$species == target, ]
      names(tar_sizes) <- c("tarchr", "size", "species", "xstart", "xend")
      ref_sizes <- compiled.size[compiled.size$species == reference, ]
      names(ref_sizes) <- c("refchr", "size", "species", "xstart", "xend")
      tar_add <- tar_sizes[as.character(tar_sizes$tarchr) == as.character(tar_chr), ]$xstart
      ref_add <- ref_sizes[as.character(ref_sizes$refchr) == as.character(ref_chr), ]$xstart
      tar_y <- tar.y
      ref_y <- ref.y
      tar_xstart <- data[i, "tarstart"] + tar_add
      tar_xend <- data[i, "tarend"] + tar_add
      ref_xstart <- data[i, "refstart"] + ref_add
      ref_xend <- data[i, "refend"] + ref_add
      
      inverted <- grepl("-", dir, fixed = TRUE)
      if (inverted == TRUE) {
        df <- data.frame(
          x = c(tar_xstart, tar_xend, ref_xstart, ref_xend), y = c(tar_y, tar_y, ref_y, ref_y),
          fill = ref_chr, group = paste0("s", i), ref = reference, tar = target
        )
      } else {
        df <- data.frame(
          x = c(tar_xstart, ref_xstart, ref_xend, tar_xend), y = c(tar_y, ref_y, ref_y, tar_y),
          fill = ref_chr, group = paste0("s", i), ref = reference, tar = target
        )
      }
      synteny <- rbind(synteny, df)
    }
    return(synteny)
  }
  
  colours.default <- c(
    "1" = "#BFD73B", "2" = "#39ACE2", "3" = "#F16E8A",
    "4" = "#2DB995", "5" = "#855823", "6" = "#A085BD",
    "7" = "#2EB560", "8" = "#D79128", "9" = "#FDBB63",
    "10" = "#AFDFE5", "11" = "#BF1E2D", "12" = "purple4",
    "13" = "#B59F31", "14" = "#F68B1F", "15" = "#EF374B",
    "16" = "#D376FF", "17" = "#009445", "18" = "#CE4699",
    "19" = "#7C9ACD", "20" = "#84C441", "21" = "#404F23",
    "22" = "#607F4B", "23" = "#EBB4A9", "24" = "#F6EB83",
    "25" = "#915F6D", "26" = "#602F92", "27" = "#81CEC6",
    "28" = "#F8DA04", "29" = "peachpuff2", "30" = "gray85", "33" = "peachpuff3",
    "W" = "#9590FF", "Z" = "#666666", "Y" = "#9590FF", "X" = "#666666",
    "LGE22" = "grey", "LGE64" = "gray64",
    "1A" = "pink", "1B" = "dark blue", "4A" = "light green",
    "Gap" = "white", "LG2" = "black", "LG5" = "#CC99CC"
  )
  
  xstart <- xend <- refchr <- tarchr <- x <- y <- group <- fill <- chromosome <- species <- NULL
  sizes <- utils::read.delim(sizefile, header = FALSE) # to be consistent with naming in EH
  names(sizes) <- c("chromosome", "size", "species")
  sizes$size <- as.numeric(gsub(",", "", sizes$size))
  
  count <- 0
  compiled.size <- data.frame()
  for (i in unique(sizes$species)) {
    size.intermediate <- sizes[sizes$species == i, ]
    for (x in c(1:nrow(size.intermediate))) {
      if (x == 1) {
        total_start <- 1
        total_end <- size.intermediate[x, "size"]
      } else {
        total_start <- total_end + 20000
        total_end <- total_start + size.intermediate[x, "size"]
      }
      size.intermediate[x, "xstart"] <- total_start
      size.intermediate[x, "xend"] <- total_end
    }
    compiled.size <- rbind(compiled.size, size.intermediate)
  }
  
  for (z in unique(compiled.size$species)) {
    compiled.size$y[compiled.size$species == z] <- count
    count <- count + 2
  }
  
  list.of.files <- list()
  for (i in list(...)) {
    list.of.files[[i]] <- i
  }
  
  listsynt <- list()
  for (i in 1:length(list.of.files)) {
    num <- i
    file <- list.of.files[[num]]
    dataTMP <- utils::read.delim(file, header = FALSE)
    data2 <- dataTMP[, c(4, 5, 6, 1, 2, 3, 7, 8, 9)]
    colnames(data2) <- c("tarchr", "tarstart", "tarend", "refchr", "refstart", "refend", "dir", "ref.species", "tar.species")
    data2$tarstart <- as.numeric(gsub(",", "", data2$tarstart))
    data2$tarend <- as.numeric(gsub(",", "", data2$tarend))
    data2$refstart <- as.numeric(gsub(",", "", data2$refstart))
    data2$refend <- as.numeric(gsub(",", "", data2$refend))
    reference <- data2[1, "ref.species"]
    target <- data2[1, "tar.species"]
    ref_y <- compiled.size[compiled.size$species == reference, "y"]
    tar_y <- compiled.size[compiled.size$species == target, "y"]
    if (tar_y[1] > ref_y[1]){
      ref_y <- ref_y[1] + 0.1
      tar_y <- tar_y[1]
    } else{
      ref_y <- ref_y[1]
      tar_y <- tar_y[1] + 0.1
    }
    x <- synteny.data.reframing(data2, tar_y, ref_y, compiled.size)
    x$fill <- as.factor(x$fill)
    listsynt[[i]] <- x
  }
  
  compiled.size$chromosome <- as.factor(compiled.size$chromosome)
  
  p <- ggplot2::ggplot()
  
  for (i in 1:length(listsynt)) {
    data <- listsynt[[i]]
    reference <- data[1, "ref"]
    target <- data[1, "tar"]
    ref_sizes <- compiled.size[compiled.size$species == reference, ]
    tar_sizes <- compiled.size[compiled.size$species == target, ]
    p <- p + ggplot2::geom_rect(
      data = ref_sizes, mapping = ggplot2::aes(xmin = xstart, xmax = xend, ymin = y, ymax = y + 0.10, fill = chromosome),
      color = "black", alpha = 0.85, linewidth = 0.2
    ) +
      #ggplot2::geom_text(data = ref_sizes, ggplot2::aes(x = (xstart + xend) / 2, y = y + 0.2, label = chromosome), size = 2, angle = 45) +
      #ggplot2::geom_text(data = ref_sizes, mapping = ggplot2::aes(x = 2, y = y, label = species), size = 3, hjust = 1) +
      ggplot2::geom_rect(
        data = tar_sizes, mapping = ggplot2::aes(xmin = xstart, xmax = xend, ymin = y, ymax = y + 0.10), fill = "grey85",
        color = "black", alpha = 0.85, linewidth = 0.2
      ) +
      #ggplot2::geom_text(data = tar_sizes, ggplot2::aes(x = (xstart + xend) / 2, y = y + 0.2, label = chromosome), size = 2, angle = 45) +
      #ggplot2::geom_text(data = tar_sizes, mapping = ggplot2::aes(x = 2, y = y+1, label = species), size = 2, hjust = 1.5) +
      ggplot2::geom_polygon(data = data, alpha = opacity, ggplot2::aes(x = x, y = y, group = group, fill = fill))
  }
  
  p <- p + ggplot2::scale_fill_manual(values = colours) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "none"
    )
  
  spp <- compiled.size %>% group_by(species) %>% summarise(y=unique(y)) %>% 
    mutate(id = c("CBS138", "1-1", "74", "30-1","88", "107"))
  
  chrs <- compiled.size %>% filter(y==0) %>% 
    mutate(short = LETTERS[1:13])
  
  p <- p +
    geom_text(data = spp, aes(x = 2, y = y + 0.07, label = id), hjust = 1.1, color = "black", size = 3) +
    geom_text(data = chrs, aes(x = (xstart + xend)/2, y = y-0.5, label = short), color = "black", size = 3) +
    expand_limits(x = -400000)
    
  
  save(listsynt, file = "listsynt.rda")
  save(sizes, file= "sizes.rda")
  save(compiled.size, file="compiledsize.rda")
  save(p, file="p.rda")
  message(paste0("Saving linear image to ", directory))
  print(p)
  ggplot2::ggsave(paste0(directory,"/",output, ".", fileformat), p, device = fileformat, width = w, height = h, units = u)
}
