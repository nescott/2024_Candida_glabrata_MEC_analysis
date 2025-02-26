chrs <- c("ChrA" = "#BFD73B", "ChrB" = "#39ACE2", "ChrC" = "#F16E8A",
          "ChrD" = "#2DB995", "ChrE" = "#855823", "ChrF" = "#A085BD",
          "ChrG" = "#2EB560", "ChrH" = "#D79128", "ChrI" = "#FDBB63",
          "ChrJ" = "#AFDFE5", "ChrK" = "#BF1E2D", "ChrL" = "purple4",
          "ChrM" = "#B59F31", "mito" = "grey")


draw.linear(paste0(Sys.Date(),"_Cglabrata_MEC_synteny"), 
            "~/umn/data/assemblies/sorted_assembly_lengths.txt", 
            "~/umn/data/assemblies/filtered_align1.txt",
            "~/umn/data/assemblies/filtered_align2.txt", 
            "~/umn/data/assemblies/filtered_align3.txt", 
            "~/umn/data/assemblies/filtered_align4.txt", 
            "~/umn/data/assemblies/filtered_align5.txt", 
            directory = "~/umn/thesis",
            fileformat = "pdf", colours = chrs, w=9, h=6)


chr_vals <- c("ChrA" = "#1F77B4FF", "ChrB" = "#FF7F0EFF", "ChrC" = "#2CA02CFF",
              "ChrD" = "#D62728FF", "ChrE" = "#9467BDFF", "ChrF" = "#8C564BFF",
              "ChrG" = "#E377C2FF", "ChrH" = "#7F7F7FFF", "ChrI" = "#BCBD22FF",
              "ChrJ" = "#17BECFFF", "ChrK" = "#AEC7E8FF", "ChrL" = "#FFBB78FF",
              "ChrM" = "#98DF8AFF", "mito" = "lightgrey")
