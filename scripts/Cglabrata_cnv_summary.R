dels <- read_excel("ch4/Cglabrata_MEC_deletion_cov_reviewed.xlsx")
dups <- read_excel("ch4/Cglabrata_MEC_duplication_cov_reviewed.xlsx")
samples <- read_xlsx("~/Google Drive/My Drive/dissertation/04 Genomic plasticity of Candida glabrata bloodstream isolates/Ch4_supp_data.xlsx")

dels <- dels %>% 
  left_join(samples %>% select(`Study code`,`Sequence type`, Isolate), by=c("sample"="Study code")) %>% 
  mutate(mean_rel_depth = round(mean_rel_depth, digits = 2))
del_sts <- dels %>% group_by(ID, `Sequence type`) %>% count()

dups <- dups %>% 
  left_join(samples %>% select(`Study code`,`Sequence type`, Isolate), by=c("sample"="Study code"))%>% 
  mutate(mean_rel_depth = round(mean_rel_depth, digits = 2))
dup_sts <- dups %>% group_by(ID, `Sequence type`) %>% count()


del_summary <- dels %>% group_by(ID) %>% count()
del_summary <- del_summary %>% 
  inner_join(dels %>% group_by(ID) %>% summarize(gene_size = unique(width)))
del_summary <- del_summary %>% 
  left_join(dels %>% group_by(ID) %>% summarize(min_cov = min(mean_rel_depth), 
                                                max_cov = max(mean_rel_depth), 
                                                avg_cov = mean(mean_rel_depth))) %>% 
  left_join(dels %>% select(ID,Gene)) %>% 
  distinct()

dup_summary <- dups %>% group_by(ID) %>% count()
dup_summary <- dup_summary %>% 
  inner_join(dups %>% group_by(ID) %>% summarize(gene_size = unique(width)))
dup_summary <- dup_summary %>% 
  left_join(dups %>% group_by(ID) %>% summarize(min_cov = min(mean_rel_depth), 
                                                max_cov = max(mean_rel_depth), 
                                                avg_cov = mean(mean_rel_depth))) %>% 
  left_join(dups %>% select(ID, Gene)) %>% 
  distinct()

write_xlsx(dels, "Cglabrata_MEC_dels_for_supp_table.xlsx")
write_xlsx(del_summary, "Cglabrata_MEC_del_summary_writeup.xlsx")
write_xlsx(dups, "Cglabrata_MEC_dups_for_supp_table.xlsx")
write_xlsx(dup_summary, "Cglabrata_MEC_dup_summary_writeup.xlsx")
