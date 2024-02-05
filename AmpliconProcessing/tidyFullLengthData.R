#
library(tidyr)

# Organize ASV abundance table
ASV_abundance_full_mosj <- seq_tab_nochim_full_mosj %>% t() %>% as.data.frame()
ASV_abundance_full_mosj$Sequence <- rownames(ASV_abundance_full_mosj)
rownames(ASV_abundance_full_mosj) <- NULL

# Silva taxonomy (full)
ASV_taxonomy_silva_full_mosj_df <- data.frame(ASV_taxonomy_silva_full_mosj)
ASV_taxonomy_silva_full_mosj_df$Sequence <- rownames(ASV_taxonomy_silva_full_mosj_df)
rownames(ASV_taxonomy_silva_full_mosj_df) <- NULL

ASV_abundance_taxonomy_full_mosj_df_silva <- ASV_abundance_full_mosj %>% 
  left_join(ASV_taxonomy_silva_full_mosj_df, by = "Sequence") %>% 
  mutate(Sequencing_technology = "PacBio",
         Database = "Silva") %>% 
  mutate(Species = ifelse(is.na(Species), NA, paste(Genus, Species)))

# GTDB taxonomy (full)
ASV_taxonomy_gtdb_full_mosj_df <- data.frame(ASV_taxonomy_gtdb_full_mosj)
ASV_taxonomy_gtdb_full_mosj_df$Sequence <- rownames(ASV_taxonomy_gtdb_full_mosj_df)
rownames(ASV_taxonomy_gtdb_full_mosj_df) <- NULL

ASV_abundance_taxonomy_full_mosj_df_gtdb <- ASV_abundance_full_mosj %>% 
  left_join(ASV_taxonomy_gtdb_full_mosj_df, by = "Sequence") %>% 
  mutate(Sequencing_technology = "PacBio",
         Database = "GTDB")

#
ASV_full_length_mosj <- ASV_abundance_taxonomy_full_mosj_df_silva %>% 
  rbind(ASV_abundance_taxonomy_full_mosj_df_gtdb)

#
mosj_2019_pacbio_tidy <- ASV_full_length_mosj %>% 
  #mutate(Species = ifelse(is.na(Species), NA, paste(Genus, Species))) %>% 
  pivot_longer(cols = contains("M19"), names_to = "Sample", values_to = "Abundance") %>% 
  pivot_longer(cols = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), 
               names_to = "Taxonomic_level", values_to = "Taxonomy")

write.csv(mosj_2019_pacbio_tidy, file = "./data/mosj_2019_pacbio_tidy")
