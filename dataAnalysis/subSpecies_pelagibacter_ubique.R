## pelagibacter ubique
pelagibacter_ubique_asvs <- 
  ASVs_merged_good_samples %>% 
  filter(Taxonomic_level == "Species", 
         str_detect(Taxonomy, "Pelagibacter ubique"), 
         Abundance > 0) %>% 
  select(Sequence, Taxonomy) %>% 
  distinct() %>% 
  mutate(ID = paste("ASV", rownames(.), sep = "_")) %>% 
  dplyr::rename(seq = Sequence,
         name = ID) %>% 
  as.data.frame()



#library(ShortRead)

write.fasta(as.list(pelagibacter_ubique_asvs$seq),
            pelagibacter_ubique_asvs$name,
            file = "./pelagibacter_ubique_asvs.fasta")
