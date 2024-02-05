## sub-species diversity correction ##

## note: this is only for PacBio + GTDB
## corrected, not rarefied

corrected_subspecies <- 
  ASVs_merged_good_samples %>% 
    select(-X) %>% 
    filter(!is.na(Taxonomy)) %>% 
  filter(Abundance > 0) %>% 
  filter(Taxonomic_level == "Species", 
         Database == "GTDB", 
         Sequencing_technology == "PacBio") %>% 
  group_by(Sample, Taxonomy, Station, Depth_m) %>% 
  summarise(correctedAbundance = sum(Abundance)) %>% 
  ungroup() %>% 
  group_by(Sample, Station, Depth_m) %>% 
  summarise(correctedRichness = length(unique(Taxonomy))) %>% 
  mutate(correctedSubSpecies = "Yes")
  

uncorrected_subspecies <- ASVs_merged_good_samples %>% 
  select(-X) %>% 
  filter(!is.na(Taxonomy)) %>% 
  filter(Abundance > 0) %>% 
  filter(Taxonomic_level == "Species", 
         Database == "GTDB", 
         Sequencing_technology == "PacBio") %>% 
  group_by(Sample, Station, Depth_m) %>% 
  summarise(correctedRichness = length(unique(Taxonomy))) %>% 
  mutate(correctedSubSpecies = "No")


correct_and_uncorrect_subspecies <-
  corrected_subspecies %>% 
  rbind(uncorrected_subspecies)


##################

sub_species_diversity_summary <- 
  ASVs_merged_good_samples %>%
  select(-X) %>% 
  filter(Abundance > 0) %>% 
  filter(Taxonomic_level == "Species", 
         Database == "GTDB", 
         Sequencing_technology == "PacBio") %>%
  mutate(Species_name = ifelse(str_detect(Taxonomy,"sp"),
                               "GTDB placeholder",
                               Taxonomy)) %>% 
  select(Sequence, Species_name) %>% 
  filter(!is.na(Species_name), 
         Species_name != "GTDB placeholder") %>% 
  distinct() %>% 
  pull(Species_name) %>% 
  table() %>% 
  as.data.frame() %>% 
  arrange(desc(Freq))


sub_species_diversity_summary %>% filter(Freq > 2) %>% 
  ggplot(aes(x = reorder(., Freq), 
             y = Freq)) + 
  geom_col() + 
  coord_flip() + 
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(face = "italic"),
        legend.position = "top",
        text = element_text(size = 15),
        strip.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.text.y.left = element_text(size = 8)) + 
  labs(x = "Species with binomial name",
       y = "Number of ASVs",
       title = "Species with binomial name and more than 1 ASV")
  



sub_species_diversity_summary %>% 
  filter(Freq > 1) %>% 
  write.csv("./data/subSpeciesOcurrences.csv")




