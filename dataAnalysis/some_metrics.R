## specific to GTDB and not present with Silva

all_found_with_silva <- 
  species_cleaned %>% 
  filter(Taxonomic_level == "Species",
         Abundance > 0) %>%
  mutate(Species_name = ifelse(str_detect(Taxonomy,"sp"),
                               "GTDB placeholder",
                               Taxonomy)) %>%
  mutate(Sequencing_technology = 
           ifelse(Sequencing_technology == "Illumina",
                  "V4-V5 16S rRNA gene",
                  "full-length 16S rRNA gene")) %>%  
  group_by(Sample, Database, Sequencing_technology) %>%
  mutate(RelativeAbundance = Abundance*100/sum(Abundance)) %>% 
  ungroup() %>% 
  filter(RelativeAbundance >= 0.1, !is.na(Taxonomy)) %>% 
  filter(Database == "Silva") %>% 
  pull(Species_name)



species_cleaned %>% 
  filter(Taxonomic_level == "Species",
         Abundance > 0) %>%
  mutate(Species_name = ifelse(str_detect(Taxonomy,"sp"),
                               "GTDB placeholder",
                               Taxonomy)) %>%
  mutate(Sequencing_technology = 
           ifelse(Sequencing_technology == "Illumina",
                  "V4-V5 16S rRNA gene",
                  "full-length 16S rRNA gene")) %>%  
  group_by(Sample, Database, Sequencing_technology) %>%
  mutate(RelativeAbundance = Abundance*100/sum(Abundance)) %>% 
  ungroup() %>% 
  filter(RelativeAbundance >= 0.1, !is.na(Taxonomy)) %>% 
  filter(!Species_name %in% all_found_with_silva, Database == "GTDB", Sequencing_technology == "full-length 16S rRNA gene") %>% 
  pull(Species_name) %>% 
  unique() %>% length()


###
species_cleaned %>% 
  filter(Taxonomic_level == "Species",
         Abundance > 0) %>%
  mutate(Species_name = ifelse(str_detect(Taxonomy,"sp"),
                               "GTDB placeholder",
                               Taxonomy)) %>%
  mutate(Sequencing_technology = 
           ifelse(Sequencing_technology == "Illumina",
                  "V4-V5 16S rRNA gene",
                  "full-length 16S rRNA gene")) %>%  
  group_by(Sample, Database, Sequencing_technology) %>%
  mutate(RelativeAbundance = Abundance*100/sum(Abundance)) %>% 
  ungroup() %>% 
  filter(RelativeAbundance >= 0.1, !is.na(Taxonomy)) %>% 
  pull(Sequence) %>% unique() %>% length()





