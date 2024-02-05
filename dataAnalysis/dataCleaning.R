#
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
library(vegan)
# Change based on your file location
#mosj_2019_illumina_tidy <- read.csv("./data/mosj_2019_illumina_tidy")
#mosj_2019_pacbio_tidy <- read.csv("./data/mosj_2019_pacbio_tidy")

## Merge PacBio and Illumina data in a single, tidy data.frame
ASVs_merged <-
  mosj_2019_illumina_tidy %>%
  bind_rows(mosj_2019_pacbio_tidy)

# Clean Sample ids in ASVs merged to filter for the selected samples
ASVs_merged <- ASVs_merged %>%
  mutate(Sample = str_replace(Sample, fixed("-"),"_"),
         Sample = str_extract(Sample, "M19_\\d+"))
#
# Load metadata info and clean
metadata <- read.table("./data/metadata", sep = "\t", header = TRUE)

# remove lines with only NA
metadata <- metadata %>% filter(!is.na(NGS_code))

# change according to your file location
#selected_samples_mosj_2019 <- read.table("./data/selected_samples", header = TRUE)

# select relevant samples and variables
mosj_env <- metadata %>%
  filter(NGS_code %in% selected_samples_mosj_2019$Sample) %>%
  select(Station, Depth_m, NGS_code)

# Add environmental data to ASVs_merged
ASVs_merged <-
  ASVs_merged %>%
  left_join(mosj_env, by = c("Sample" = "NGS_code"))

## Taxonomy cleaning
# remove chloroplasts, mithocondria and eukarya, if any (there are none)
ASVs_merged <- ASVs_merged %>%
  filter(!Taxonomy %in% c("Chloroplast", "Mithocondria", "Eukaryota"))

# Make sure the factors are correctly leveled
ASVs_merged <- ASVs_merged %>%
  mutate(Taxonomic_level = factor(Taxonomic_level,
                                  levels = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")))

# Samples to remove for rarefaction (<10000 reads)
samples_to_remove <-
  ASVs_merged %>%
  filter(Taxonomic_level == "Kingdom") %>%
  group_by(Sample, Sequencing_technology, Database) %>%
  summarise(total_reads = sum(Abundance)) %>%
  filter(total_reads < 10000) %>%
  pull(Sample)

## Rarefied to 10 000 reads
ASVs_with_rarefaction <-
  ASVs_merged %>%
  filter(!Sample %in% samples_to_remove) %>%
  group_by(Sample, Sequencing_technology, Database, Taxonomic_level) %>%
  nest() %>%
  mutate(Rarefied_reads = map(.x = data, ~as.data.frame(t(rrarefy(.x$Abundance, sample = 10000))))) %>%
  unnest(c(data, Rarefied_reads))%>%
  dplyr::rename(Rarefied_abundance = "V1")

ASVs_merged_good_samples <- ASVs_merged %>% 
  filter(!Sample %in% samples_to_remove)

