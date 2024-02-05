## Main analysis
library(RColorBrewer)
library(vegan)
library(tidyr)
library(stringr)
library(dplyr)
library(ggplot2)
#
qualitative_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

## some pre-processing for figures
percent_classification <- 
  ASVs_merged_good_samples %>%
  filter(Abundance > 0) %>%
  group_by(Sample, Taxonomic_level, Sequencing_technology, Database) %>%
  summarise(Percent_classified = mean(!is.na(Taxonomy))*100)
#
percent_classification_rarefied <- 
  ASVs_with_rarefaction %>%
  mutate(Abundance = Rarefied_abundance) %>% 
  filter(Abundance > 0) %>%
  group_by(Sample, Taxonomic_level, Sequencing_technology, Database) %>%
  summarise(Percent_classified = mean(!is.na(Taxonomy))*100)

#
diversity_count <- 
  ASVs_merged_good_samples %>% 
  filter(Abundance > 0) %>% 
  group_by(Sample, Taxonomic_level,Sequencing_technology,Database) %>%
  summarise(n_classified_ASVs = sum(!is.na(Taxonomy)))

#
diversity_count_rarefied <- 
  ASVs_with_rarefaction %>% 
  mutate(Abundance = Rarefied_abundance) %>%  
  filter(Abundance > 0) %>% 
  group_by(Sample, Taxonomic_level, Sequencing_technology, Database) %>%
  summarise(n_classified_ASVs = sum(!is.na(Taxonomy)))

#
species_cleaned <- 
  ASVs_merged_good_samples %>% 
  mutate(Taxonomy = str_replace(Taxonomy, "\\s*\\([^\\)]+\\)", ""),
         Taxonomy = str_remove_all(Taxonomy, "_.")) %>%
  # Additional manual curation
  mutate(Taxonomy = ifelse(Taxonomy == "Pseudomonas stutzeriB", "Pseudomonas stutzeri", Taxonomy))

#
species_cleaned_rarefied <- 
  ASVs_with_rarefaction %>% 
  filter(Rarefied_abundance > 0) %>% 
  mutate(Taxonomy = str_replace(Taxonomy, "\\s*\\([^\\)]+\\)", ""),
         Taxonomy = str_remove_all(Taxonomy, "_.")) %>%
  # Additional manual curation
  mutate(Taxonomy = ifelse(Taxonomy == "Pseudomonas stutzeriB", "Pseudomonas stutzeri", Taxonomy)) 


## Figures ##
## Figure 2
percent_classification %>% 
  mutate(Sequencing_technology = 
           ifelse(Sequencing_technology == "Illumina",
                  "V4-V5 16S rRNA gene",
                  "full-length 16S rRNA gene")) %>%  
  ggplot(aes(x = Taxonomic_level,
             y = Percent_classified,
             fill = Sequencing_technology))+
  facet_wrap(~Database) + 
  geom_boxplot(width = 0.5, outlier.shape = "cross")+
  theme_bw()+
  theme(axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        strip.background = element_blank(),
        text = element_text(size = 15),
        legend.text = element_text(size = 15),
        panel.grid.major.y = element_line(linetype = "dashed",colour = "grey"),
        #panel.grid.major.x = element_line(colour = "grey")
        )+
  labs(x = "Taxonomic level",
       y = "Percentage of ASVs classified (%)",
       fill = "Sequencing of ",
       subtitle = paste("n =", length(unique(pull(percent_classification,Sample))))) +
  scale_fill_manual(values = qualitative_colors[c(1:2)])


## Supplementary Figure S2
percent_classification_rarefied %>% 
  mutate(Sequencing_technology = 
           ifelse(Sequencing_technology == "Illumina",
                  "V4-V5 16S rRNA gene",
                  "full-length 16S rRNA gene")) %>%  
  ggplot(aes(x=Taxonomic_level,y = Percent_classified,
             fill=Sequencing_technology))+
  geom_boxplot(outlier.shape = "cross", width = 0.5)+
  facet_wrap(~Database)+
  theme_bw()+
  theme(axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        panel.grid.major.y = 
          element_line(linetype = "dashed", colour = "grey"),
        strip.text = element_text(size = 15),
        text = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.background = element_blank())+
  labs(x = "Taxonomic level",
       y = "Percentage of ASVs classified (%)",
       fill = "Sequencing of ",
       subtitle = "Samples rarefied to 10 000 reads \nn = 13") + 
  scale_fill_manual(values = qualitative_colors[c(1,2)])
#

## Figure 3
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
  ggplot(aes(reorder(Species_name, RelativeAbundance), 
             RelativeAbundance, 
             col = Sequencing_technology))+
  geom_jitter(width = 0.1, size = 1)+
  facet_wrap(~Database)+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(face = "italic"),
        legend.position = "top",
        text = element_text(size = 15),
        strip.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.text.y.left = element_text(size = 8))+
  labs(col = "Sequencing of ",
       y = "Relative abundance (%) in Log10 scale",
       x = "Species",
       subtitle = "Abundant ASVs (>0.1% relative abundance)"
       )+
  scale_y_log10()+
  scale_color_manual(values = qualitative_colors[c(1,2)])


## alpha diversity
gridExtra::grid.arrange(
  diversity_count %>% 
    mutate(Sequencing_technology = 
             ifelse(Sequencing_technology == "Illumina",
                    "V4-V5 16S rRNA gene",
                    "full-length 16S rRNA gene")) %>%  
    filter(Taxonomic_level == "Species", n_classified_ASVs > 0) %>% 
    left_join(mosj_env, by = c("Sample" = "NGS_code")) %>% 
    mutate(Station = factor(Station, levels = c("KB3", "KB6", "KB0", "V12", "V6", "HG-IV"))) %>% 
    mutate(Depth_factor = case_when(Depth_m < 10 ~ "<10",
                                    Depth_m < 100 ~ "10-100",
                                    TRUE ~ ">100")) %>% 
    mutate(Depth_factor = factor(Depth_factor, 
                                 levels = c("<10", "10-100", ">100"))) %>% 
    ggplot(aes(Station, n_classified_ASVs))+
    geom_point(aes(col = Depth_factor)) +
    theme_bw()+ 
    theme(panel.grid = element_blank(),
          legend.position = "top",
          strip.background = element_blank(),
          text = element_text(size = 15),
          legend.text = element_text(size = 15),
          panel.grid.major.y = element_line(linetype = "dashed",colour = "grey"),
          axis.text.x = element_text(size = 10)) + 
    facet_grid(facets = c("Database", "Sequencing_technology")) +
    scale_y_log10() + 
    labs(y = "Number of ASVs (Log10)",
         col = "Depth (m):", 
         title = "a",
         subtitle = "ASVs classified at species level") + 
    scale_color_manual(values = brewer.pal(n = 9, "Blues")[c(3,6,9)]),
  
  ## Rare biosphere section
  species_cleaned %>% 
    filter(!is.na(Taxonomy),
           Abundance > 1) %>% 
    group_by(Sample, Database, Sequencing_technology, Taxonomic_level, Depth_m) %>% 
    mutate(RelativeAbundance = Abundance * 100/ sum(Abundance)) %>% 
    mutate(rarity = ifelse(RelativeAbundance >= 0.9, "Abundant", "Rare")) %>% 
    ungroup() %>% 
    select(Sequencing_technology, Database, Sample, Taxonomic_level, rarity, Depth_m) %>% 
    group_by(Sequencing_technology, Database, rarity, Sample, Taxonomic_level, Depth_m) %>% 
    summarise(counts = n()) %>% 
    group_by(Sequencing_technology, Database, rarity, Taxonomic_level, Depth_m) %>% 
    filter(Taxonomic_level %in% c("Species")) %>% 
    mutate(Sequencing_technology = 
             ifelse(Sequencing_technology == "Illumina",
                    "V4-V5 16S rRNA gene",
                    "full-length 16S rRNA gene"))  %>% 
    mutate(Depth_factor = case_when(Depth_m < 10 ~ "<10",
                                    Depth_m < 100 ~ "10-100",
                                    TRUE ~ ">100")) %>% 
    mutate(Depth_factor = factor(Depth_factor, 
                                 levels = c("<10", "10-100", ">100"))) %>% 
    ggplot(aes(rarity, counts))+
    geom_boxplot(position = "dodge", outlier.shape = "cross") + 
    geom_jitter(aes(col = Depth_factor)) + 
    facet_grid(facets = c("Database", "Sequencing_technology")) + 
    theme_bw() + 
    labs(title = "b",
         x = "Abundance classification",
         y = "",
         col = "Depth (m):",
         subtitle = "ASVs classified at species level") + 
    theme(panel.grid = element_blank(),
          legend.position = "top",
          strip.background = element_blank(),
          text = element_text(size = 15),
          legend.text = element_text(size = 15),
          panel.grid.major.y = element_line(linetype = "dashed",colour = "grey"))+
    scale_y_log10() + 
    scale_color_manual(values = brewer.pal(n = 9, "Blues")[c(3,6,9)])
  ,
  layout_matrix = rbind(c(1,2),c(1,2))
)

# Sup Figure S4
## alpha diversity, stations, rarefied, species level
gridExtra::grid.arrange(
  diversity_count_rarefied %>% 
    mutate(Sequencing_technology = 
             ifelse(Sequencing_technology == "Illumina",
                    "V4-V5 16S rRNA gene",
                    "full-length 16S rRNA gene")) %>%  
    filter(Taxonomic_level == "Species", n_classified_ASVs > 0) %>% 
    left_join(mosj_env, by = c("Sample" = "NGS_code")) %>% 
    mutate(Station = factor(Station, levels = c("KB3", "KB6", "KB0", "V12", "V6", "HG-IV"))) %>% 
    mutate(Depth_factor = case_when(Depth_m < 10 ~ "<10",
                                    Depth_m < 100 ~ "10-100",
                                    TRUE ~ ">100")) %>% 
    mutate(Depth_factor = factor(Depth_factor, 
                                 levels = c("<10", "10-100", ">100"))) %>% 
    ggplot(aes(Station, n_classified_ASVs))+
    geom_point(aes(col = Depth_factor)) +
    theme_bw()+ 
    theme(panel.grid = element_blank(),
          legend.position = "top",
          strip.background = element_blank(),
          text = element_text(size = 15),
          legend.text = element_text(size = 15),
          panel.grid.major.y = element_line(linetype = "dashed",colour = "grey"),
          axis.text.x = element_text(size = 10)) + 
    facet_grid(facets = c("Database", "Sequencing_technology")) +
    scale_y_log10() + 
    labs(y = "Number of ASVs (Log10)",
         col = "Depth (m):", 
         title = "A",
         subtitle = "ASVs classified at species level\nrarefied") + 
    scale_color_manual(values = brewer.pal(n = 9, "Blues")[c(3,6,9)]),
  
  ## Rare biosphere section
  species_cleaned_rarefied %>% 
    filter(!is.na(Taxonomy),
           Abundance > 1) %>% 
    group_by(Sample, Database, Sequencing_technology, Taxonomic_level, Depth_m) %>% 
    mutate(RelativeAbundance = Abundance * 100/ sum(Abundance)) %>% 
    mutate(rarity = ifelse(RelativeAbundance >= 0.9, "Abundant", "Rare")) %>% 
    ungroup() %>% 
    select(Sequencing_technology, Database, Sample, Taxonomic_level, rarity, Depth_m) %>% 
    group_by(Sequencing_technology, Database, rarity, Sample, Taxonomic_level, Depth_m) %>% 
    summarise(counts = n()) %>% 
    group_by(Sequencing_technology, Database, rarity, Taxonomic_level, Depth_m) %>% 
    filter(Taxonomic_level %in% c("Species")) %>% 
    mutate(Sequencing_technology = 
             ifelse(Sequencing_technology == "Illumina",
                    "V4-V5 16S rRNA gene",
                    "full-length 16S rRNA gene"))  %>% 
    mutate(Depth_factor = case_when(Depth_m < 10 ~ "<10",
                                    Depth_m < 100 ~ "10-100",
                                    TRUE ~ ">100")) %>% 
    mutate(Depth_factor = factor(Depth_factor, 
                                 levels = c("<10", "10-100", ">100"))) %>% 
    ggplot(aes(rarity, counts))+
    geom_boxplot(position = "dodge", outlier.shape = "cross") + 
    geom_jitter(aes(col = Depth_factor)) + 
    facet_grid(facets = c("Database", "Sequencing_technology")) + 
    theme_bw() + 
    labs(title = "B",
         x = "Abundance classification",
         y = "",
         col = "Depth (m):",
         subtitle = "ASVs classified at species level\nrarefied") + 
    theme(panel.grid = element_blank(),
          legend.position = "top",
          strip.background = element_blank(),
          text = element_text(size = 15),
          legend.text = element_text(size = 15),
          panel.grid.major.y = element_line(linetype = "dashed",colour = "grey"))+
    scale_y_log10() + 
    scale_color_manual(values = brewer.pal(n = 9, "Blues")[c(3,6,9)])
  ,
  layout_matrix = rbind(c(1,2),c(1,2))
)

## Phylum level alpha diversity (for supplementary)
# Supplementary Figure S5
gridExtra::grid.arrange(
  diversity_count %>% 
    mutate(Sequencing_technology = 
             ifelse(Sequencing_technology == "Illumina",
                    "V4-V5 16S rRNA gene",
                    "full-length 16S rRNA gene")) %>%  
    filter(Taxonomic_level == "Phylum", n_classified_ASVs > 0) %>% 
    left_join(mosj_env, by = c("Sample" = "NGS_code")) %>% 
    mutate(Station = factor(Station, levels = c("KB3", "KB6", "KB0", "V12", "V6", "HG-IV"))) %>% 
    mutate(Depth_factor = case_when(Depth_m < 10 ~ "<10",
                                    Depth_m < 100 ~ "10-100",
                                    TRUE ~ ">100")) %>% 
    mutate(Depth_factor = factor(Depth_factor, 
                                 levels = c("<10", "10-100", ">100"))) %>% 
    ggplot(aes(Station, n_classified_ASVs))+
    geom_point(aes(col = Depth_factor)) +
    theme_bw()+ 
    theme(panel.grid = element_blank(),
          legend.position = "top",
          strip.background = element_blank(),
          text = element_text(size = 15),
          legend.text = element_text(size = 15),
          panel.grid.major.y = element_line(linetype = "dashed",colour = "grey"),
          axis.text.x = element_text(size = 10)) + 
    facet_grid(facets = c("Database", "Sequencing_technology")) +
    scale_y_log10() + 
    labs(y = "Number of ASVs (Log10)",
         col = "Depth (m):", 
         title = "a",
         subtitle = "ASVs classified at phylum level") + 
    scale_color_manual(values = brewer.pal(n = 9, "Blues")[c(3,6,9)]),
  
  ## Rare biosphere section
  species_cleaned %>% 
    filter(!is.na(Taxonomy),
           Abundance > 1) %>% 
    group_by(Sample, Database, Sequencing_technology, Taxonomic_level, Depth_m) %>% 
    mutate(RelativeAbundance = Abundance * 100/ sum(Abundance)) %>% 
    mutate(rarity = ifelse(RelativeAbundance >= 0.9, "Abundant", "Rare")) %>% 
    ungroup() %>% 
    select(Sequencing_technology, Database, Sample, Taxonomic_level, rarity, Depth_m) %>% 
    group_by(Sequencing_technology, Database, rarity, Sample, Taxonomic_level, Depth_m) %>% 
    summarise(counts = n()) %>% 
    group_by(Sequencing_technology, Database, rarity, Taxonomic_level, Depth_m) %>% 
    filter(Taxonomic_level %in% c("Phylum")) %>% 
    mutate(Sequencing_technology = 
             ifelse(Sequencing_technology == "Illumina",
                    "V4-V5 16S rRNA gene",
                    "full-length 16S rRNA gene"))  %>% 
    mutate(Depth_factor = case_when(Depth_m < 10 ~ "<10",
                                    Depth_m < 100 ~ "10-100",
                                    TRUE ~ ">100")) %>% 
    mutate(Depth_factor = factor(Depth_factor, 
                                 levels = c("<10", "10-100", ">100"))) %>% 
    ggplot(aes(rarity, counts))+
    geom_boxplot(position = "dodge", outlier.shape = "cross") + 
    geom_jitter(aes(col = Depth_factor)) + 
    facet_grid(facets = c("Database", "Sequencing_technology")) + 
    theme_bw() + 
    labs(title = "b",
         x = "Abundance classification",
         y = "",
         col = "Depth (m):",
         subtitle = "ASVs classified at phylum level") + 
    theme(panel.grid = element_blank(),
          legend.position = "top",
          strip.background = element_blank(),
          text = element_text(size = 15),
          legend.text = element_text(size = 15),
          panel.grid.major.y = element_line(linetype = "dashed",colour = "grey"))+
    scale_y_log10() + 
    scale_color_manual(values = brewer.pal(n = 9, "Blues")[c(3,6,9)])
  ,
  layout_matrix = rbind(c(1,2),c(1,2))
)

# Supplementary Figure S6
gridExtra::grid.arrange(
  diversity_count_rarefied %>% 
    mutate(Sequencing_technology = 
             ifelse(Sequencing_technology == "Illumina",
                    "V4-V5 16S rRNA gene",
                    "full-length 16S rRNA gene")) %>%  
    filter(Taxonomic_level == "Phylum", n_classified_ASVs > 0) %>% 
    left_join(mosj_env, by = c("Sample" = "NGS_code")) %>% 
    mutate(Station = factor(Station, levels = c("KB3", "KB6", "KB0", "V12", "V6", "HG-IV"))) %>% 
    mutate(Depth_factor = case_when(Depth_m < 10 ~ "<10",
                                    Depth_m < 100 ~ "10-100",
                                    TRUE ~ ">100")) %>% 
    mutate(Depth_factor = factor(Depth_factor, 
                                 levels = c("<10", "10-100", ">100"))) %>% 
    ggplot(aes(Station, n_classified_ASVs))+
    geom_point(aes(col = Depth_factor)) +
    theme_bw()+ 
    theme(panel.grid = element_blank(),
          legend.position = "top",
          strip.background = element_blank(),
          text = element_text(size = 15),
          legend.text = element_text(size = 15),
          panel.grid.major.y = element_line(linetype = "dashed",colour = "grey"),
          axis.text.x = element_text(size = 10)) + 
    facet_grid(facets = c("Database", "Sequencing_technology")) +
    scale_y_log10() + 
    labs(y = "Number of ASVs (Log10)",
         col = "Depth (m):", 
         title = "a",
         subtitle = "ASVs classified at phylum level\nrarefied") + 
    scale_color_manual(values = brewer.pal(n = 9, "Blues")[c(3,6,9)]),
  
  ## Rare biosphere section
  species_cleaned_rarefied %>% 
    filter(!is.na(Taxonomy),
           Abundance > 1) %>% 
    group_by(Sample, Database, Sequencing_technology, Taxonomic_level, Depth_m) %>% 
    mutate(RelativeAbundance = Abundance * 100/ sum(Abundance)) %>% 
    mutate(rarity = ifelse(RelativeAbundance >= 0.9, "Abundant", "Rare")) %>% 
    ungroup() %>% 
    select(Sequencing_technology, Database, Sample, Taxonomic_level, rarity, Depth_m) %>% 
    group_by(Sequencing_technology, Database, rarity, Sample, Taxonomic_level, Depth_m) %>% 
    summarise(counts = n()) %>% 
    group_by(Sequencing_technology, Database, rarity, Taxonomic_level, Depth_m) %>% 
    filter(Taxonomic_level %in% c("Phylum")) %>% 
    mutate(Sequencing_technology = 
             ifelse(Sequencing_technology == "Illumina",
                    "V4-V5 16S rRNA gene",
                    "full-length 16S rRNA gene"))  %>% 
    mutate(Depth_factor = case_when(Depth_m < 10 ~ "<10",
                                    Depth_m < 100 ~ "10-100",
                                    TRUE ~ ">100")) %>% 
    mutate(Depth_factor = factor(Depth_factor, 
                                 levels = c("<10", "10-100", ">100"))) %>% 
    ggplot(aes(rarity, counts))+
    geom_boxplot(position = "dodge", outlier.shape = "cross") + 
    geom_jitter(aes(col = Depth_factor)) + 
    facet_grid(facets = c("Database", "Sequencing_technology")) + 
    theme_bw() + 
    labs(title = "b",
         x = "Abundance classification",
         y = "",
         col = "Depth (m):",
         subtitle = "ASVs classified at Phylum level\nrarefied") + 
    theme(panel.grid = element_blank(),
          legend.position = "top",
          strip.background = element_blank(),
          text = element_text(size = 15),
          legend.text = element_text(size = 15),
          panel.grid.major.y = element_line(linetype = "dashed",colour = "grey"))+
    scale_y_log10() + 
    scale_color_manual(values = brewer.pal(n = 9, "Blues")[c(3,6,9)])
  ,
  layout_matrix = rbind(c(1,2),c(1,2))
)


## Alpha diversity along depth
## summary metrics for alpha diversity across depth
# not rarefied
summary_metrics <- 
  species_cleaned %>% 
  group_by(Sample, Taxonomic_level, Sequencing_technology, Database) %>%
  filter(!is.na(Taxonomy)) %>% 
  nest() %>% 
  mutate(n_ASV = purrr::map(.x = data, 
                            .f = ~vegan::specnumber(.x$Abundance)),
         Shannon = purrr::map(.x = data, 
                              .f = ~vegan::diversity(.x$Abundance)),
         Simpson = purrr::map(.x = data, 
                              .f = ~vegan::diversity(.x$Abundance, index = "simpson"))) %>% 
  unnest(cols = c("n_ASV" , "Shannon", "Simpson"))

# rarefied
summary_metrics_rarefied <- 
  species_cleaned_rarefied %>% 
  group_by(Sample, Taxonomic_level, Sequencing_technology, Database) %>%
  filter(!is.na(Taxonomy)) %>% 
  nest() %>% 
  mutate(n_ASV = purrr::map(.x = data, 
                            .f = ~vegan::specnumber(.x$Rarefied_abundance)),
         Shannon = purrr::map(.x = data, 
                              .f = ~vegan::diversity(.x$Rarefied_abundance)),
         Simpson = purrr::map(.x = data, 
                              .f = ~vegan::diversity(.x$Rarefied_abundance, index = "simpson"))) %>% 
  unnest(cols = c("n_ASV" , "Shannon", "Simpson"))


## phylum level, not rarefied
gridExtra::grid.arrange(
  summary_metrics %>%
    mutate(log_ASV_richness = log(n_ASV)) %>% 
    pivot_longer(cols = c("log_ASV_richness", "Shannon", "Simpson"),
                 values_to = "Score",
                 names_to = "Alpha_index") %>%
    left_join(mosj_env, by = c("Sample" = "NGS_code")) %>% 
    mutate(Station = factor(Station, 
                            levels = c("KB3", "KB6", "KB0", "V12", "V6", "HG-IV"))) %>%
    mutate(Alpha_index = ifelse(Alpha_index == "log_ASV_richness", 
                                "ASV richness (log10)", Alpha_index)) %>% 
    filter(Taxonomic_level == "Phylum") %>% 
    mutate(Sequencing_technology = ifelse(Sequencing_technology == "Illumina", 
                                          "V4-V5 16S rRNA gene",
                                          "full-length 16S rRNA gene")) %>% 
    # start of ggplot
    ggplot(aes(x = Depth_m, y = Score, fill = Sequencing_technology)) + 
    facet_grid(facets = c("Alpha_index", "Database"), scales = "free_y") + 
    geom_point(shape = 21) +
    geom_smooth(se = FALSE, 
                method = "lm", 
                lty = "dashed", 
                aes(col = Sequencing_technology)) +
    theme_bw()+
    theme(axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "top",
          strip.background = element_blank(),
          text = element_text(size = 15),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          panel.grid.major.y = element_line(linetype = "dashed",colour = "grey"))+
    scale_fill_manual(values = qualitative_colors[c(1:2)]) + 
    scale_color_manual(values = qualitative_colors[c(1:2)]) +
    scale_x_log10() +
    labs(x = "Log10 Depth (m)",
         y = "Score",
         col = "Sequencing of ",
         title = "a",
         subtitle = "Phylum level") +
    guides(fill = "none")
  ,
  summary_metrics %>%
    mutate(log_ASV_richness = log(n_ASV)) %>% 
    pivot_longer(cols = c("log_ASV_richness", "Shannon", "Simpson"),
                 values_to = "Score",
                 names_to = "Alpha_index") %>%
    left_join(mosj_env, by = c("Sample" = "NGS_code")) %>% 
    mutate(Station = factor(Station, 
                            levels = c("KB3", "KB6", "KB0", "V12", "V6", "HG-IV"))) %>%
    mutate(Alpha_index = ifelse(Alpha_index == "log_ASV_richness", 
                                "ASV richness (log10)", Alpha_index)) %>% 
    filter(Taxonomic_level == "Species") %>% 
    mutate(Sequencing_technology = ifelse(Sequencing_technology == "Illumina", 
                                          "V4-V5 16S rRNA gene",
                                          "full-length 16S rRNA gene")) %>% 
    ggplot(aes(x = Depth_m, y = Score, fill = Sequencing_technology)) + 
    facet_grid(facets = c("Alpha_index", "Database"), scales = "free_y") + 
    geom_point(shape = 21) + 
    geom_smooth(se = FALSE, 
                method = "lm", 
                lty = "dashed", 
                aes(col = Sequencing_technology)) +
    theme_bw()+
    theme(axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "top",
          strip.background = element_blank(),
          text = element_text(size = 15),
          legend.text = element_text(size = 10),
          panel.grid.major.y = element_line(linetype = "dashed",colour = "grey"), 
          legend.title = element_text(size = 10))+
    scale_fill_manual(values = qualitative_colors[c(1:2)]) +
    scale_color_manual(values = qualitative_colors[c(1:2)]) +
    scale_x_log10() +
    labs(x = "Log10 Depth (m)",
         y = "Score",
         col = "Sequencing of ",
         title = "b",
         subtitle = "Species level") + 
    guides(fill = "none"),
  ncol=2)


## supplementary figure S4
gridExtra::grid.arrange(
  summary_metrics_rarefied %>%
    mutate(log_ASV_richness = log(n_ASV)) %>% 
    pivot_longer(cols = c("log_ASV_richness", "Shannon", "Simpson"),
                 values_to = "Score",
                 names_to = "Alpha_index") %>%
    left_join(mosj_env, by = c("Sample" = "NGS_code")) %>% 
    mutate(Station = factor(Station, 
                            levels = c("KB3", "KB6", "KB0", "V12", "V6", "HG-IV"))) %>%
    filter(Taxonomic_level == "Phylum") %>% 
    mutate(Sequencing_technology = ifelse(Sequencing_technology == "Illumina", 
                                          "V4-V5 16S rRNA gene",
                                          "full-length 16S rRNA gene")) %>%
    mutate(Alpha_index = ifelse(Alpha_index == "log_ASV_richness", 
                                "ASV richness (log10)", Alpha_index)) %>% 
    ggplot(aes(x = Depth_m, y = Score, fill = Sequencing_technology)) + 
    facet_grid(facets = c("Alpha_index", "Database"), scales = "free_y") + 
    geom_point(shape = 21) +
    geom_smooth(se = FALSE, 
                method = "lm", 
                lty = "dashed", 
                aes(col = Sequencing_technology)) +
    theme_bw()+
    theme(axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "top",
          strip.background = element_blank(),
          text = element_text(size = 15),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          panel.grid.major.y = element_line(linetype = "dashed",colour = "grey"))+
    scale_fill_manual(values = qualitative_colors[c(1:2)]) + 
    scale_color_manual(values = qualitative_colors[c(1:2)]) +
    scale_x_log10() +
    labs(x = "Log10 Depth (m)",
         y = "Score",
         col = "Sequencing of ",
         title = "a",
         subtitle = "Phylum level") + 
    guides(fill = "none")
  ,
  summary_metrics_rarefied %>%
    mutate(log_ASV_richness = log(n_ASV)) %>% 
    pivot_longer(cols = c("log_ASV_richness", "Shannon", "Simpson"),
                 values_to = "Score",
                 names_to = "Alpha_index") %>%
    left_join(mosj_env, by = c("Sample" = "NGS_code")) %>% 
    mutate(Station = factor(Station, 
                            levels = c("KB3", "KB6", "KB0", "V12", "V6", "HG-IV"))) %>%
    filter(Taxonomic_level == "Species") %>% 
    mutate(Sequencing_technology = ifelse(Sequencing_technology == "Illumina", 
                                          "V4-V5 16S rRNA gene",
                                          "full-length 16S rRNA gene")) %>%
    mutate(Alpha_index = ifelse(Alpha_index == "log_ASV_richness", 
                                "ASV richness (log10)", Alpha_index)) %>% 
    ggplot(aes(x = Depth_m, y = Score, fill = Sequencing_technology)) + 
    facet_grid(facets = c("Alpha_index", "Database"), scales = "free_y") + 
    geom_point(shape = 21) + 
    geom_smooth(se = FALSE, 
                method = "lm", 
                lty = "dashed", 
                aes(col = Sequencing_technology)) +
    theme_bw()+
    theme(axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "top",
          strip.background = element_blank(),
          text = element_text(size = 15),
          legend.text = element_text(size = 10),
          panel.grid.major.y = element_line(linetype = "dashed",colour = "grey"), 
          legend.title = element_text(size = 10))+
    scale_fill_manual(values = qualitative_colors[c(1:2)]) +
    scale_color_manual(values = qualitative_colors[c(1:2)]) +
    scale_x_log10() +
    labs(x = "Log10 Depth (m)",
         y = "Score",
         col = "Sequencing of ",
         title = "b",
         subtitle = "Species level") + 
    guides(fill = "none"),
  ncol=2)

###### Beta diversity below ######
ASV_full_length_mosj <- readRDS("./AmpliconProcessing/ASV_full_length_mosj")
#
pacbio_gtdb_samples <- 
  ASV_full_length_mosj %>% 
  select(contains("M19")) %>% 
  select(!all_of(samples_to_remove)) %>% 
  colnames() %>% 
  paste0("_GTDB")

pacbio_silva_samples <- 
  ASV_full_length_mosj %>% 
  select(contains("M19")) %>% 
  select(!all_of(samples_to_remove)) %>% 
  colnames() %>% 
  paste0("_Silva")

#
###

### Replace NAs with zeros

## GTDB 
pac_gtdb_classified_species <- 
  ASV_full_length_mosj %>%
  filter(Database == "GTDB", 
         !is.na(Species)) %>% 
  select(contains("M19"))
#
pac_gtdb_na_species <- 
  ASV_full_length_mosj %>%
  filter(Database == "GTDB", 
         is.na(Species)) %>% 
  select(contains("M19"))

pac_gtdb_na_species[pac_gtdb_na_species > 0] <- 0


pac_gtdb_na_species_zeroed <- 
  rbind(pac_gtdb_classified_species,
        pac_gtdb_na_species)


## Silva 
pac_silva_classified_species <- 
  ASV_full_length_mosj %>%
  filter(Database == "Silva", 
         !is.na(Species)) %>% 
  select(contains("M19"))
#
pac_silva_na_species <- 
  ASV_full_length_mosj %>%
  filter(Database == "Silva", 
         is.na(Species)) %>% 
  select(contains("M19"))

pac_silva_na_species[pac_silva_na_species > 0] <- 0


pac_silva_na_species_zeroed <- 
  rbind(pac_silva_classified_species,
        pac_silva_na_species)

###

colnames(pac_gtdb_na_species_zeroed) <- pacbio_gtdb_samples 
colnames(pac_silva_na_species_zeroed) <- pacbio_silva_samples


##
pacbio_silva_vs_gtdb_nas_zeroed_t <- 
  pac_gtdb_na_species_zeroed %>% 
  cbind(pac_silva_na_species_zeroed) %>% t()

## Remove empty rows
pacbio_silva_vs_gtdb_nas_zeroed_t <- 
  pacbio_silva_vs_gtdb_nas_zeroed_t[-rowSums(pacbio_silva_vs_gtdb_nas_zeroed_t) != 0, ]

###  removed samples:
#M19_70_Silva
#M19_28_Silva
#M19_5_Silva

###

mds_pacbio_silva_vs_gtdb_nas_zeroed_t <- 
  pacbio_silva_vs_gtdb_nas_zeroed_t %>% metaMDS()

beta_plot_parameters <- readRDS("./data/beta_plot_parameters")

env_pac_gtdb <- beta_plot_parameters %>% 
  mutate(Sample_code = paste0(Sample, "_GTDB"),
         Database = "GTDB",
         Database_col = qualitative_colors[2])

env_pac_silva <- beta_plot_parameters %>% 
  mutate(Sample_code = paste0(Sample, "_Silva"),
         Database = "Silva",
         Database_col = qualitative_colors[3])

env_all_pac <- env_pac_gtdb %>% 
  rbind(env_pac_silva)

## remove removed samples from env data
env_all_pac_species <- env_all_pac %>% 
  filter(!Sample_code %in% c("M19_70_Silva", "M19_28_Silva", "M19_5_Silva"))


## plots
row.names(env_all_pac) <- env_all_pac$Sample_code
env_all_pac_species$Station <- as.factor(env_all_pac_species$Station)


## Plot
# Full-length 16S rRNA gene: GTDB vs Silva)
plot(mds_pacbio_silva_vs_gtdb_nas_zeroed_t,
     display = "sites", type="p", 
     main = "Full-length 16S rRNA gene: GTDB vs Silva \nSpecies level")
points(mds_pacbio_silva_vs_gtdb_nas_zeroed_t,
       display = "sites",
       bg = env_all_pac$Database_col,
       pch = 21, 
       col = "grey", 
       cex = 2)
with(env_all_pac_species,
     ordiellipse(mds_pacbio_silva_vs_gtdb_nas_zeroed_t,
                 Database, kind="se", conf=0.95, col="grey"))
with(env_all_pac_species,
     ordihull(mds_pacbio_silva_vs_gtdb_nas_zeroed_t,
              Database, col = "grey"))
with(env_all_pac_species,
     ordispider(mds_pacbio_silva_vs_gtdb_nas_zeroed_t,
                Database, label = TRUE, col= "grey"))

# beta pacbio at phylum level
#
pacbio_gtdb_samples <- 
  ASV_full_length_mosj %>% 
  select(contains("M19")) %>% 
  select(!all_of(samples_to_remove)) %>% 
  colnames() %>% 
  paste0("_GTDB")

pacbio_silva_samples <- 
  ASV_full_length_mosj %>% 
  select(contains("M19")) %>% 
  select(!all_of(samples_to_remove)) %>% 
  colnames() %>% 
  paste0("_Silva")

#
###

### Replace NAs with zeros

## GTDB 
# filter out NAs from species table
pac_gtdb_classified_phylum <- 
  ASV_full_length_mosj %>%
  filter(Database == "GTDB", 
         !is.na(Phylum)) %>% 
  select(contains("M19"))

# get the NAs from species table
pac_gtdb_na_phylum <- 
  ASV_full_length_mosj %>%
  filter(Database == "GTDB", 
         is.na(Phylum)) %>% 
  select(contains("M19"))
# transform NAs to zero
pac_gtdb_na_phylum[pac_gtdb_na_phylum > 0] <- 0

# Get the NAs back, but as zeroes
pac_gtdb_na_phyla_zeroed <- 
  rbind(pac_gtdb_classified_phylum,
        pac_gtdb_na_phylum)


## Silva 
pac_silva_classified_phylum <- 
  ASV_full_length_mosj %>%
  filter(Database == "Silva", 
         !is.na(Phylum)) %>% 
  select(contains("M19"))
#
pac_silva_na_phylum <- 
  ASV_full_length_mosj %>%
  filter(Database == "Silva", 
         is.na(Phylum)) %>% 
  select(contains("M19"))

pac_silva_na_phylum[pac_silva_na_phylum > 0] <- 0


pac_silva_na_phyla_zeroed <- 
  rbind(pac_silva_classified_phylum,
        pac_silva_na_phylum)

###
# change colnames
colnames(pac_gtdb_na_phyla_zeroed) <- pacbio_gtdb_samples 
colnames(pac_silva_na_phyla_zeroed) <- pacbio_silva_samples


# merge gtdb and silva
pacbio_silva_vs_gtdb_nas_zeroed_t_phyla <- 
  pac_gtdb_na_phyla_zeroed %>% 
  cbind(pac_silva_na_phyla_zeroed) %>% t()

## Remove empty rows
pacbio_silva_vs_gtdb_nas_zeroed_t <- 
  pacbio_silva_vs_gtdb_nas_zeroed_t[-rowSums(pacbio_silva_vs_gtdb_nas_zeroed_t) != 0, ]

###

mds_pacbio_silva_vs_gtdb_nas_zeroed_t_phyla <- 
  pacbio_silva_vs_gtdb_nas_zeroed_t_phyla %>% metaMDS()

env_pac_gtdb_phyla <- beta_plot_parameters %>% 
  mutate(Sample_code = paste0(Sample, "_GTDB"),
         Database = "GTDB",
         Database_col = qualitative_colors[2])

env_pac_silva_phyla <- beta_plot_parameters %>% 
  mutate(Sample_code = paste0(Sample, "_Silva"),
         Database = "Silva",
         Database_col = qualitative_colors[3])

env_all_pac_phyla <- env_pac_gtdb_phyla %>% 
  rbind(env_pac_silva_phyla)


## Beta diversity phylum level
# Full-length 16S rRNA gene: GTDB vs Silva
plot(mds_pacbio_silva_vs_gtdb_nas_zeroed_t_phyla,
     display = "sites", type="p", 
     main = "Full-length 16S rRNA gene: GTDB vs Silva \nPhylum level")
points(mds_pacbio_silva_vs_gtdb_nas_zeroed_t_phyla,
       display = "sites",
       bg = env_all_pac$Database_col,
       pch = 21, 
       col = "grey", 
       cex = 2)
with(env_all_pac_phyla,
     ordiellipse(mds_pacbio_silva_vs_gtdb_nas_zeroed_t_phyla,
                 Database, kind="se", conf=0.95, col="grey"))
with(env_all_pac_phyla,
     ordihull(mds_pacbio_silva_vs_gtdb_nas_zeroed_t_phyla,
              Database, col = "grey"))
with(env_all_pac_phyla,
     ordispider(mds_pacbio_silva_vs_gtdb_nas_zeroed_t_phyla,
                Database, label = TRUE, col= "grey"))

#
mosj_2019_illumina_gtdb <- readRDS("./data/mosj_2019_illumina_gtdb")
mosj_2019_illumina_silva <- readRDS("./data/mosj_2019_illumina_silva")

illumina_gtdb_samples <- 
  mosj_2019_illumina_gtdb %>% 
  select(contains("M19")) %>% 
  select(!all_of(samples_to_remove)) %>% 
  colnames() %>% 
  paste0("_GTDB")

illumina_silva_samples <- 
  mosj_2019_illumina_silva %>% 
  select(contains("M19")) %>% 
  select(!all_of(samples_to_remove)) %>% 
  colnames() %>% 
  paste0("_Silva")

#
###

### Replace NAs with zeros

## GTDB 
ilu_gtdb_classified_species <- 
  mosj_2019_illumina_gtdb %>%
  filter(Database == "GTDB", 
         !is.na(Species)) %>% 
  select(contains("M19"))
#
ilu_gtdb_na_species <- 
  mosj_2019_illumina_gtdb %>%
  filter(Database == "GTDB", 
         is.na(Species)) %>% 
  select(contains("M19"))

ilu_gtdb_na_species[ilu_gtdb_na_species > 0] <- 0


ilu_gtdb_na_species_zeroed <- 
  rbind(ilu_gtdb_classified_species,
        ilu_gtdb_na_species)


## Silva 
ilu_silva_classified_species <- 
  mosj_2019_illumina_silva %>%
  filter(Database == "Silva", 
         !is.na(Species)) %>% 
  select(contains("M19"))
#
ilu_silva_na_species <- 
  mosj_2019_illumina_silva %>%
  filter(Database == "Silva", 
         is.na(Species)) %>% 
  select(contains("M19"))

ilu_silva_na_species[ilu_silva_na_species > 0] <- 0


ilu_silva_na_species_zeroed <- 
  rbind(ilu_silva_classified_species,
        ilu_silva_na_species)

###

colnames(ilu_gtdb_na_species_zeroed) <- illumina_gtdb_samples 
colnames(ilu_silva_na_species_zeroed) <- illumina_silva_samples


##
illumina_silva_vs_gtdb_nas_zeroed_t <- 
  ilu_gtdb_na_species_zeroed %>% 
  cbind(ilu_silva_na_species_zeroed) %>% t()

## Remove empty rows
illumina_silva_vs_gtdb_nas_zeroed_t <- 
  illumina_silva_vs_gtdb_nas_zeroed_t[-rowSums(illumina_silva_vs_gtdb_nas_zeroed_t) != 0, ]

###

mds_illumina_silva_vs_gtdb_nas_zeroed_t <- 
  illumina_silva_vs_gtdb_nas_zeroed_t %>% metaMDS()

env_ilu_gtdb <- beta_plot_parameters %>% 
  mutate(Sample_code = paste0(Sample, "_GTDB"),
         Database = "GTDB",
         Database_col = qualitative_colors[2])

env_ilu_silva <- beta_plot_parameters %>% 
  mutate(Sample_code = paste0(Sample, "_Silva"),
         Database = "Silva",
         Database_col = qualitative_colors[3])

env_all_ilu <- env_ilu_gtdb %>% 
  rbind(env_ilu_silva)

## remove removed samples from env data
#env_all_ilu <- env_all_ilu %>% filter(!Sample_code %in% c("M19_70_Silva", "M19_28_Silva", "M19_5_Silva"))


## plots
row.names(env_all_ilu) <- env_all_ilu$Sample_code
env_all_ilu$Station <- as.factor(env_all_ilu$Station)

# Full-length 16S rRNA gene: GTDB vs Silva)
plot(mds_illumina_silva_vs_gtdb_nas_zeroed_t,
     display = "sites", type="p", 
     main = "Full-length 16S rRNA gene: GTDB vs Silva \nSpecies level")
points(mds_illumina_silva_vs_gtdb_nas_zeroed_t,
       display = "sites",
       bg = env_all_ilu$Database_col,
       pch = 21, 
       col = "grey", 
       cex = 2)
with(env_all_ilu,
     ordiellipse(mds_illumina_silva_vs_gtdb_nas_zeroed_t,
                 Database, kind="se", conf=0.95, col="grey"))
with(env_all_ilu,
     ordihull(mds_illumina_silva_vs_gtdb_nas_zeroed_t,
              Database, col = "grey"))
with(env_all_ilu,
     ordispider(mds_illumina_silva_vs_gtdb_nas_zeroed_t,
                Database, label = TRUE, col= "grey"))
# illumina phylum level beta diversity
### Replace NAs with zeros

## GTDB 
ilu_gtdb_classified_phylum <- 
  mosj_2019_illumina_gtdb %>%
  filter(Database == "GTDB", 
         !is.na(Phylum)) %>% 
  select(contains("M19"))
#
ilu_gtdb_na_phylum <- 
  mosj_2019_illumina_gtdb %>%
  filter(Database == "GTDB", 
         is.na(Phylum)) %>% 
  select(contains("M19"))

ilu_gtdb_na_phylum[ilu_gtdb_na_phylum > 0] <- 0


ilu_gtdb_na_phyla_zeroed <- 
  rbind(ilu_gtdb_classified_phylum,
        ilu_gtdb_na_phylum)


## Silva 
ilu_silva_classified_phylum <- 
  mosj_2019_illumina_silva %>%
  filter(Database == "Silva", 
         !is.na(Phylum)) %>% 
  select(contains("M19"))
#
ilu_silva_na_phylum <- 
  mosj_2019_illumina_silva %>%
  filter(Database == "Silva", 
         is.na(Phylum)) %>% 
  select(contains("M19"))

ilu_silva_na_phylum[ilu_silva_na_phylum > 0] <- 0


ilu_silva_na_phyla_zeroed <- 
  rbind(ilu_silva_classified_phylum,
        ilu_silva_na_phylum)

###

colnames(ilu_gtdb_na_phyla_zeroed) <- illumina_gtdb_samples 
colnames(ilu_silva_na_phyla_zeroed) <- illumina_silva_samples


##
illumina_silva_vs_gtdb_nas_zeroed_t_phyla <- 
  ilu_gtdb_na_phyla_zeroed %>% 
  cbind(ilu_silva_na_phyla_zeroed) %>% t()

## Remove empty rows
illumina_silva_vs_gtdb_nas_zeroed_t_phyla <- 
  illumina_silva_vs_gtdb_nas_zeroed_t_phyla[-rowSums(illumina_silva_vs_gtdb_nas_zeroed_t_phyla) != 0, ]

###

mds_illumina_silva_vs_gtdb_nas_zeroed_t_phylum <- 
  illumina_silva_vs_gtdb_nas_zeroed_t_phyla %>% metaMDS()

env_ilu_gtdb <- beta_plot_parameters %>% 
  mutate(Sample_code = paste0(Sample, "_GTDB"),
         Database = "GTDB",
         Database_col = qualitative_colors[2])

env_ilu_silva <- beta_plot_parameters %>% 
  mutate(Sample_code = paste0(Sample, "_Silva"),
         Database = "Silva",
         Database_col = qualitative_colors[3])

env_all_ilu <- env_ilu_gtdb %>% 
  rbind(env_ilu_silva)

## remove removed samples from env data
#env_all_ilu <- env_all_ilu %>% filter(!Sample_code %in% c("M19_70_Silva", "M19_28_Silva", "M19_5_Silva"))


## plots
row.names(env_all_ilu) <- env_all_ilu$Sample_code
env_all_ilu$Station <- as.factor(env_all_ilu$Station)

# Full-length 16S rRNA gene: GTDB vs Silva)
plot(mds_illumina_silva_vs_gtdb_nas_zeroed_t_phylum,
     display = "sites", type="p", 
     main = "Full-length 16S rRNA gene: GTDB vs Silva \nPhylum level")
points(mds_illumina_silva_vs_gtdb_nas_zeroed_t,
       display = "sites",
       bg = env_all_ilu$Database_col,
       pch = 21, 
       col = "grey", 
       cex = 2)
with(env_all_ilu,
     ordiellipse(mds_illumina_silva_vs_gtdb_nas_zeroed_t_phylum,
                 Database, kind="se", conf=0.95, col="grey"))
with(env_all_ilu,
     ordihull(mds_illumina_silva_vs_gtdb_nas_zeroed_t_phylum,
              Database, col = "grey"))
with(env_all_ilu,
     ordispider(mds_illumina_silva_vs_gtdb_nas_zeroed_t_phylum,
                Database, label = TRUE, col= "grey"))

####

### Merge beta diversity plots

par(mfrow = c(2,2))
### Phylum level
## pacbio
plot(mds_pacbio_silva_vs_gtdb_nas_zeroed_t_phyla,
     display = "sites", type="p", 
     main = "Full-length 16S rRNA gene: GTDB vs Silva \nPhylum level")
points(mds_pacbio_silva_vs_gtdb_nas_zeroed_t_phyla,
       display = "sites",
       bg = env_all_pac$Database_col,
       pch = 21, 
       col = "grey", 
       cex = 2)
with(env_all_pac_phyla,
     ordiellipse(mds_pacbio_silva_vs_gtdb_nas_zeroed_t_phyla,
                 Database, kind="se", conf=0.95, col="grey"))
with(env_all_pac_phyla,
     ordihull(mds_pacbio_silva_vs_gtdb_nas_zeroed_t_phyla,
              Database, col = "grey"))
with(env_all_pac_phyla,
     ordispider(mds_pacbio_silva_vs_gtdb_nas_zeroed_t_phyla,
                Database, label = TRUE, col= "grey"))

## illumina
plot(mds_illumina_silva_vs_gtdb_nas_zeroed_t_phylum,
     display = "sites", type="p", 
     main = "Full-length 16S rRNA gene: GTDB vs Silva \nPhylum level")
points(mds_illumina_silva_vs_gtdb_nas_zeroed_t_phylum,
       display = "sites",
       bg = env_all_ilu$Database_col,
       pch = 21, 
       col = "grey", 
       cex = 2)
with(env_all_ilu,
     ordiellipse(mds_illumina_silva_vs_gtdb_nas_zeroed_t_phylum,
                 Database, kind="se", conf=0.95, col="grey"))
with(env_all_ilu,
     ordihull(mds_illumina_silva_vs_gtdb_nas_zeroed_t_phylum,
              Database, col = "grey"))
with(env_all_ilu,
     ordispider(mds_illumina_silva_vs_gtdb_nas_zeroed_t_phylum,
                Database, label = TRUE, col= "grey"))
### Species level
## pacbio
# Full-length 16S rRNA gene: GTDB vs Silva)
plot(mds_pacbio_silva_vs_gtdb_nas_zeroed_t,
     display = "sites", type="p", 
     main = "Full-length 16S rRNA gene: GTDB vs Silva \nSpecies level")
points(mds_pacbio_silva_vs_gtdb_nas_zeroed_t,
       display = "sites",
       bg = env_all_pac$Database_col,
       pch = 21, 
       col = "grey", 
       cex = 2)
with(env_all_pac_species,
     ordiellipse(mds_pacbio_silva_vs_gtdb_nas_zeroed_t,
                 Database, kind="se", conf=0.95, col="grey"))
with(env_all_pac_species,
     ordihull(mds_pacbio_silva_vs_gtdb_nas_zeroed_t,
              Database, col = "grey"))
with(env_all_pac_species,
     ordispider(mds_pacbio_silva_vs_gtdb_nas_zeroed_t,
                Database, label = TRUE, col= "grey"))
## illumina
plot(mds_illumina_silva_vs_gtdb_nas_zeroed_t,
     display = "sites", type="p", 
     main = "Full-length 16S rRNA gene: GTDB vs Silva \nSpecies level")
points(mds_illumina_silva_vs_gtdb_nas_zeroed_t,
       display = "sites",
       bg = env_all_ilu$Database_col,
       pch = 21, 
       col = "grey", 
       cex = 2)
with(env_all_ilu,
     ordiellipse(mds_illumina_silva_vs_gtdb_nas_zeroed_t,
                 Database, kind="se", conf=0.95, col="grey"))
with(env_all_ilu,
     ordihull(mds_illumina_silva_vs_gtdb_nas_zeroed_t,
              Database, col = "grey"))
with(env_all_ilu,
     ordispider(mds_illumina_silva_vs_gtdb_nas_zeroed_t,
                Database, label = TRUE, col= "grey"))

## Permanovas

# full - phylo
adonis2(pacbio_silva_vs_gtdb_nas_zeroed_t_phyla ~ Database, 
        data = env_all_pac_phyla)

# full - species
adonis2(pacbio_silva_vs_gtdb_nas_zeroed_t ~ Database, 
        data = env_all_pac_species)

# short - phylo
adonis2(illumina_silva_vs_gtdb_nas_zeroed_t_phyla ~ Database,
        data = env_all_ilu)

# short - species
adonis2(illumina_silva_vs_gtdb_nas_zeroed_t ~ Database,
        data = env_all_ilu)

####### END #####

