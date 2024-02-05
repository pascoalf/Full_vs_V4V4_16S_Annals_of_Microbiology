# Rarefaction curves PacBio vs Illumina
library(tidyverse)
library(vegan)
library(ggplot2)
#
set.seed("1234")
par(mfrow = c(1,2))
## Make rarefaction curves for selected samples
mosj_2019_illumina_silva %>%
  select(contains("M19")) %>% 
  t() %>% 
  rarecurve(step = 10, main = "a \nV4-V5 16S rRNA gene", ylab = "ASVs")
#
colnames(ASV_full_length_mosj) <- colnames(ASV_full_length_mosj) %>% 
  str_remove(".hifi_reads.fastq")
ASV_full_length_mosj %>% 
  select(contains("M19")) %>% 
  t() %>% 
  rarecurve(step = 10, main = "b \nfull-length 16S rRNA gene", ylab = "ASVs")
