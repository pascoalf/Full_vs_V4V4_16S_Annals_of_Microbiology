library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)

# get read processing summary from DADA2
pacbio_track <- read.csv("./data/tracker_full_mosj")
illumina_track <- read.csv("./data/tracker_mosj_short")

## PacBio part
#
pacbio_track <- 
  pacbio_track %>% 
  as.data.frame() %>%
  dplyr::rename(Sample = X) %>% 
  mutate(Sample = str_remove(Sample, ".hifi_reads.fastq")) %>% 
  filter(Sample %in% selected_samples_mosj_2019$Sample) 

# Summary of sequencing results
pacbio_track <- 
  pacbio_track %>% 
  filter(!Sample %in% samples_to_remove) %>% 
  mutate(Percentage_hq = nonchim*100/ccs)

#
pacbio_track %>% 
  ggplot(aes(Sample, Percentage_hq))+
  geom_col()+
  theme_bw()+
  labs(y = "Percentage of final, high quality reads",
       title = "Percentage of final, high quality reads (PacBio)")

#
pacbio_track %>% 
  pivot_longer(cols = c("ccs", "primers","filtered","denoised","nonchim"),
               names_to = "Step",
               values_to = "Reads") %>%
  mutate(Step = factor(Step, levels = c("ccs", "primers","filtered","denoised","nonchim"))) %>% 
  ggplot(aes(Sample, Reads, fill = Step))+
  geom_col(position = "dodge")+
  theme_bw()+
  scale_fill_manual(values = RColorBrewer::brewer.pal(5, "Set1"))+
  scale_y_log10()+
  labs(y = "Number of reads (Log10 scale)",
       title = "Reads at each processing step (PacBio)")+
  theme(axis.text.x = element_text(angle = 90))

##
(summaryStatisticsPacBio <- 
  pacbio_track %>% 
  summarise(max_raw = max(ccs),
            min_raw = min(ccs),
            mean_raw = round(mean(ccs), 2),
            sd_raw = round(sd(ccs), 2),
            median_raw = median(ccs),
            iqr_raw = IQR(ccs),
            max_hq = max(nonchim),
            min_hq = min(nonchim),
            mean_hq = round(mean(nonchim), 2),
            sd_hq = round(sd(nonchim), 2),
            median_hq = median(nonchim),
            iqr_hq = IQR(nonchim),
            n = n()) %>% 
  mutate(Strategy = "PacBio") %>% 
    pivot_longer(cols = contains("_"),
                 names_to = "Statistic",
                 values_to = "Reads") %>% 
    mutate(Step = ifelse(str_detect(Statistic,"raw"), "ccs","high quality"),
           Statistic = str_remove(Statistic, "_raw"),
           Statistic = str_remove(Statistic, "_hq"))
)

## Illumina part 
#
illumina_track <- illumina_track %>% 
  as.data.frame() %>%
  filter(!X %in% samples_to_remove) %>% 
  rename(Sample = X) %>% 
  mutate(Sample = str_replace(Sample,"-","_"),
         Sample = str_extract(Sample, "M19_\\d+")) %>% 
  filter(Sample %in% selected_samples_mosj_2019$Sample) 

# Summary of sequencing results
illumina_track <- 
  illumina_track %>% 
  mutate(Percentage_hq = nonchim*100/input)
#
illumina_track %>% 
  ggplot(aes(Sample, Percentage_hq))+
  geom_col()+
  theme_bw()+
  labs(y = "Percentage of final, high quality reads",
       title = "Percentage of final, high quality reads (illumina)")
#
illumina_track %>% 
  pivot_longer(cols = names(illumina_track)[2:7],
               names_to = "Step",
               values_to = "Reads") %>%
  mutate(Step = factor(Step, levels = names(illumina_track)[2:7])) %>% 
  ggplot(aes(Sample, Reads, fill = Step))+
  geom_col(position = "dodge")+
  theme_bw()+
  scale_fill_manual(values = RColorBrewer::brewer.pal(6, "Set1"))+
  labs(y = "Number of reads (Log10 scale)",
       title = "Reads at each processing step (illumina)")+
  theme(axis.text.x = element_text(angle = 90))


##
(summaryStatisticsillumina <- 
    illumina_track %>% 
    summarise(max_raw = max(input),
              min_raw = min(input),
              mean_raw = round(mean(input), 2),
              sd_raw = round(sd(input), 2),
              median_raw = median(input),
              iqr_raw = IQR(input),
              max_hq = max(nonchim),
              min_hq = min(nonchim),
              mean_hq = round(mean(nonchim), 2),
              sd_hq = round(sd(nonchim), 2),
              median_hq = median(nonchim),
              iqr_hq = IQR(nonchim),
              n = n()) %>% 
    mutate(Strategy = "illumina") %>% 
    pivot_longer(cols = contains("_"),
                 names_to = "Statistic",
                 values_to = "Reads") %>% 
    mutate(Step = ifelse(str_detect(Statistic,"raw"), "input","high quality"),
           Statistic = str_remove(Statistic, "_raw"),
           Statistic = str_remove(Statistic, "_hq"))
)

# Combine Illumina and PacBio info
(summaryStatisticsCombined <- 
  summaryStatisticsPacBio %>% 
  rbind(summaryStatisticsillumina))
#
write.csv(summaryStatisticsCombined, "data/summaryStatisticsCombined.csv")

