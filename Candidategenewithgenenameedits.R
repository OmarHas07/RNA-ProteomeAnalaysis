#WD

setwd("/Users/omarhasannin/Documents/Research Projects/Chapter 1/RNA-SeqData/LogFPKMGraphs")

#Packages
install.packages(c("ggraph", "RColorBrewer","viridis", "patchwork" ))
library(dplyr)
library(tidyverse)
library(igraph)
library(ggraph)

library(readxl)
library(readr)
library(patchwork)
library(RColorBrewer)
library(viridis)


set.seed(666)
###Files
Exp_table_long <- read_csv("Exp_Table_long.csv")
Metadata <- read_csv("MetadataR.csv")




# Make the PCA input data numeric matrix - trouble with this next line
Exp_table_log_wide <- Exp_table_long %>% 
  select(gene_id, library, logFPKM) %>% 
  pivot_wider(names_from = library, values_from = logFPKM) %>% 
  dplyr::group_by(gene_id)

head(Exp_table_log_wide)

# prcomp() performs PCA for you, given a numeric matrix, which is just the transposed Exp_table_log_wide, but without the gene ID column. as.data.frame(t(summary(my_pca)$importance)) saves the standard deviation and proportion of variance into a data table. In this case, the 1st PC accounts for 43% of the variance in this experiment. The 2nd PC accounts for 10% of the variance.
my_pca <- prcomp(t(Exp_table_log_wide[, -1]))



pc_importance <- as.data.frame(t(summary(my_pca)$importance))
head(pc_importance, 20)

# To make a PCA plot, we will graph the data stored in my_pca$x, which stores the coordinates of each library in PC space. Let's pull that data out and annotate them (with metadata).
PCA_coord <- my_pca$x[, 1:10] %>% 
  as.data.frame() %>% 
  mutate(Sample = row.names(.)) %>% 
  full_join(Metadata %>% 
              select(Sample, Treatment, Timepoint, Group), by = "Sample")
head(PCA_coord)

PCA_coord <- PCA_coord %>%
  mutate(Timepoint = factor(Timepoint, levels = c(
    "T2", 
    "T48",
    "T96",
    "T144"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels= c(
    "Control",
    "tZ",
    "iP",
    "DZ",
    "cZ",
    "tZ7G",
    "iP7G",
    "DZ7G",
    "tZ9G",
    "iP9G",
    "DZ9G",
    "cZ9G"
  )))







module_lines_color_strip <- expand.grid(
  Treatment = unique(Metadata$Treatment),
  Timepoint = unique(Metadata$Timepoint), 
  stringsAsFactors = F
) %>% 
  mutate(order_x = case_when(
    str_detect(Timepoint, "T2") ~ 1,
    str_detect(Timepoint, "T48") ~ 2,
    str_detect(Timepoint, "T96") ~ 3,
    str_detect(Timepoint, "T144") ~ 4
  )) %>% 
  mutate(stage = factor(Timepoint, levels = c(
    "T2",	"T48",	"T96",	"T144"
  ))) %>% 
  mutate(Timepoint = reorder(Timepoint, order_x)) %>% 
  ggplot(aes(x = Timepoint, y = 1)) +
  facet_grid(. ~ Treatment) +
  geom_tile(aes(fill = Timepoint)) +
  scale_fill_manual(values = brewer.pal(4, "BuPu"))  +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    panel.spacing = unit(1, "lines")
  )

### TargetGenes
gene_mapping_df <- data.frame(
  gene_id = c("AT3G26060", "AT3G06050", "AT1G07890", "AT4G35000", "AT1G08830","AT1G20620"),
  gene_name = c("PRXQ", "PRXIIF", "APX1", "APX3", "SOD1","CAT3")
)

####PCA CORD

# Assuming gene_mapping_df has columns "gene_id" and "gene_name"
TF_TPM <- Exp_table_long %>% 
  filter(gene_id %in% c("AT3G26060", "AT3G06050", "AT1G07890", "AT4G35000", "AT1G08830","AT1G20620")) %>% 
  inner_join(PCA_coord, by = c("library"="Sample")) %>% 
  left_join(gene_mapping_df, by = c("gene_id" = "gene_id")) %>% 
  mutate(order_x = case_when(
    str_detect(Timepoint, "T2") ~ 1,
    str_detect(Timepoint, "T48") ~ 2,
    str_detect(Timepoint, "T96") ~ 3,
    str_detect(Timepoint, "T144") ~ 4
  )) %>% 
  mutate(Timepoint = reorder(Timepoint, order_x)) %>% 
  ggplot(aes(x = Timepoint, y = logFPKM)) +
  facet_grid(gene_name ~ Treatment, scales = "free_y", labeller = labeller(gene_name = label_parsed)) +
  geom_point(aes(fill = Treatment), color = "white", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(group = gene_name), 
               fun = mean, alpha = 0.8, size = 1.1, color = "grey20") +
  scale_fill_manual(values = viridis(12)) +
  labs(x = NULL,
       y = "log10(FPKM)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 14)  # Adjust text size if needed
  )

wrap_plots(TF_TPM, module_lines_color_strip, 
           nrow = 2, heights = c(1, 0.05))


ggsave("Results/Candidate_genes_FPKM.png", height = 7, width = 10, bg = "white")
