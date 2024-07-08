# library dependencies
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


setwd("/Users/omarhasannin/Documents/Research Projects/Chapter 1/RNA-SeqData/AverageExpressionSelectedProcesses")
# read in TPM OR FPKM expression matrix
Exp_table <- read_csv("dataset.csv", col_types = cols())
head(Exp_table)
dim(Exp_table)
# 33548; 143

# read in metadata file
Metadata <- read_csv("MetadataR.csv", col_names = TRUE)
head(Metadata)
dim(Metadata)
#  142;  7

# You could do it with bait genes, or not depending on the your project. 

Baits <- read_csv("target.csv") #101, aftre adding NAC, and WRKYs 133
head(Baits)
dim(Baits)


unique(Baits$Process)

# read in bait genes
#Baits <- read_delim("Baitgenes.txt", 
#  delim = "\t", escape_double = FALSE, 
# col_names = FALSE, trim_ws = TRUE)
#head(Baits)



# Setting up the major factors in the experiment 
Metadata %>% 
  dplyr::group_by(Treatment) %>% 
  dplyr::count()
# A tibble: 12 × 12

Metadata %>% 
  dplyr::group_by(Replicate) %>% 
  dplyr::count()
# A tibble: 3 × 2

Metadata %>% 
  dplyr::group_by(Timepoint) %>% 
  dplyr::count()
# A tibble: 35 × 35 x 36 x 36

Metadata%>%
  group_by(Treatment, Timepoint) %>%
  count()

Metadata %>% 
  group_by(Treatment, Timepoint) %>% 
  count()
Metadata %>% 
  group_by(Replicate) %>% 
  count()
# A tibble: 48 × 3
# Groups:   Treatment, Timepoint [48]

## go from wide to tidy format
Exp_table_long <- Exp_table %>% 
  pivot_longer(cols = !gene_id, names_to = "library", values_to = "FPKM") %>% 
  mutate(logFPKM = log10(FPKM + 1)) 

head(Exp_table_long)




# Make the PCA input data numeric matrix - trouble with this next line
Exp_table_log_wide <- Exp_table_long %>% 
  dplyr::select(gene_id, library, logFPKM) %>% 
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
              dplyr::select(Sample, Treatment, Timepoint, Type), by = "Sample")
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





#Now that we have a better handle on the data we will do the co-expression analysis
#Average across the reps
# We start from the long (tidy) table we made earlier. I also pulled the metadata as well to guide the averaging process. 
# By = c("library"="Run) inside full_join() deals with the fact that the library ID is called library 
# in the long table, but Run in the metadata. 
# group_by() followed by summarise(mean = ...) - this takes each gene, tissue, and dev_stage, and computes the mean. 
# The elegance of a tidyverse based workflow is that you do not have to do loops! 
# You let group_by() do the heavy lifting. This could take a moment. This step is doing a lot of mean calculations

Exp_table_long_averaged <- Exp_table_long %>% 
  full_join(PCA_coord %>% 
              dplyr::select(Sample, Treatment, Type, Timepoint), 
            by = c("library"="Sample")) %>% 
  group_by(gene_id,Treatment, Type, Timepoint) %>% 
  summarise(mean.logFPKM= mean(logFPKM)) %>% 
  ungroup()  

head(Exp_table_long_averaged)

# Calculate the z-score once you average across the replicates
# The z-score is the difference from mean over standard deviation and standardizes expression patterns (mean = 1, sd = 1)
# In this step, we are grouping by gene. 
# Tissue-stages with higher expression will have a higher z score and vice versa. 
# Note that this is completely relative to each gene itself. 

Exp_table_long_averaged_z <- Exp_table_long_averaged %>% 
  group_by(gene_id) %>% 
  mutate(z.score = (mean.logFPKM - mean(mean.logFPKM))/sd(mean.logFPKM)) %>% 
  ungroup()

head(Exp_table_long_averaged_z)


# Gene selection
# The next step is correlating each gene to every other gene. However, we have almost 57k genes in this dataset. 
# The number of correlations scales to the square of number of genes. 
# To make things faster and less cumbersome, we can select only the high variance genes. 
# The underlying rationale is if a gene is expressed at a similar level across all samples, it is unlikely that is specifically involved in the biology in a particular stage or tissue.
# There are multiple ways to selecting for high variance genes, and multiple cutoffs. 
# For example, you can calculate the gene-wise variance of logTPM for all genes, and take the upper third. You can only take genes with a certain expression level (say > 5 tpm across all tissues), then take high variance gene. These are arbitrary. You do you.
# Using log(TPM) reduces the bias towards highly expressed genes
# Then I filtered for top 33% high var genes.

#To get the excel file of top 5000 high variance genes
high_var_genes <- Exp_table_long_averaged_z %>% 
  group_by(gene_id) %>% 
  summarise(var = var(mean.logFPKM)) %>% 
  ungroup() %>% 
  filter(var > quantile(var, 0.667))

head(high_var_genes)
dim(high_var_genes)
# 11172     2


# The above chunk just listed the high var genes, now we need to filter those out in the long table that contains the z-scores.
# For the sake of this example, let's just take top 5000 genes with highest var as a quick exercise
high_var_genes5000 <- high_var_genes %>% 
  slice_max(order_by = var, n = 10000) 

head(high_var_genes5000)

## Check if the bait genes exist in the gene name column 

high_var_genes5000 %>% 
  filter(gene_id %in% Baits$gene_id) %>% 
  nrow() ##124


## IT isn't clearly if the bait genes are in the top 5000 but You could look at the file "high_var_genes5000" and check by hand

# The %in% operator filters geneIDss that are present in high_var_genes5000$geneIDs, thus retaining only high var genes.
Exp_table_long_averaged_z_high_var <- Exp_table_long_averaged_z %>% 
  filter(gene_id %in% high_var_genes5000$gene_id)

head(Exp_table_long_averaged_z_high_var)


####

Data <- Exp_table_long_averaged_z_high_var %>%  
  filter(gene_id %in% Baits$gene_id)

data_fun <- Data %>% 
  inner_join(Baits, by = "gene_id")

process_mean_z <- data_fun %>% 
  group_by(Process, Timepoint, Treatment, Type) %>% 
  summarise(mean.z = mean(z.score)) %>% 
  ungroup()

head(process_mean_z)

process_mean_z

max(process_mean_z$mean.z)
min(process_mean_z$mean.z)

##### LIne ploting 


module_line_plot <- data_fun %>% 
  mutate(order_x = case_when(
    str_detect(Timepoint, "T2") ~ 1,
    str_detect(Timepoint, "T48") ~ 2,
    str_detect(Timepoint, "T96") ~ 3,
    str_detect(Timepoint, "T144") ~ 4
    
  )) %>% 
  mutate(Timepoint = reorder(Timepoint, order_x)) %>% 
  filter(#Process == "CCGs"| 
    Process == "ROS" | 
      Process == "CK-associated" |
      Process == "Senescence" |
    Process == "LHC" ) %>% 
  ggplot(aes(x = Timepoint, y = z.score)) +
  #do it by treatment
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  facet_grid(Process ~ Treatment) +
  geom_line(aes(group = gene_id), alpha = 0.3, color = "grey70") +
  geom_line(
    data = process_mean_z %>% 
      filter(
              # Process ==   "CCGs"|
       
        Process == "ROS" | 
          Process == "CK-associated" |
          Process == "Senescence" |
          Process == "LHC" 
           ) %>% 
      mutate(order_x = case_when(
        str_detect(Timepoint, "T2") ~ 1,
        str_detect(Timepoint, "T48") ~ 2,
        str_detect(Timepoint, "T96") ~ 3,
        str_detect(Timepoint, "T144") ~ 4
      )) %>% 
      mutate(Timepoint = reorder(Timepoint, order_x)),
    aes(y = mean.z, group = Process), 
    size = 1.1, alpha = 0.8
  ) +
  labs(x = NULL,
       y = "Mean z-scores") +
  theme_classic() +  # Starting from a classic theme
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    panel.spacing = unit(1, "line"),
    panel.grid.major = element_line(color = "gray80", size = 0.5), # Add major grid lines
    panel.grid.minor = element_line(color = "gray90", size = 0.25) # Add minor grid lines
  )

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
  theme_void() +  # Originally removes most theme elements
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_line(color = "gray80", size = 0.5), # Assuming a different base theme
    panel.grid.minor = element_line(color = "gray90", size = 0.25)
  )

wrap_plots(module_line_plot, module_lines_color_strip,
           nrow = 2, heights = c(1, 0.08))

ggsave("Results/CK-assoc.-LHC-ROS-Senescence.png", height = 10, width = 16, dpi = 300, bg = "white")

