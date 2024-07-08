## WD
setwd("/Users/omarhasannin/Documents/Research Projects/Chapter 1/Proteome data /Co-expression")
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
library(limma)


set.seed(666)


# read in TPM OR FPKM expression matrix
Exp_table <- read_csv("abundancesfourrepsallTs.csv", col_types = cols())
head(Exp_table)
dim(Exp_table)
# 4245;, 193

# Remove anything after the dot in the 'GeneID' column
Exp_table <- Exp_table %>%
  mutate(GeneID = sub("\\..*$", "", GeneID))



Baits <- read_xlsx("Genes of interests.xlsx", col_names = TRUE)
colnames(Baits) <- c("GeneID", "gene_name", "process")
head(Baits)


# read in metadata file
Metadata <- read_csv("MetadataR.csv", col_names = TRUE)
head(Metadata)
dim(Metadata)
#  192;  3


# Setting up the major factors in the experiment 
Metadata %>% 
  dplyr::group_by(Treatment) %>% 
  dplyr::count()
# A tibble: 12 × 12


Metadata %>% 
  dplyr::group_by(Timepoint) %>% 
  dplyr::count()
# A tibble: 4XX2

Metadata%>%
  group_by(Treatment, Timepoint) %>%
  count()

colnames(Exp_table)
## go from wide to tidy format
Exp_table_long <- Exp_table %>% 
  pivot_longer(cols = !GeneID, names_to = "library", values_to = "Abundance") %>% 
  mutate(logAbundance = log10(Abundance + 1)) 

head(Exp_table_long)
## Check if it's originally numberic
is.numeric(Exp_table_long$logAbundance) ## True

# Make the PCA input data numeric matrix - trouble with this next line ## it worked just created uniuqe identified row for each librabry first 
## before using pivot wider. 
Exp_table_log_wide <- Exp_table_long %>% 
  select(GeneID, library, logAbundance) %>%  
  group_by(library) %>%  
  mutate(row= row_number()) %>% 
  pivot_wider(names_from = library, values_from = logAbundance) %>% 
  select(-row) %>%  
  dplyr::group_by(GeneID)


#### 

head(Exp_table_log_wide)

# prcomp() performs PCA for you, given a numeric matrix, which is just the transposed Exp_table_log_wide, but without the gene ID column. as.data.frame(t(summary(my_pca)$importance)) saves the standard deviation and proportion of variance into a data table. In this case, the 1st PC accounts for 43% of the variance in this experiment. The 2nd PC accounts for 10% of the variance.
my_pca <- prcomp(t(Exp_table_log_wide[, -1]))

pc_importance <- as.data.frame(t(summary(my_pca)$importance))
head(pc_importance, 20)
####



# To make a PCA plot, we will graph the data stored in my_pca$x, which stores the coordinates of each library in PC space. Let's pull that data out and annotate them (with metadata).
PCA_coord <- my_pca$x[, 1:10] %>% 
  as.data.frame() %>% 
  mutate(Sample = row.names(.)) %>% 
  full_join(Metadata %>% 
              select(Sample, Treatment, Timepoint, Type, Group), by = "Sample")
head(PCA_coord)

PCA_coord$Group

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
    "tZ9G",
    "iP7G",
    "iP9G",
    "DZ7G",
    "DZ9G",
    "cZ9G"
  ))) %>%  
  mutate(Group = factor(Group, levels=c(
    "Control", 
    "BaseForm", 
    "N7-Glucoside", 
    "N9-Glucoside"
  )))



# defining factors - stopped here for PCA plots, need to edit ###############
#Then plotting th output from the PCA analysis with ggplot this is PC1 versus PC2 by replicate
PCA_by_Timepoint <- PCA_coord %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Group, shape = Timepoint), size = 3, alpha = 0.8) +  # Map color to Type and shape to Timepoint
  scale_shape_manual(values = c(T2 = 17, T48 = 15, T96 = 19, T144 = 9)) +  # Assign specific shapes to each time point
  scale_color_brewer(palette = "Dark2") +  # Use a Brewer color palette for Type
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),  
       color = "Group", 
       shape = "Timepoint") +  
  theme_bw() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    legend.position = "right"
  )

# Print the plot as output to view
PCA_by_Timepoint

wrap_plots(PCA_by_Timepoint)
ggsave("PCA_by_Group&TimePoint.png", height = 4, width = 8.5, bg = "white")



### PCA by treatment
PCA_by_Treatment <- PCA_coord %>% 
  ggplot(aes(x = PC2, y = PC3)) +
  geom_point(aes(fill = Treatment), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  labs(x = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC3 (", pc_importance[3, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = "Treatment") +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )


#Print the plot as output to view
PCA_by_Treatment

wrap_plots(PCA_by_Treatment)
ggsave("PCAbyTreatment.png", height = 4.5, width = 8.5, bg = "white")

### PCA by Group
PCA_by_Group <- PCA_coord %>% 
  ggplot(aes(x = PC2, y = PC3)) +
  geom_point(aes(fill = Group), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  labs(x = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC3 (", pc_importance[3, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = "Group") +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )


#Print the plot as output to view
PCA_by_Group
wrap_plots(PCA_by_Group)
ggsave("PCAbyGroup.png", height = 4.5, width = 8.5, bg = "white")

#Print the plot as output to view
PCA_by_Timepoint2

# This is PC1 versus PC2 colored by tissue type
PCA_by_timepoint2 <- PCA_coord %>% filter(Timepoint == "T2") %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Treatment), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  #scale_fill_manual(values = brewer.pal(11, "Set3")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = "Treatment after 2 Hours") +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )
PCA_by_timepoint2

#saving the plot



wrap_plots(PCA_by_timepoint2)
ggsave("../Results/PCA_by_Different timepoint seperate.png", height = 4, width = 8.5, bg = "white")

PCA_by_timepoint48 <- PCA_coord %>% filter(Timepoint == "T48") %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Treatment), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  #scale_fill_manual(values = brewer.pal(11, "Set3")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = "Treatment after 48 Hours") +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )
PCA_by_timepoint48

#saving the plot

wrap_plots(PCA_by_timepoint48)
ggsave("../Results/PCA_by_Different timepoint seperate.png", height = 4, width = 8.5, bg = "white")

PCA_by_timepoint96 <- PCA_coord %>% filter(Timepoint == "T96") %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Treatment), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  #scale_fill_manual(values = brewer.pal(11, "Set3")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = "Treatment after 96 Hours") +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )
PCA_by_timepoint96

#saving the plot

wrap_plots(PCA_by_timepoint96)
ggsave("../Results/PCA_by_Different timepoint seperate.png", height = 4, width = 8.5, bg = "white")

PCA_by_timepoint144 <- PCA_coord %>% filter(Timepoint == "T144") %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Treatment), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  #scale_fill_manual(values = brewer.pal(11, "Set3")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = "Treatment after 144 Hours") +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )
PCA_by_timepoint144

wrap_plots(PCA_by_timepoint2, PCA_by_timepoint48,PCA_by_timepoint96,PCA_by_timepoint144)
ggsave("../Results/PCA_by_Different timepoint seperate.png", height = 8, width = 12, bg = "white")


###################
###################



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
              select(Sample, Treatment, Type, Timepoint), 
            by = c("library"="Sample")) %>% 
  group_by(GeneID,Treatment, Type, Timepoint) %>% 
  summarise(mean.logAbundance= mean(logAbundance)) %>% 
  ungroup()  



head(Exp_table_long_averaged)


# Calculate the z-score once you average across the replicates
# The z-score is the difference from mean over standard deviation and standardizes expression patterns (mean = 1, sd = 1)
# In this step, we are grouping by gene. 
# Tissue-stages with higher expression will have a higher z score and vice versa. 
# Note that this is completely relative to each gene itself. 

Exp_table_long_averaged_z <- Exp_table_long_averaged %>% 
  group_by(GeneID) %>% 
  mutate(z.score = (mean.logAbundance - mean(mean.logAbundance))/sd(mean.logAbundance)) %>% 
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
  group_by(GeneID) %>% 
  summarise(var = var(mean.logAbundance)) %>% 
  ungroup() %>% 
  filter(var > quantile(var, 0.000000001))## it was 0.667
head(high_var_genes)
dim(high_var_genes) ##2119 2
# 4236     2
### with 0.000001 genes are 4237
##0.2 is 3392
### I do not think that I need to filter based on the highly variable genes because there
### already very small number of IDs in the dataset

## Check if the bait genes exist in the gene name column 

high_var_genes5000 %>% 
  filter(GeneID %in% Baits$GeneID) %>% 
  nrow() ##12

BaitPresent <- Baits %>%  
  filter(GeneID %in% high_var_genes5000$GeneID)

# The above chunk just listed the high var genes, now we need to filter those out in the long table that contains the z-scores.
# For the sake of this example, let's just take top 5000 genes with highest var as a quick exercise
high_var_genes5000 <- high_var_genes %>% 
  slice_max(order_by = var, n = 4237) 

head(high_var_genes5000)
nrow(high_var_genes5000)### 1412 ### new 4237
## Check if the bait genes exist in the gene name column 

#high_var_genes5000 %>% 
#  filter(Accession %in% Baits$gene_id) %>% 
 # nrow() ##124 out of 141. 

# The %in% operator filters geneIDss that are present in high_var_genes5000$geneIDs, thus retaining only high var genes.
Exp_table_long_averaged_z_high_var <- Exp_table_long_averaged_z %>% 
  filter(GeneID %in% high_var_genes5000$GeneID)
head(Exp_table_long_averaged_z_high_var)

Exp_table_long_averaged_z_high_var %>% 
  group_by(GeneID) %>% 
  count() %>% 
  nrow() #4236 ### 4237

Exp_table_long_averaged_z_high_var_count <- ncol(Exp_table_long_averaged_z_high_var) - 1
Exp_table_long_averaged_z_high_var_count # this is solution to above code - 6 rows
#5
all_var_and_ranks <- Exp_table_long_averaged_z %>% 
  group_by(GeneID) %>% 
  summarise(var = var(mean.logAbundance)) %>% 
  ungroup() %>% 
  mutate(rank = rank(var, ties.method = "average")) 


bait_var <- all_var_and_ranks %>% 
  mutate(Gene_ID2 = str_sub(GeneID, start = 1, end = 19)) %>% 
  filter(Gene_ID2 %in% Baits$GeneID) %>% 
  group_by(Gene_ID2) %>% 
  slice_max(n = 1, order_by = var)

bait_var
#bait_var ## ### 41 out of 141

all_var_and_ranks %>% 
  ggplot(aes(x = var, y = rank)) +
  geom_rect( 
    xmax = max(high_var_genes5000$var), 
    xmin = min(high_var_genes5000$var),
    ymax = nrow(all_var_and_ranks),
    ymin = nrow(all_var_and_ranks) - 5000,
    fill = "dodgerblue2", alpha = 0.2
  ) +
  geom_line(size = 1.1) +
  geom_hline(
    data = bait_var, aes(yintercept = rank),
    color = "tomato1", size = 0.8, alpha = 0.5
  ) +
  geom_vline(
    data = bait_var, aes(xintercept = var), 
    color = "tomato1", size = 0.8, alpha = 0.5
  ) + 
  labs(y = "rank",
       x = "var(log10(TPM))",
       caption = "Blue box = top 10000 high var genes.\nRed lines = bait genes.") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    plot.caption = element_text(hjust = 0)
  )




# Gene-wise correlation
# We will use the cor() function in R. 
# But the cor() only take vector or matrix as input, so we need to go from long to wide again.

z_score_wide <- Exp_table_long_averaged_z_high_var %>% 
  select(GeneID, Type, z.score) %>% 
  pivot_wider(names_from = Type, values_from = z.score) %>% 
  as.data.frame()

row.names(z_score_wide) <- z_score_wide$GeneID
head(z_score_wide)

#### Exporting z-score file 

write_csv(z_score_wide, "z_score_proteome2024.csv")
# The GroupName column contains info for both treatment,cultivar and tissue, which we can recall using the metadata. 
# After long to wide transformation, the GroupName column now becomes the column name of this wide table. 
# Then we produce the correlation matrix. 
# The underlying math here is R takes each column of a matrix and correlates it to every other columns. 
# To get this to work on our wide table, we remove the geneIDs column, transpose it, and feed it into cor().

cor_matrix <- cor(t(z_score_wide[,-1]))
dim(cor_matrix)
# 1414 1414
##4237 4237

# Edge selection
# This is building up to a network analysis, where each gene is node, and each correlation is an edge. I have two ways to do this.
# 1. t distribution approximation
# 2. Empirical determination using rank distribution

# For the t distribution you need to calculate the number of tissue, cultivar and treatment combos
number_of_timepoint_tissue_trt <- ncol(z_score_wide) - 1
number_of_timepoint_tissue_trt
## [1] 48
## OR you can find it this way
PCA_coord %>% 
  dplyr::group_by(Timepoint, Treatment) %>% 
  dplyr::count() %>% 
  nrow()
## [1] 48

# Since the correlation matrix is symmetrical along the diagonal we can eliminate redundant data.
cor_matrix_upper_tri <- cor_matrix
cor_matrix_upper_tri[lower.tri(cor_matrix_upper_tri)] <- NA


# Now we can compute a t statistic from r and compute a p value using the t distribution. 
# This chunk converts the correlation matrix into a data table. 
# Then it goes from wide to long using pivot_longer(). 
# After that, everything is normal dyplr verbs, such as mutate() and filter(). 
# P values are computed using the t distribution. 
# Depending on the sign of t, the upper of lower tail probability is taken. 
# Finally, the p values are adjusted for multiple comparisons using FDR. 
# This step can take a while. Turning a large wide table to a long table always takes a while. Your computer may not have enough memory to run this step if you put in many genes. In this case we only used 5000 genes, so no problem.

edge_table <- cor_matrix_upper_tri %>% 
  as.data.frame() %>% 
  mutate(from = row.names(cor_matrix)) %>% 
  pivot_longer(cols = !from, names_to = "to", values_to = "r") %>% 
  filter(is.na(r) == F) %>% 
  filter(from != to) %>% 
  mutate(t = r*sqrt((number_of_timepoint_tissue_trt-2)/(1-r^2))) %>% 
  mutate(p.value = case_when(
    t > 0 ~ pt(t, df = number_of_timepoint_tissue_trt-2, lower.tail = F),
    t <=0 ~ pt(t, df = number_of_timepoint_tissue_trt-2, lower.tail = T)
  )) %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr")) 
head(edge_table)

# You can look at various adjusted p value cutoffs and the corresponding r value before proceeding. 
# Let's say we just look at positively correlated genes.
edge_table %>% 
  filter(r > 0) %>% 
  filter(FDR < 0.05) %>% 
  slice_min(order_by = abs(r), n = 10) 

## A tibble: 10 × 6


edge_table %>% 
  filter(r > 0) %>% 
  filter(FDR < 0.01) %>% 
  slice_min(order_by = abs(r), n = 10)

## 10X6



edge_table %>% 
  slice_sample(n = 5000) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(color = "white", bins = 100) +
  geom_vline(xintercept = 0.8, color = "tomato1", size = 1.2) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )


ggsave("r_histogram.png", height = 3.5, width = 5, bg = "white")

edge_table_select <- edge_table %>% 
  filter(r>0.7)
  #filter(FDR < 0.05)
dim(edge_table_select)
## 189703 6
#FDR <0.01 218,6546
#r > 0.7 154,320


Bait_cor_by_Treatment <- z_score_wide %>% 
  filter(GeneID == "AT4G29740" |
           GeneID == "AT5G05860") %>% 
  select(-GeneID) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(Type = row.names(.)) %>% 
  inner_join(PCA_coord, by = "Type") %>% 
  ggplot(aes(x = AT4G29740,
             y = AT5G05860)) +
  geom_point(aes(fill = Treatment), color = "grey20", 
             size = 2, alpha = 0.8, shape = 21) +
  scale_fill_manual(values = viridis(12, option = "D")) +
  labs(x = "CKX4 z score",
       y = "UGT76C2 z score") + 
  theme_classic() +
  theme(
    legend.position = c(0.2, 0.8),
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )
Bait_cor_by_Treatment

### Baits across time points

Bait_cor_by_Time <- z_score_wide %>% 
  filter(GeneID == "AT4G29740" |
           GeneID == "AT5G05860") %>% 
  select(-GeneID) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(Type = row.names(.)) %>% 
  inner_join(PCA_coord, by = "Type") %>% 
  ggplot(aes(x = AT4G29740,
             y = AT5G05860)) +
  geom_point(aes(fill = Timepoint), color = "grey20", 
             size = 2, alpha = 0.8, shape = 21) +
  scale_fill_manual(values = viridis(4, option = "D")) +
  labs(x = "CKX4 z score",
       y = "UGT76C2 z score") + 
  theme_classic() +
  theme(
    legend.position = c(0.2, 0.8),
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )
Bait_cor_by_Time


wrap_plots(Bait_cor_by_Treatment, Bait_cor_by_Time, ncol = 2)
ggsave("CKX4&UGT76C2Correlation.png", height = 10, width = 12, bg = "white")

#1472 2
# Module detection
# We will be using the Leiden algorithm to detect module, which is a graph based clustering method
# We will be using igraph to do some of the downstream analyses. 
# While you can get Leiden as a standalone package, Leiden is also part of the igraph package. 
# The first thing to do is producing a graph object, also known as a network object.

# To make a graph object, you need a edge table. We already made that, which is edge_table_select, a edge table that we filtered based on some kind of r cutoff. 
# Optionally, we can also provide a node table, which contains information about all the notes present in this network. We can make that.
# We need to two things.
# Non-redundant gene IDs from the edge table.
# Functional annotation, which I downloaded.
# Functional annotation
#funct_anno <- read_delim("DM_1-3_516_R44_potato.v6.1.func_annot_transcripts.txt", 
#delim = "\t", escape_double = FALSE, 
#  trim_ws = TRUE)


funct_anno <- read_csv("FunctionalAnnotation.csv", col_names = F)

head(funct_anno)

node_table <- data.frame(
  geneIDs = c(edge_table_select$from, edge_table_select$to) %>%
    unique()) %>% 
  left_join(funct_anno, by = c("geneIDs"="X1")) %>% 
  rename(functional_annotation = X2) %>% 
  unique()
head(node_table)
dim(node_table)
duplicated_rows <- duplicated(node_table$geneIDs)
node_table <- node_table[!duplicated_rows, ]
dim(node_table)
# 3176 2
### 4233 2
write_csv(node_table, "func_anno_node.csv")


##THIS CODE IS WORKING 
##the error that you was getting before was because you need 
##to reference a unique node data in vertices. 
my_network <- graph_from_data_frame(
  edge_table_select,
  vertices = node_table,
  directed = F
)

# Graph-based clustering
# cluster_leiden() runs the Leiden algorithm for you. resolution_parameter controls how many clusters you will get. 
# The larger it is, the more clusters. You can play around with the resolution and see what you get. 

#I made the resolution parameter to be 1.5, that generated 16 different modules at the end, based 
# on you dataset you might need to change that number. 
### Visualize the network
plot.igraph(my_network, vertex.size = 3, vertex.label = NA)
######
### Tesing another alog

modules <- cluster_leiden(my_network, resolution_parameter = 1.5, 
                          objective_function = "modularity", )

### Exporting My Network Modules 
write_csv(my_network_modules, "ModulesadjustedFDR<0.01.csv")

# use heuristics to find optimal number of modules
# Here I wrote a function to detect module, pull out number of modules that have >= 5 genes, and count number of genes contained in modules that have >= 5 genes. All in one function.

optimize_resolution <- function(network, resolution){
  modules = network %>% 
    cluster_leiden(resolution_parameter = resolution,
                   objective_function = "modularity")
  
  parsed_modules = data.frame(
    geneIDs = names(membership(modules)),
    module = as.vector(membership(modules)) 
  )
  
  num_module_5 = parsed_modules %>% 
    group_by(module) %>% 
    count() %>% 
    arrange(-n) %>% 
    filter(n >= 10) %>% 
    nrow() %>% 
    as.numeric()
  
  num_genes_contained = parsed_modules %>% 
    group_by(module) %>% 
    count() %>% 
    arrange(-n) %>% 
    filter(n >= 10) %>% 
    ungroup() %>% 
    summarise(sum = sum(n)) %>% 
    as.numeric()
  
  c(num_module_5, num_genes_contained)
  
}

# Then I can test a list of resolutions in this function. Let's test a range of resolution from 0.25 to 5, in steps of 0.25.

optimization_results <- purrr::map_dfc(
  .x = seq(from = 0.01, to = 3, by = 0.2),
  .f = optimize_resolution, 
  network = my_network
) %>% 
  t() %>% 
  cbind(
    resolution = seq(from = 0.01, to = 3, by = 0.2)
  ) %>% 
  as.data.frame() %>% 
  rename(num_module = V1,
         num_contained_gene = V2)
head(optimization_results)

# We have the results organized into one tidy data table. We can graph it.
Optimize_num_module <- optimization_results %>% 
  ggplot(aes(x = resolution, y = num_module)) +
  geom_line(size = 1.1, alpha = 0.8, color = "dodgerblue2") +
  geom_point(size = 3, alpha = 0.7) +
  geom_vline(xintercept = 2, size = 1, linetype = 4) +
  labs(x = "resolution parameter",
       y = "num. modules\nw/ >=10 genes") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

Optimize_num_gene <- optimization_results %>% 
  ggplot(aes(x = resolution, y = num_contained_gene)) +
  geom_line(size = 1.1, alpha = 0.8, color = "violetred2") +
  geom_point(size = 3, alpha = 0.7) +
  geom_vline(xintercept = 2, size = 1, linetype = 4) +
  labs(x = "resolution parameter",
       y = "num. genes in\nmodules w/ >=10 genes") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

wrap_plots(Optimize_num_module, Optimize_num_gene, nrow = 2)
ggsave("Num.Modules.png", height = 4, width = 8.5, bg = "white")



# Let's say we move on with module detection using a resolution of 2. Next, we need to link the module membership to the gene IDs.

my_network_modules <- data.frame(
  geneIDs = names(membership(modules)),
  module = as.vector(membership(modules)) 
) %>% 
  inner_join(node_table, by = "geneIDs")

my_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 10) 
##9
##13
my_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 10) %>% 
  ungroup() %>% 
  summarise(sum = sum(n))
##3121
##3855
# Moving forward we will only use modules that have 5 or more genes.
module_5 <- my_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 10)

### edit #####
## each module at least contains 10 IDs

my_network_modules <- my_network_modules %>% 
  filter(module %in% module_5$module)
head(my_network_modules)

my_network_modules %>% 
  count(my_network_modules$module) %>% 
  unique()



write_csv(my_network_modules, "9ModulesMayth12024.csv")
# Module treatment correspondence
# the essence of this workflow is simple, so we will use a simple method: peak expression.

#unify columns names first. 

names(Exp_table_long_averaged_z_high_var)
names(Exp_table_long_averaged_z_high_var)[1] <- "geneIDs"


Exp_table_long_averaged_z_high_var_modules <- Exp_table_long_averaged_z_high_var %>% 
  inner_join(my_network_modules, by = "geneIDs")

head(Exp_table_long_averaged_z_high_var_modules)

# Now we can produce summary statistics for each cluster and look at their expression pattern using mean.
modules_mean_z <- Exp_table_long_averaged_z_high_var_modules %>% 
  group_by(module, Timepoint, Treatment, Type) %>% 
  summarise(mean.z = mean(z.score)) %>% 
  ungroup()

head(modules_mean_z)


# Then we look at at which developmental stage and tissue is each module most highly expressed.
module_peak_exp <- modules_mean_z %>% 
  group_by(module) %>% 
  slice_max(order_by = mean.z, n = 1)

module_peak_exp

# You can also QC the clusters via a line graph. It will be too much to look at if graph all the modules, so let's just pick 2.
# This code chunk is very long, because a few things:
# I reordered x-axis to reflect the biological time sequence.
# Overlaid the average of clusters.
# Added a color strip at the bottom to annotate stages, which reduces the amount of text on the figure.
#"Control",	"cZ",	"cZ9G",	"DZ",	"DZ7G",	"DZ9G",	"iP","iP7G","iP9G",	"tZ","tZ7G","tZ9G"
######
#Works
#####

# Module quality control
my_network_modules %>% 
  filter( geneIDs == "AT4G29740" |
            geneIDs == "AT5G05860") ### CKX4 prsent in module 7, while UGT76C2 in Module 1. 

my_network_modules %>% 
  filter( geneIDs == "ATCG00280" |
            geneIDs == "ATCG00680") ### PSBC prsent in module 5, while PSBB in Module 3. 




my_network_modules %>% 
  filter( geneIDs == "AT3G54890" |
            geneIDs == "AT3G61470" | 
            geneIDs == "AT3G47470" |
          geneIDs == "AT5G54270" | 
            geneIDs == "AT3G08940" ) ### LHCs prsent in module 5 and 6


module_line_plot <- Exp_table_long_averaged_z_high_var_modules %>% 
  mutate(order_x = case_when(
    str_detect(Treatment, "Control") ~ 1,
    str_detect(Treatment, "tZ") ~ 2,
    str_detect(Treatment, "DZ") ~ 3,
    str_detect(Treatment, "iP") ~ 4,
    str_detect(Treatment, "cZ") ~ 5,
    str_detect(Treatment, "tZ7G") ~ 6,
    str_detect(Treatment, "tZ9G") ~ 7,
    str_detect(Treatment, "DZ7G") ~ 8,
    str_detect(Treatment, "DZ9G") ~ 9,
    str_detect(Treatment, "iP7G") ~ 10,
    str_detect(Treatment, "iP9G") ~ 11,
    str_detect(Treatment, "cZ9G") ~ 12,
    
  )) %>% 
  mutate(Treatment = reorder(Treatment, order_x)) %>% 
  filter(module == "5" |
           module == "7") %>% 
  ggplot(aes(x = Treatment, y = z.score)) +
  #do it by treatment
  facet_grid(module ~ Timepoint) +
  geom_line(aes(group = geneIDs), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "5" |
               module == "7" ) %>% 
      mutate(order_x = case_when(
        str_detect(Treatment, "Control") ~ 1,
        str_detect(Treatment, "tZ") ~ 2,
        str_detect(Treatment, "DZ") ~ 3,
        str_detect(Treatment, "iP") ~ 4,
        str_detect(Treatment, "cZ") ~ 5,
        str_detect(Treatment, "tZ7G") ~ 6,
        str_detect(Treatment, "tZ9G") ~ 7,
        str_detect(Treatment, "DZ7G") ~ 8,
        str_detect(Treatment, "DZ9G") ~ 9,
        str_detect(Treatment, "iP7G") ~ 10,
        str_detect(Treatment, "iP9G") ~ 11,
        str_detect(Treatment, "cZ9G") ~ 12,
      )) %>% 
      mutate(Treatment = reorder(Treatment, order_x)),
    aes(y = mean.z, group = module), 
    size = 1.1, alpha = 0.8
  ) +
  labs(x = NULL,
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_lines_color_strip <- expand.grid(
  Treatment = unique(Metadata$Treatment),
  Timepoint = unique(Metadata$Timepoint), 
  stringsAsFactors = F
) %>% 
  mutate(order_x = case_when(
    str_detect(Treatment, "Control") ~ 1,
    str_detect(Treatment, "tZ") ~ 2,
    str_detect(Treatment, "DZ") ~ 3,
    str_detect(Treatment, "iP") ~ 4,
    str_detect(Treatment, "cZ") ~ 5,
    str_detect(Treatment, "tZ7G") ~ 6,
    str_detect(Treatment, "tZ9G") ~ 7,
    str_detect(Treatment, "DZ7G") ~ 8,
    str_detect(Treatment, "DZ9G") ~ 9,
    str_detect(Treatment, "iP7G") ~ 10,
    str_detect(Treatment, "iP9G") ~ 11,
    str_detect(Treatment, "cZ9G") ~ 12,
  )) %>% 
  mutate(stage = factor(Treatment, levels = c(
    "Control",	"cZ",	"cZ9G",	"DZ",	"DZ7G",	"DZ9G",	"iP","iP7G","iP9G",	"tZ","tZ7G","tZ9G"
  ))) %>% 
  mutate(Treatment = reorder(Treatment, order_x)) %>% 
  ggplot(aes(x = Treatment, y = 1)) +
  facet_grid(. ~ Timepoint) +
  geom_tile(aes(fill = Treatment)) +
  scale_fill_manual(values = viridis(12, option = "D")) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    panel.spacing = unit(1, "lines")
  )

wrap_plots(module_line_plot, module_lines_color_strip,
           nrow = 2, heights = c(1, 0.08))
ggsave("Modules32.png", height = 5.5, width = 12,dpi=300, bg = "white")



## heat map - works?

modules_mean_z$mean.z %>% summary()
quantile(modules_mean_z$mean.z, 0.95)

# The 95th percentile of averaged z score is 1.3. We can probably roughly clipped the z-scores at 1.5 or -1.5
modules_mean_z <- modules_mean_z %>% 
  mutate(mean.z.clipped = case_when(
    mean.z > 0.9 ~ 0.9,
    mean.z < -0.9 ~ -0.9,
    T ~ mean.z
  ))

modules_mean_z <- modules_mean_z %>% 
  mutate(order_x = case_when(
    str_detect(Timepoint, "T2") ~ 1,
    str_detect(Timepoint, "T48") ~ 2,
    str_detect(Timepoint, "T96") ~ 3,
    str_detect(Timepoint, "T144") ~ 4,
  )) %>%  
  mutate(Timeopint = case_when(
    str_detect(Timepoint, "T2|T48|T96|T144") ~ str_sub(Treatment, start = 1, end = 2),
    T ~ Timepoint
  )) %>% 
  mutate(Timepoint = factor(Timepoint, levels = c(
    "T2",
    "T48",
    "T96",
    "T144"
  ))) %>% 
  mutate(Timepoint = reorder(Timepoint, order_x)) 

head(modules_mean_z)

#########################################################
#########################################################
####HeatMapCode for Treatment and timepoint (16) modules###############


module_heatmap <- modules_mean_z_reorded %>% 
  ggplot(aes(x = Treatment, y = as.factor(module))) +
  facet_grid(.~ Timepoint, scales = "free", space = "free") +
  geom_tile(aes(fill = mean.z.clipped), color = "grey80") +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), limits = c(-0.9, 0.9),
                       breaks = c(-0.9, 0, 0.9), labels = c("< -0.9", "0", "> 0.9")) +
  labs(x = NULL,
       y = "Module",
       fill = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_blank(),
    legend.position = "top",
    panel.spacing = unit(0.5, "lines") 
  )


heatmap_color_strip1 <- expand.grid(
  Treatment = unique(Metadata$Treatment),
  Timepoint = unique(Metadata$Timepoint), 
  stringsAsFactors = F
) %>% 
  mutate(order_x = case_when(
    str_detect(Treatment, "Control") ~ 1,
    str_detect(Treatment, "tZ") ~ 2,
    str_detect(Treatment, "iP") ~ 3,
    str_detect(Treatment, "DZ") ~ 4,
    str_detect(Treatment, "cZ") ~ 5,
    str_detect(Treatment, "tZ7G") ~ 6,
    str_detect(Treatment, "iP7G") ~ 7,
    str_detect(Treatment, "DZ7G") ~ 8,
    str_detect(Treatment, "tZ9G") ~ 9,
    str_detect(Treatment, "iP9G") ~ 10,
    str_detect(Treatment, "DZ9G") ~ 11,
    str_detect(Treatment, "cZ9G") ~ 12
  )) %>% 
  mutate(Treatment = factor(Treatment, levels = c(
    "Control",	"tZ","iP",	"DZ",	"cZ",	"tZ7G","iP7G","DZ7G",	"tZ9G","iP9G","DZ9G", "cZ9G"
    
  ))) %>% 
  mutate(Timepoint = reorder(Timepoint, order_x)) %>% 
  ggplot(aes(x = Treatment, y = 1)) +
  facet_grid(.~ Timepoint, scales = "free", space = "free") +
  geom_tile(aes(fill = Treatment)) +
  scale_fill_manual(values = brewer.pal(12, "Set3")) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    panel.spacing = unit(0.5, "lines"),
    legend.key.height = unit(0.75, "lines")
  )


heatmap_color_strip2 <- expand.grid(
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
  mutate(Timepoint = factor(Timepoint, levels = c(
    "T2",
    "T48",
    "T96",
    "T144"
  ), ordered = TRUE)) %>% 
  ggplot(aes(x = Treatment, y = 1)) +
  facet_grid(.~ Timepoint, scales = "free", space = "free") +
  geom_tile(aes(fill = Timepoint)) +
  scale_fill_manual(values = viridis(4, option = "D")) +
  labs(fill = "Timepoint") +
  guides(fill = guide_legend(nrow = 1)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    panel.spacing = unit(0.5, "lines"),
    legend.key.height = unit(0.75, "lines")
  )



wrap_plots(module_heatmap, heatmap_color_strip1, heatmap_color_strip2, 
           nrow = 3, heights = c(1, 0.08, 0.08), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )

ggsave("module_heatmap.png", height = 4.8, width = 12, bg = "white")
