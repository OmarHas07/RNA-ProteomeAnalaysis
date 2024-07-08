##WD 

setwd("/Users/omarhasannin/Documents/Research Projects/Chapter 1/RNA-SeqData/OPLS/CK-responseive Genes")


install.packages(c("caret", "tidyverse", "mixOmics"))
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("mixOmics")

## Librairies
library(caret)
library(tidyverse)
library(mixOmics)

### Loading Files
sample_info <- read_csv("Sample_info.csv")
wide_data <- read_csv("CK.df.csv")
dim(wide_data) # 5959 So this is just SAGs. 
wide_data <- as.data.frame(wide_data)


# Set the first column as row names
rownames(wide_data) <- wide_data[,1]

# Remove the first column from the data frame
wide_data <- wide_data[,-1]  # This drops the first column
# Check the structure of your data

# Ensure the sample order matches between the data and sample information
sample_info <- sample_info[match(colnames(wide_data), sample_info$Sample), ]

transposed_data <- t(wide_data)
transposed_data <- as.data.frame(transposed_data)

transposed_data$Sample <- rownames(transposed_data)

# Merge the data frames based on the "Sample" column
merged_df <- merge(transposed_data, sample_info, by = "Sample")

merged_df <- merged_df %>%  
  filter(Treatment != "Control")

merged_df$Treatment

T2 <- merged_df %>%
  filter(Timepoint == "T2") %>%
  column_to_rownames(var = "Treatment") %>% 
  dplyr::select(-c(Sample, Group, Timepoint))

T48 <- merged_df %>%  
  filter(Timepoint == "T48")  %>%
  column_to_rownames(var = "Treatment") %>% 
  dplyr::select(-c(Sample, Group, Timepoint))
T96 <- merged_df %>%  
  filter(Timepoint == "T96")  %>%
  column_to_rownames(var = "Treatment") %>% 
  dplyr::select(-c(Sample, Group, Timepoint))
T144 <- merged_df %>%  
  filter(Timepoint == "T144")  %>%
  column_to_rownames(var = "Treatment") %>% 
  dplyr::select(-c(Sample, Group, Timepoint))

### Making a combined DF for all times
# Rename columns to include time points
colnames(T2) <- paste0(colnames(T2), "_T2")
colnames(T48) <- paste0(colnames(T48), "_T48")
colnames(T96) <- paste0(colnames(T96), "_T96")
colnames(T144) <- paste0(colnames(T144), "_T144")

# Combine columns into one data frame
combined_df <- cbind(T2, T48, T96, T144)

# Convert row names to a column
combined_df <- combined_df %>% 
  rownames_to_column("Treatment")

# Add the "group" column
combined_df <- combined_df %>%
  mutate(group = case_when(
    #Treatment == "Control" ~ "Control",
    Treatment %in% c("cZ", "tZ", "DZ", "iP") ~ "Base forms",
    Treatment %in% c("tZ7G", "iP7G", "DZ7G") ~ "N7-Glucosides",
    Treatment %in% c("tZ9G", "iP9G", "DZ9G", "cZ9G") ~ "N9-Glucosides"
  ))


# Set the first column as row names
rownames(combined_df) <- combined_df[,1]

# Remove the first column from the data frame
combined_df <- combined_df[,-1]  # This drops the first column
# Check the structure of your data
## Making sure Group columns is there. 
combined_df$group


# Create a named list of data frames for each of the time points
df <- list("2 Hours" = T2, "48 Hours" = T48, "96 Hours" = T96, "144 Hours" = T144)

# set up a full design where every block is connected
design.df = matrix(1, ncol = length(df), nrow = length(df),
                   dimnames = list(names(df), names(df)))
diag(design.df) =  0
design.df

# set number of component per data set
ncomp = c(2)


####### Removing variables in each of the data frames 

# Identify zero variance features
nzv <- nearZeroVar(df$`2 Hours`)

# Print the zero variance features



# Get the indices of zero variance features
zero_var_indices <- nzv$Position

# Subset the data frame to include only non-zero variance features
df$`2 Hours` <- df$`2 Hours`[, -zero_var_indices]


### Checking zero var for T48 
nzv_t48 <- nearZeroVar(df$`48 Hours`)
zero_t48 <- nzv_t48$Position
df$`48 Hours` <- T48[, -zero_t48]

### Chekcing zero var for %96 
nzv_t96 <- nearZeroVar(df$`96 Hours`)
zero_t96 <- nzv_t96$Position
df$`96 Hours`<- T96[, -zero_t96]



## Checking var for T144 
nzv_t144 <- nearZeroVar(df$`144 Hours`)
zero_t144 <- nzv_t144$Position
df$`144 Hours` <- T144[, -zero_t144]

#####################



### Conducting the the analysis Group 

df.block.splsda = block.splsda(X = df, Y = combined_df$group,
                               ncomp = ncomp, design = design.df)




### Plotting the results 
png(filename = "sPLS-DA(wo.Control)---CK-Responsive.png", width = 12, height = 10, units = "in",  res = 300)
plotIndiv(df.block.splsda, ind.names = TRUE, ellipse = TRUE, legend = TRUE, style = "lattice")
dev.off()

png(filename = "Arrow-Plot-sPLS-DA(wo.Control).png", width = 12, height = 10, units = "in",  res = 300)
plotArrow(df.block.splsda, ind.names = TRUE, legend = TRUE)
dev.off()

df.block.splsda$loadings$T144
# illustrates coefficient weights in each block
plotLoadings(df.block.splsda, ncomp = 2, contrib = 'max')
plotVar(df.block.splsda, style = 'graphics', legend = TRUE)
network(df.block.splsda)
### Error zero varience 

### Cannot be done for some reason it's taking along time as well as the network plots 
circosPlot(df.block.splsda, cutoff = 0.7, line = TRUE, 
           color.blocks= c('darkorchid', 'brown1', 'lightgreen', 'blue'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)


df.block.splsda$variates
# Load necessary library
