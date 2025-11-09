rm(list = ls(all.names = TRUE)) # clears global environ.

# Load packages
library(tidyverse)
library(dplyr) 
library(biomaRt)
library(broom)
library(tibble)
library(edgeR)
library(factoextra)
library(biomaRt)
library(dendextend)

#########################
### Ingestion of data ### 
#########################

# # Create a list of the files from your target directory
# file_list <- list.files(path = "C:/---")
# 
# # Create a blank list for storage of .csv data
# list <- list()
# 
# # Loop through the length of the files and read that .csv and save it at location i within the list
# for (i in 1:length(file_list)){
# 
#   list[[i]] <- read.csv(paste0('C:/---/',
#                                file_list[i]))
# 
# }
# 
# # Merge all data frames together and remove duplicate samples
# RNAseq <- Reduce(function(x, y) merge(x, y, all = TRUE), list)

########################################################
### Select genes using edgeR function filterByExpr() ### 
########################################################

# Set seed
set.seed(1234)

# Create a synthetic data set with 100 genes, 5 samples per group, and 400/500 genes differentially expressed
RNAseq <- compcodeR::generateSyntheticData(dataset = "mydata", 
                                           n.vars = 500, 
                                           samples.per.cond = 20, 
                                           n.diffexp = 400)@count.matrix

# Shuffle around the synthetic data to create 4 random groups with 10 samples each
Control <- RNAseq[sample(nrow(RNAseq)), sample(ncol(RNAseq))][, 1:10] %>% 
  `rownames<-` (paste0('g', seq(1, 500, 1))) 

Treatment_1 <- data.frame(RNAseq[sample(nrow(RNAseq)), sample(ncol(RNAseq))][, 1:10]) %>% 
  `rownames<-` (paste0('g', seq(1, 500, 1))) 

Treatment_2 <- data.frame(RNAseq[sample(nrow(RNAseq)), sample(ncol(RNAseq))][, 1:10]) %>% 
  `rownames<-` (paste0('g', seq(1, 500, 1))) 

Treatment_3 <- data.frame(RNAseq[sample(nrow(RNAseq)), sample(ncol(RNAseq))][, 1:10]) %>% 
  `rownames<-` (paste0('g', seq(1, 500, 1)))

# Combine synthetic data groups into one data frame
RNAseq <- data.frame(Control,
                     Treatment_1,
                     Treatment_2,
                     Treatment_3) %>% t() 

# Lets get a list of real Ensembl IDs to randomly assign them to our synthetic data
Ensembl_ID <- read.csv('Ensembl_IDs.csv')

# Randomly sample genes to be included into the data set and set samples from 1 to 40
colnames(RNAseq) <- Ensembl_ID[sample(nrow(Ensembl_ID)), ][1:ncol(RNAseq)]
rownames(RNAseq) <- paste0('Sample_', seq(1, nrow(RNAseq), 1))

# Lets create 5 more dummy genes with lower counts to demonstrate filtering
RNAseq <- data.frame(Biological_Sample = factor(rep(seq(1, 10, 1), 4)), 
                     Treatment = factor(rep(c('Control', 
                                              'Treatment_1', 
                                              'Treatment_2', 
                                              'Treatment_3'), 
                                            each = 10)),
                     RNAseq, 
                     ENSG00000000061 = sample(1000, size = 40, replace = TRUE), 
                     ENSG00000000062 = sample(100, size = 40, replace = TRUE),
                     ENSG00000000063 = sample(10, size = 40, replace = TRUE),
                     ENSG00000000064 = sample(1, size = 40, replace = TRUE),
                     ENSG00000000065 = sample(0, size = 40, replace = TRUE))

head(RNAseq[1:5, 1:5])

# Create DGEList object
# DGElist accepts a matrix of nrow genes and ncol samples
# The group designation must be an equal length to ncol samples
d0 <- DGEList(as.matrix(t(RNAseq[3:ncol(RNAseq)])), 
              group = RNAseq$Biological_Sample)

# Calculate normalization factor
d0 <- calcNormFactors(d0)

# Use EdgeR function to automatically filter genes
# This gives logical values to keep or remove
keep.exprs <- filterByExpr(d0, 
                           group = RNAseq$Biological_Sample)

# Removed genes due to insufficient counts across samples
remove <- colnames(RNAseq[, !keep.exprs])

dim(RNAseq)

# Remove filtered genes from the dataset
RNAseq <- RNAseq[, keep.exprs]

dim(RNAseq)


#################################
### Check for outlier samples ### 
#################################

# Convert into a matrix
counts <- as.matrix(t(RNAseq[3:ncol(RNAseq)]))

# Create DGEList object
dataset <- DGEList(counts, 
                   group = RNAseq$Biological_Sample)

# Calculate normalization factor
dataset <- calcNormFactors(dataset)

# Calculate logCPM of filtered gene set
lcpm <- cpm(dataset, log = TRUE)

# Get count of number of samples in each group
counts <- RNAseq %>% group_by(Biological_Sample) %>% 
  summarise(total_count = n(),
            .groups = 'drop')

# Set group variable to a factor
Biological_Sample <- factor(RNAseq$Biological_Sample)

# Create a variable setting group to a number for MDS plot
col <- as.numeric(Biological_Sample)

# Create an MDS plot of log CPM data, labeled based on the group number set in col above
plotMDS(lcpm, 
        col = col,
        xlim = c(-4, 4))

# Run PCA
pca.res <- prcomp(t(lcpm), center = TRUE, scale = TRUE)

# Set theme
theme_set(theme_bw())

# Set colors object to the length of colors you need for below
colors <- c("forestgreen", 
            "steelblue4", 
            "darkorange4", 
            "red")

# Look up this table if you want different/additional symbols
shapes <- c(18, 9, 15, 12, 17)

# All groups cluster plot
fviz_pca_ind(pca.res,
             label = "none",
             habillage = RNAseq$Treatment, # What are ellipses drawn around?
             palette = colors, # Set palette to colors defined above
             addEllipses = TRUE, # Set to FALSE if no ellipse desired
             ellipse.type = 'convex' # Look up variable for options
) + 
  
  ggtitle("PCA - Clustering") + # Set title
  
  # Define axis information
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        panel.border = element_rect(fill = NA, color = "black", size = 0.3),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

########################################################
### Determine outliers through dendrogram clustering ### 
########################################################

# Add sample information back to log cpm data
lcpm <- data.frame(RNAseq[1:2], t(lcpm))

list <- as.list(unique(as.character(lcpm$Treatment)))

# Create dendrograms based on sample correlations from each group
dend <- lapply(list, function(x){
  
  tmp <- lcpm %>% subset(Treatment == x)
  names <- rownames(tmp)
  tmp_series <- tmp$series_id
  
  tmp <- tmp %>% 
    subset(select = -c(1:2)) %>% 
    sapply(as.numeric) %>% 
    as.data.frame(row.names = names)
  
  tmp_dend <- cor(t(tmp))
  tmp_dend <- as.dendrogram(hclust(as.dist(1 - tmp_dend)))
  
  useries = unique(tmp_series)
  
  # Match unique series to series list and record row location of each series
  series_match = useries[match(tmp_series, useries)]
  
  # Set colors
  colos <- colorspace::rainbow_hcl(length(tmp_series), 
                                   c = 160, 
                                   l  = 50)
  
  # Set colors to series ID
  names(colos) = tmp_series
  
  # Set matched colors
  series_color <- colos[series_match]
  
  # Create clusters
  clu = cutree(tmp_dend, 
               h = 0.25) # Height for cut (1 - correlation coefficient)
  
  # Set colors of labels
  labels_colors(tmp_dend) <- series_color[order.dendrogram(tmp_dend)]
  
  # Create dendrograms
  tmp_dend <- color_branches(tmp_dend, 
                             h = 0.25) # Height of cut (1 - correlation coefficient)
  
  par(mar = c(4,1,1,12))
  
  plot(tmp_dend,
       main = as.character(x),
       horiz = TRUE) +
    
    abline(v = 0.25, lty = 2)
  
  colored_bars(cbind(clu, 
                     series_color), 
               tmp_dend, 
               rowLabels = c("Cluster", "Series"), 
               horiz = TRUE)
  
  tmp <- lcpm %>% subset(Treatment == x)
  
  tmp <- cbind(clu, tmp)
  colnames(tmp)[1] <- 'cluster'
  tmp <- tmp %>% dplyr::select(Treatment, everything())
  
  largest_cluster = names(rev(sort(table(tmp$cluster))))[1]
  ww = which(tmp$cluster == largest_cluster)
  tmp[ww, ]
  
})

# Combine all selected samples into one place and remove dendrogram cluster information
Filtered_df <- do.call(rbind, dend) %>% subset(select = -c(2))

# Remove filtered samples from raw count data frame
RNAseq <- RNAseq[rownames(Filtered_df), ]

############################
### Export selected data ### 
############################

save(RNAseq, 
     file = 'Selected_Raw_Data.RDATA')

save(Filtered_df, 
     file = 'Selected_lcpm_Data.RDATA')
