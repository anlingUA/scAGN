library(Matrix)
library(stringr)
library(pracma)
library(RSpectra)
library(umap)
library(Rtsne)
library(bigstatsr)
library(data.table)
library(cluster)

source('save_processed_data.R')

# Folder count data and labels are
main_directory <- 'scData/GSE99254'
file <- 'GSE99254_NSCLC.TCell.S12346.count.txt.gz'
labelfile <- 'GSE99254_labels.txt'

# read count data
dataGSE99254 <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep="\t", header = TRUE)
rownames(dataGSE99254) <- make.names(dataGSE99254$symbol, unique = TRUE) 
dataGSE99254 <- dataGSE99254[, -c(1, 2)]
dataGSE99254 <- dataGSE99254[,sort(colnames(dataGSE99254))]


labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = "\t", header = TRUE)
rownames(labels) <- labels$UniqueCell_ID
labels<- labels[sort(rownames(labels)),]
labels <- as.data.frame(labels[, -c(1, 2, 3)])
colnames(labels) <- c("type")

countGSE99254.list <- list(dataGSE99254)
labelGSE99254.list <- list(labels)

save_processed_data(countGSE99254.list, labelGSE99254.list, folder_name='dataGSE99254')