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
main_directory <- 'scData/GSE108989'
file <- 'GSE108989_CRC.TCell.S11138.count.txt.gz'
labelfile <- 'cell_labels.txt'

# read count data
data108989 <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep="\t", header = TRUE)
rownames(data108989) <- make.names(data108989$symbol, unique = TRUE) 
data108989 <- data108989[, -c(1, 2)]
data108989 <- data108989[,sort(colnames(data108989))]


labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = "\t", header = TRUE)
rownames(labels) <- labels$UniqueCell_ID
labels<- labels[sort(rownames(labels)),]
labels <- as.data.frame(labels[, -c(1, 2, 3)])
colnames(labels) <- c("type")


count108989.list <- list(data108989)
label108989.list <- list(labels)

save_processed_data(count108989.list, label108989.list, folder_name='data108989')

# Dividing data into reference and query
ref_cols <- sample(ncol(data108989), size=round(ncol(data108989)/2) )
ref_data <- data108989[, ref_cols, drop=FALSE]
ref_labels <- as.data.frame(labels[ref_cols,])
colnames(ref_labels) <- c("type")

query_data <- data108989[, -ref_cols, drop=FALSE]
query_labels <- as.data.frame(labels[-ref_cols, ])
colnames(query_labels) <- c("type")

refcount108989.list <- list(ref_data)
reflabel108989.list <- list(ref_labels)

querycount108989.list <- list(query_data)
querylabel108989.list <- list(query_labels)

save_processed_data(refcount108989.list, reflabel108989.list, folder_name='data108989_reference')
save_processed_data(querycount108989.list, querylabel108989.list, folder_name='data108989_query')
