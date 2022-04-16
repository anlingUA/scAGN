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
main_directory <- 'scData/GSE72056'
file <- 'GSE72056_melanoma_single_cell_revised_v2.txt.gz'

# read count data
dataGSE72056 <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep="\t", header = TRUE)
rownames(dataGSE72056) <- make.names(dataGSE72056[,1], unique = TRUE)
rows <- rownames(dataGSE72056)
rows[2] <- "malignant_status"
rows[3] <- "non_malignant_celltype"
rownames(dataGSE72056) <- rows
dataGSE72056 <- dataGSE72056[, -c(1)]


tdataGSE72056 <- as.data.frame( t(dataGSE72056))
# remove any cells where malignant type is unresolved
tdataGSE72056 <- tdataGSE72056[tdataGSE72056$malignant_status != 0,]
tdataGSE72056$malignant_status <- as.factor(tdataGSE72056$malignant_status)
tdataGSE72056$celltype <- as.factor(tdataGSE72056$non_malignant_celltype)
#levels(tdataGSE72056$celltype) <- c("Malignant", "T", "B", "Macro", "Endo", "CAF", "NK")
labels <- tdataGSE72056$celltype
labels <- as.data.frame(labels)
colnames(labels) <- c("type")
tdataGSE72056 <- tdataGSE72056[, -c(1,2,3)]
dataGSE72056 <- t(tdataGSE72056)

countGSE72056.list <- list(dataGSE72056)
labelGSE72056.list <- list(labels)

label_mnemonics <-  c("Malignant", "T", "B", "Macro", "Endo", "CAF", "NK")
label_number <- c(0, 1, 2, 3, 4, 5, 6)
label_mnemonics.data <- data.frame(label_mnemonics, label_number)
#' create folder
dir.create('dataGSE72056'); 
write.csv(label_mnemonics.data,file=sprintf('%s/label_id.csv', 'dataGSE72056'),quote=F,row.names=T)

save_processed_data(countGSE72056.list, labelGSE72056.list, folder_name='dataGSE72056')