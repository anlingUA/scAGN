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
main_directory <- 'scData/GSE118389'
file <- 'GSE118389_counts_rsem.txt.gz'

# read count data
dataGSE118389 <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep="\t", header = TRUE)

splitdot <- function (s)
{
    n <- gregexpr("[_]", s)[[1]][1]
    if (n == -1)
    {
        return (s)
    }
    t <- substr(s, 1, n-1)
    return(t)
}

cellnames <- colnames(dataGSE118389) 

labels <- as.data.frame(splitdot(cellnames))
colnames(labels) <- c("type")
labels$type <- as.factor(labels$type)


countGSE118389.list <- list(dataGSE118389)
labelGSE118389.list <- list(labels)

save_processed_data(countGSE118389.list, labelGSE118389.list, folder_name='dataGSE118389')