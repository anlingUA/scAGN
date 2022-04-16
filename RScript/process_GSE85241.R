library(Matrix)
library(stringr)
library(pracma)
library(RSpectra)
library(umap)
library(Rtsne)
library(bigstatsr)
library(data.table)
library(cluster)
library(dplyr)

source('save_processed_data.R')

# Folder count data and labels are
main_directory <- 'scData/GSE85241'
file <- 'GSE85241_cellsystems_dataset_4donors_updated.csv.gz'
Split <- 0.75
# read count data
dataGSE85241 <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep="\t", header = TRUE)

splitdot <- function (s)
{
    n <- gregexpr("[.]", s)[[1]][1]
    if (n == -1)
    {
        return (s)
    }
    t <- substr(s, 1, n-1)
    return(t)
}

cellnames <- colnames(dataGSE85241) 

labels <- as.data.frame(splitdot(cellnames))
colnames(labels) <- c("type")
labels$type <- as.factor(labels$type)

tdataGSE85241 <- as.data.frame( t(dataGSE85241))

rows <- row.names(tdataGSE85241)
cols <- colnames(tdataGSE85241)
rownames(labels) <- rows

refdf <- data.frame(matrix(ncol = length(colnames(tdataGSE85241)), nrow = 0))
reflabel <- data.frame(matrix(ncol = 1, nrow = 0))
querydf <- data.frame(matrix(ncol = length(colnames(tdataGSE85241)), nrow = 0))
querylabel <- data.frame(matrix(ncol = 1, nrow = 0))
L<- split(labels, labels$type)
unique_types <- names(L)
for ( i in seq(1, length(L)))
{
    fprintf("Number of samples for type [%s] is %d\n", unique_types[i], length(L[[i]]$type))
    
    if ( length(L[[i]]$type) < 50)
    {
        print("Skipping .. Next")
        next
    }
    
    print("Moving Ahead ...")
    
    #data_list <- append(data_list, data[rownames(data) %in% rownames(L[[i]]),])
    temporary_rownames <- rownames(L[[i]])
    sample_index <- sample( 1:length(temporary_rownames), round(length(temporary_rownames)*Split))
    
    sampled_rows <- temporary_rownames[ sample_index  ]
    ref_data <- tdataGSE85241[sampled_rows,]
    ref_labels <- as.data.frame(labels[sampled_rows,])
    rownames(ref_labels) <- sampled_rows
    
    refdf <- rbind(refdf, ref_data)
    reflabel<-rbind(reflabel, ref_labels)
    
    sampled_rows <- temporary_rownames[ -sample_index  ]
    query_data <- tdataGSE85241[sampled_rows,]
    query_labels <- as.data.frame(labels[sampled_rows,])
    rownames(query_labels) <- sampled_rows
    
    querydf <- rbind(querydf, query_data)
    querylabel<-rbind(querylabel, query_labels)
}

colnames(reflabel) <- c("type")
colnames(querylabel) <- c("type")
K<- split(reflabel, reflabel$type)
Q <- split(querylabel, querylabel$type)

unique_types <- names(K)
for ( i in seq(1, length(K)))
{
    fprintf("Number of samples for type [%s] is %d\n", unique_types[i], length(K[[i]]$type))
}
unique_types <- names(Q)
for ( i in seq(1, length(Q)))
{
    fprintf("Number of samples for type [%s] is %d\n", unique_types[i], length(Q[[i]]$type))
}


RR = rownames(refdf)
QQ = rownames(querydf)
identical(sort(c(QQ, RR)), sort(rows))

## Shuffle
shuffled_row_index <- sample(nrow(refdf))
rows <- rownames(refdf)
shuffled_rows <- rows[shuffled_row_index]
refdf <- refdf[shuffled_rows,]
reflabel<- as.data.frame(reflabel[shuffled_rows,])
rownames(reflabel) <- shuffled_rows
colnames(reflabel) <- c("type")

shuffled_row_index <- sample(nrow(querydf))
rows <- rownames(querydf)
shuffled_rows <- rows[shuffled_row_index]
querydf <- querydf[shuffled_rows,]
querylabel<- as.data.frame(querylabel[shuffled_rows,])
rownames(querylabel) <- shuffled_rows
colnames(querylabel) <- c("type")

querylabel$type <- as.character(querylabel$type)
reflabel$type <- as.character(reflabel$type)

refdf2 <- refdf[ , colSums(is.na(refdf)) < nrow(refdf)]
querydf2 <- querydf[ , colSums(is.na(querydf)) < nrow(querydf)]

tref = transpose(refdf2)
rownames(tref) <- colnames(refdf2)
colnames(tref) <- rownames(refdf2)

tquery = transpose(querydf2)
rownames(tquery) <- colnames(querydf2)
colnames(tquery) <- rownames(querydf2)

mycount.list <- list(tref,tquery)
mylabel.list <- list(reflabel,querylabel)


countGSE85241.list <- list(dataGSE85241)
labelGSE85241.list <- list(labels)

save_processed_data(countGSE85241.list, labelGSE85241.list, folder_name='dataGSE85241')