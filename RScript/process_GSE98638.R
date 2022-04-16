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
main_directory <- 'scData/GSE98638'
file <- 'GSE98638_HCC.TCell.S5063.count.txt.gz'
labelfile <- 'GSE98638_labels.txt'
Split <- 0.75

# read count data
data98638 <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep="\t", header = TRUE)
rownames(data98638) <- make.names(data98638$symbol, unique = TRUE) 
data98638 <- data98638[, -c(1, 2)]
data98638 <- data98638[,sort(colnames(data98638))]


labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = "\t", header = TRUE)
rownames(labels) <- labels$UniqueCell_ID
labels<- labels[sort(rownames(labels)),]
labels <- as.data.frame(labels[, -c(1, 2, 3)])
colnames(labels) <- c("type")

count98638.list <- list(data98638)
label98638.list <- list(labels)

save_processed_data(count98638.list, label98638.list, folder_name='dataGSE98638')


Split <- 0.75
main_directory <- 'scData/GSE98638'
fprintf("%s\n",main_directory)
file <- 'GSE98638_HCC.TCell.S5063.count.txt.gz'
labelfile <- 'GSE98638_labels.txt'

# read count data
data <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep="\t", header = TRUE)

data <- data[data$symbol != "NA",] # Remove Gene/feature whose name is
rownames(data) <- make.names(data$symbol, unique = TRUE) 
data <-  subset(data, select = -c(geneID,symbol) ) #data[, -c(1, 2)]
data <- data[,sort(colnames(data))]
data <- as.data.frame(t(data))
fprintf("Shape of dataset is [%d x %d]\n", dim(data)[1],dim(data)[2] )

labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = "\t", header = TRUE)

labels <- labels[order(labels$UniqueCell_ID),]
rownames(labels) <- labels$UniqueCell_ID
# labels$majorCluster <- recode(labels$majorCluster, "CD4_C05-CXCR6" = "CD4+" , 
#                               "CD4_C11-IL10" = "CD4+" , "diverse.DN" = "diverse.DN" ,
#                               "diverse.other" = "diverse.other" , 
#                               "filtered" = "filtered" , "CD4_C06-CXCR5" = "CD4+" , 
#                               "CD4_C04-TCF7" = "CD4+" , "CD4_C02-ANXA1" = "CD4+" ,
#                               "CD8_C04-GZMK" = "CD8+" , "CD8_C05-CD6" = "CD8+", 
#                               "CD8_C06-CD160" = "CD8+"  , "CD8_C07-LAYN" = "CD8+" , 
#                               "CD8_C03-CX3CR1" = "CD8+" , "CD8_C01-LEF1" = "CD8+", 
#                               "CD8_C02-GPR183" = "CD8+", "diverse.DP" = "diverse.DP" , 
#                               "CD8_C08-SLC4A10" = "CD8+", "CD4_C12-CTLA4" = "CD4+" ,
#                               "CD4_C07-GZMK" = "CD4+", "CD4_C01-CCR7" = "CD4+", 
#                               "iNKT.CD4" = "iNKT.CD4", "CD4_C09-CXCL13" = "CD4+",
#                               "CD4_C10-FOXP3" = "CD4+", "MAIT.CD4" = "MAIT.CD4", 
#                               "CD4_C08-IL23R" = "CD4+", "CD4_C03-GNLY" = "CD4+", 
#                               "MAIT.other" = "MAIT.other", "MAIT.DN" = "MAIT.DN", 
#                               "iNKT.DN" = "iNKT.DN")


# labels$majorCluster <- recode(labels$majorCluster, "C01_CD8-LEF1" = "CD8+" ,
#                               "C02_CD8-CX3CR1" = "CD8+" , "C03_CD8-SLC4A10" = "CD8+",
#                               "C04_CD8-LAYN" = "CD8+", "C05_CD8-GZMK" = "CD8+",
#                               "C06_CD4-CCR7" = "CD4+", "C07_CD4-FOXP3" = "CD4+",
#                               "C08_CD4-CTLA4" = "CD4+", "C09_CD4-GZMA" = "CD4+",
#                               "C10_CD4-CXCL13" = "CD4+", "C11_CD4-GNLY" = "CD4+",
#                               "unknown" = "uknown"
#                               )

rows_labels <- rownames(labels)
rows_labels <- str_replace_all(rows_labels, '-', '.')
rows_data <- rownames(data)
rows_data <- str_replace_all(rows_data, '-', '.')

if (identical(rows_labels, rows_data) )
{
    fprintf("IDs are identical for labels and data\n")    
} else{
    fprintf("IDs are not identical for labels and data\n")    
}
rownames(data) <- rows_data
rownames(labels) <- rows_data

print(head(labels))
labels <-  subset(labels, select = -c(UniqueCell_ID,Patient, majorCluster) )
#labels <-  subset(labels, select = -c(UniqueCell_ID,Patient, sampleType) )
colnames(labels) <- c("type")

print(head(labels))
fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )

# Split the label data frame by category


refdf <- data.frame(matrix(ncol = length(colnames(data)), nrow = 0))
reflabel <- data.frame(matrix(ncol = 1, nrow = 0))
querydf <- data.frame(matrix(ncol = length(colnames(data)), nrow = 0))
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
    ref_data <- data[sampled_rows,]
    ref_labels <- as.data.frame(labels[sampled_rows,])
    rownames(ref_labels) <- sampled_rows
    
    refdf <- rbind(refdf, ref_data)
    reflabel<-rbind(reflabel, ref_labels)
    
    sampled_rows <- temporary_rownames[ -sample_index  ]
    query_data <- data[sampled_rows,]
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
identical(sort(c(QQ, RR)), rownames(data))

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