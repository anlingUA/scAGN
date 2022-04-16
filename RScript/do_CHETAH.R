suppressMessages(library(CHETAH))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(stringr))
suppressMessages(library(pracma))
suppressMessages(library(RSpectra))
suppressMessages(library(umap))
suppressMessages(library(Rtsne))
suppressMessages(library(bigstatsr))
suppressMessages(library(data.table))
suppressMessages(library(cluster))
suppressMessages(library("optparse"))
suppressMessages(library(dplyr))
suppressMessages(library(Dict));
options(future.globals.maxSize = 64000 * 4096^4)
source('data_preprocess_utility.R')
options(future.globals.maxSize = 64000 * 4096^4)
setwd(sprintf("%s/VersionControl/gene-cn/RScript",path.expand("~")))
set.seed(100)

CHETAH_prediction_inter <- function(refdf, querydf, reflabel,  querylabel, folder_to_save, dataset_name)
{
    refdf <- data.frame(matrix(ncol = length(colnames(data)), nrow = 0))
    reflabel <- data.frame(matrix(ncol = 1, nrow = 0))
    querydf <- data.frame(matrix(ncol = length(colnames(data)), nrow = 0))
    querylabel <- data.frame(matrix(ncol = 1, nrow = 0))
    L<- split(labels, labels$type)
    unique_types <- names(L)
    for ( i in seq(1, length(L)))
    {
        fprintf("Number of samples for type [%s] is %d\n", unique_types[i], length(L[[i]]$type))
        
        if ( length(L[[i]]$type) < type_threshold)
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
    print(identical(sort(c(QQ, RR)), rownames(data)))
    #stopifnot(identical(sort(c(QQ, RR)), rownames(data)))
    
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
    
    print("Reference Snippet:")
    print(tref[1:4, 1:4])
    print("Query Snippet:")
    print(tquery[1:4, 1:4])
    mycount.list <- list(tref,tquery)
    mylabel.list <- list(reflabel,querylabel)

    start_time <- Sys.time()	
    sce <- SingleCellExperiment(list(normcounts = mycount.list[[1]]), colData = data.frame(celltypes = mylabel.list[[1]]$type) )
    sce_test <- SingleCellExperiment(list(normcounts = mycount.list[[2]] ), colData = data.frame(celltypes = mylabel.list[[2]]$type  ))
    test <- CHETAHclassifier(input = sce_test, ref_cells = sce)
    end_time <- Sys.time()
    total_time_taken <- as.numeric(difftime(end_time,start_time,units = 'secs'))

    True_Labels_CHETAH <- mylabel.list[[2]]$type
    Pred_Labels_CHETAH <- test$celltype_CHETAH

    write.csv(mylabel.list[[2]], file = sprintf("%s/%s/true_labels_Chetah.csv", folder, dataset_name) )
    write.csv(as.data.frame(test$celltype_CHETAH), file = sprintf("%s/%s/predicted_labels_Chetah.csv", folder, dataset_name) )
    
    prediction_match = True_Labels_CHETAH == Pred_Labels_CHETAH
    final_accuracy <- (sum(prediction_match)/length(prediction_match))
    fprintf("Dataset = %s, Final Prediction is = %f\n\n", dataset_name, final_accuracy, file = sprintf("%s/Chetah_predict.txt",folder_to_save ),append=TRUE)
    fprintf("%s,%f\n", dataset_name, final_accuracy, file = sprintf("%s/Chetah_metric.txt",folder_to_save ),append=TRUE)
    fprintf("%s,%f\n", dataset_name, total_time_taken, file = sprintf("%s/Chetah_time_taken.txt",folder_to_save ),append=TRUE)

    print("Done!")
    
}


CHETAH_prediction_intra <- function(data, labels, folder_to_save, dataset_name)
{
    refdf <- data.frame(matrix(ncol = length(colnames(data)), nrow = 0))
    reflabel <- data.frame(matrix(ncol = 1, nrow = 0))
    querydf <- data.frame(matrix(ncol = length(colnames(data)), nrow = 0))
    querylabel <- data.frame(matrix(ncol = 1, nrow = 0))
    L<- split(labels, labels$type)
    unique_types <- names(L)
    for ( i in seq(1, length(L)))
    {
        fprintf("Number of samples for type [%s] is %d\n", unique_types[i], length(L[[i]]$type))
        
        if ( length(L[[i]]$type) < type_threshold)
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
    print(identical(sort(c(QQ, RR)), rownames(data)))
    #stopifnot(identical(sort(c(QQ, RR)), rownames(data)))
    
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
    
    print("Reference Snippet:")
    print(tref[1:4, 1:4])
    print("Query Snippet:")
    print(tquery[1:4, 1:4])
    mycount.list <- list(tref,tquery)
    mylabel.list <- list(reflabel,querylabel)
  	
    start_time <- Sys.time()
    sce <- SingleCellExperiment(list(normcounts = mycount.list[[1]]), colData = data.frame(celltypes = mylabel.list[[1]]$type) )
    sce_test <- SingleCellExperiment(list(normcounts = mycount.list[[2]] ), colData = data.frame(celltypes = mylabel.list[[2]]$type  ))
    test <- CHETAHclassifier(input = sce_test, ref_cells = sce)
    end_time <- Sys.time()
    total_time_taken <- as.numeric(difftime(end_time,start_time,units = 'secs'))

    True_Labels_CHETAH <- mylabel.list[[2]]$type
    Pred_Labels_CHETAH <- test$celltype_CHETAH
    write.csv(mylabel.list[[2]], file = sprintf("%s/%s/true_labels_Chetah.csv", folder, dataset_name) )
    write.csv(as.data.frame(test$celltype_CHETAH), file = sprintf("%s/%s/predicted_labels_Chetah.csv", folder, dataset_name) )

    prediction_match = True_Labels_CHETAH == Pred_Labels_CHETAH
    final_accuracy <- (sum(prediction_match)/length(prediction_match))
    fprintf("Dataset = %s, Final Prediction is = %f\n\n", dataset_name, final_accuracy, file = sprintf("%s/Chetah_predict.txt",folder_to_save ),append=TRUE)
    fprintf("%s,%f\n", dataset_name, final_accuracy, file = sprintf("%s/Chetah_metric.txt",folder_to_save ),append=TRUE)
    fprintf("%s,%f\n", dataset_name, total_time_taken, file = sprintf("%s/Chetah_time_taken.txt",folder_to_save ),append=TRUE)
    print("Done!")

}

option_list = list(
    make_option(c("-d", "--dataset"), type="integer", default=NULL,
                help="Dataset ID. Search path would <folder>/<dataset>/.\n\t\t\t ID Mapping: \n\t\t\t 1 = GSE108989\n\t\t\t 2 = GSE72056\n\t\t\t 3= GSE85241\n\t\t\t 4 = GSE98638\n\t\t\t 5 = GSE99254 \n\t\t\t 6 = GSE118389 \n\t\t\t 7 = GSE84133 \n\t\t\t 6 = GSE165092", metavar="integer"),
    make_option(c("-f", "--folder"), type="character", default=NULL,
                help="main folder where data is located. Search path would <folder>/<dataset>/", metavar="character"),
    make_option(c("-s", "--split"), type="double", default=NULL,
                help="Spliting factor to obtain reference data. Should be between [0, 1]. ", metavar="double"),
    make_option(c("-t", "--threshold"), type="double", default=NULL,
                help="Minimum number of cell-types to be present for each cell-type", metavar="double")

)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

id = opt$dataset
folder = opt$folder
Split = opt$split
type_threshold = opt$threshold

Dataset_list <- c('GSE108989', 'GSE72056', 'GSE85241', 'GSE98638', 'GSE99254',
                  'GSE118389', 'GSE84133', 'GSE165092', 'GSE83139', 'GSE85241_GSE81076',
                  'GSE81076_GSE85241', 'GSE124952', 'Xin', 'TM', 'Tasic',
                  'Segerstolpe', 'Romanov', 'PBMC_68k', 'Mouse_retina', 'Klein',
                  'villani_mgh_covid19', 'vento_pbmc_covid19', 'shalek_nasal_epithelia_covid19', 'blish_pbmc_covid19', 'meyer_airway_covid19',
                  'meyer_pbmc_covid19', 'vandamme_bronchoalveolar_lavage_fluid_covid19', 'aldinger_fetal_cerebellum', 'popescu19_fetal_liver',  'litvinukova20_heart',
                  'madissoon_oesophagus', 'vento_placenta_10X', 'martin_ileum', 'park_fetal_thymus', 'guo_testis',
                  'warner_salivery_gland', 'madissoon_spleen', 'miller_fetal_lung', 'vento_placenta_SmartSeq2','habib_brain',
                  'voigt_retina', 'james_colon_immune', 'stewart_kidney', 'vallier_restricted_gallbladder','menon_retina10xV3',
                  'byrd_gingiva', 'wang_rectum', 'lukowski_retina10xV2', 'Wyler_LungCellLine_Calu3','wang_colon',
                  'henry_prostate', 'lako_cornea', 'Wyler_LungCellLine_H1299', 'cheng18_skin','wang_ileum',
                  'smillie_colon_epithelial', 'macparland18_adult_liver',
                  'GSE83139_GSE84133', 'GSE84133_GSE83139', 'GSE103892')


if (Split > 1 || Split < 0)
{
    print("Splitting factor should be between 0 and 1.")
    quit(status=1)
}

dataset = Dataset_list[id]

fprintf("Folder = %s\n", folder)
fprintf("Dataset = %s\n", dataset)
fprintf("Reference-to-Query Split Ratio (applicable for intra-dataset) = %f\n", Split)

if(dataset == "GSE165092") {
    main_directory <- sprintf("%s/%s", folder, dataset)
    fprintf("%s\n",main_directory)
    file <- 'GSE165092_sc_count.tsv.gz'
    labelfile <- 'GSE165092_sc_meta.tsv.gz'
    data <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep="\t", header = TRUE)
    data <- as.data.frame(t(data))
    fprintf("Shape of dataset is [%d x %d]\n", dim(data)[1],dim(data)[2] )
    labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = "\t", header = TRUE)
    labels <-  subset(labels, select = -c(sample) )
    colnames(labels) <- c("type")
    rownames(data) <- rownames(labels)
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = "GSE165092")
} else  if(dataset == "GSE108989") {
    main_directory <- sprintf("%s/%s", folder, dataset)
    fprintf("%s\n",main_directory)
    file <- 'GSE108989_CRC.TCell.S11138.count.txt.gz'
    labelfile <- 'cell_labels.txt' 
    
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
    labels$majorCluster <- recode(labels$majorCluster, "CD4_C05-CXCR6" = "CD4+" , 
                                  "CD4_C11-IL10" = "CD4+" , "diverse.DN" = "diverse.DN" ,
                                  "diverse.other" = "diverse.other" , 
                                  "filtered" = "filtered" , "CD4_C06-CXCR5" = "CD4+" , 
                                  "CD4_C04-TCF7" = "CD4+" , "CD4_C02-ANXA1" = "CD4+" ,
                                  "CD8_C04-GZMK" = "CD8+" , "CD8_C05-CD6" = "CD8+", 
                                  "CD8_C06-CD160" = "CD8+"  , "CD8_C07-LAYN" = "CD8+" , 
                                  "CD8_C03-CX3CR1" = "CD8+" , "CD8_C01-LEF1" = "CD8+", 
                                  "CD8_C02-GPR183" = "CD8+", "diverse.DP" = "diverse.DP" , 
                                  "CD8_C08-SLC4A10" = "CD8+", "CD4_C12-CTLA4" = "CD4+" ,
                                  "CD4_C07-GZMK" = "CD4+", "CD4_C01-CCR7" = "CD4+", 
                                  "iNKT.CD4" = "iNKT.CD4", "CD4_C09-CXCL13" = "CD4+",
                                  "CD4_C10-FOXP3" = "CD4+", "MAIT.CD4" = "MAIT.CD4", 
                                  "CD4_C08-IL23R" = "CD4+", "CD4_C03-GNLY" = "CD4+", 
                                  "MAIT.other" = "MAIT.other", "MAIT.DN" = "MAIT.DN", 
                                  "iNKT.DN" = "iNKT.DN")
    
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
    # labels <-  subset(labels, select = -c(UniqueCell_ID,Patient_ID, majorCluster) )
    labels <-  subset(labels, select = -c(UniqueCell_ID,Patient_ID, sampleType) )
    colnames(labels) <- c("type")
    print(head(labels))
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = "GSE108989")
} else if(dataset == "GSE72056"){
    # Folder count data and labels are
    main_directory <- sprintf("%s/%s", folder, dataset)
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
    tdataGSE72056$celltype <- recode(tdataGSE72056$celltype, "0"="Malignant", 
                                     "1" = "T",
                                     "2" = "B", 
                                     "3" = "Macro", 
                                     "4" = "Endo",
                                     "5" = "CAF",
                                     "6" = "NK")
    
    #levels(tdataGSE72056$celltype) <- c("Malignant", "T", "B", "Macro", "Endo", "CAF", "NK")
    
    rows <- row.names(tdataGSE72056)
    cols <- colnames(tdataGSE72056)
    labels <- tdataGSE72056$celltype
    labels <- as.data.frame(labels)
    labels <- tdataGSE72056$celltype
    labels <- as.data.frame(labels)
    colnames(labels) <- c("type")
    rownames(labels) <- rows
    tdataGSE72056 <- subset(tdataGSE72056, select = -c(tumor, malignant_status,non_malignant_celltype, celltype) )
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    CHETAH_prediction_intra(tdataGSE72056, labels, folder_to_save = "./", dataset_name = "GSE72056")
} else if(dataset == "GSE85241"){
    # Folder count data and labels are
    main_directory <- sprintf("%s/%s", folder, dataset)
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
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    CHETAH_prediction_intra(tdataGSE85241, labels, folder_to_save = "./", dataset_name = "GSE85241")
} else if(dataset == "GSE98638") {
    
    main_directory <- sprintf("%s/%s", folder, dataset)
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
    #labels <-  subset(labels, select = -c(UniqueCell_ID,Patient, majorCluster) )
    labels <-  subset(labels, select = -c(UniqueCell_ID,Patient, sampleType) )
    colnames(labels) <- c("type")
    print(head(labels))
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = "GSE98638")
    
} else if(dataset == "GSE99254") {
    main_directory <- sprintf("%s/%s", folder, dataset)
    fprintf("%s\n",main_directory)
    file <- 'GSE99254_NSCLC.TCell.S12346.count.txt.gz'
    labelfile <- 'GSE99254_labels.txt'
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
    #labels <-  subset(labels, select = -c(UniqueCell_ID,Patient, majorCluster) )
    labels <-  subset(labels, select = -c(UniqueCell_ID,Patient, sampleType) )
    colnames(labels) <- c("type")
    print(head(labels))
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = "GSE99254")
    
} else if(dataset == "GSE118389"){
    
    main_directory <- sprintf("%s/%s", folder, dataset)
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
    tdataGSE118389 <- as.data.frame( t(dataGSE118389))
    rows <- row.names(tdataGSE118389)
    cols <- colnames(tdataGSE118389)
    rownames(labels) <- rows
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    CHETAH_prediction_intra(tdataGSE118389, labels, folder_to_save = "./", dataset_name = "GSE118389")
} else if(dataset == "GSE83139"){
    main_directory <- sprintf("%s/%s", folder, dataset)
    file <- 'GSE83139_tbx-v-f-norm-ntv-cpms.csv.gz'
    # read count data
    dataGSE83139 <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep="\t", header = TRUE)
    dataGSE83139 <- subset(dataGSE83139, select = -c(chr, start, end, strand, transcript, part) )
    dataGSE83139 <- dataGSE83139[dataGSE83139$gene != "NA",] # Remove Gene/feature whose name is
    dataGSE83139 <- na.omit(dataGSE83139)
    rownames(dataGSE83139) <- dataGSE83139$gene
    dataGSE83139 <- subset(dataGSE83139, select = -c(gene))
    dataGSE83139 <- dataGSE83139[,sort(colnames(dataGSE83139))]
    dataGSE83139 <- as.data.frame(t(dataGSE83139))
    
    label_file <- "cell_annotations.csv"
    labels <- read.table(sprintf("%s/%s", main_directory, label_file),sep=",", header = TRUE)
    
    labels$Barcode <- unlist(str_split(labels$Barcode, '_'))[c(FALSE, TRUE)]
    labels <- labels[order(labels$Barcode),]
    rownames(labels) <- labels$Barcode
    labels <- subset(labels, select = -c(Barcode, Set, Platform))
    de <- merge(dataGSE83139, labels, by=0, all=TRUE)
    de <- na.omit(de)
    row.names(de) <- de$Row.names
    labels <- de$Type
    dataGSE83139 <- subset(de, select = -c(Row.names, Type))
    labels <- as.data.frame(labels)
    row.names(labels) <- row.names(dataGSE83139)
    colnames(labels) <- c("type")
    fprintf("Shape of dataset is [%d x %d]\n", dim(dataGSE83139)[1],dim(dataGSE83139)[2] )
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    CHETAH_prediction_intra(dataGSE83139, labels, folder_to_save = "./", dataset_name = "GSE83139")
} else if(dataset == "GSE85241_GSE81076"){
    
    folder <- 'scData'
    dataset_ref = "GSE85241"
    main_directory_ref <- sprintf("%s/%s", folder, dataset_ref)
    file <- 'GSE85241_cellsystems_dataset_4donors_updated.csv.gz'
    # read count data
    dataGSE85241 <- read.table(gzfile(sprintf("%s/%s",  main_directory_ref, file)),sep="\t", header = TRUE)
    dataGSE85241 <- na.omit(dataGSE85241)
    dataGSE85241 <- dataGSE85241[,sort(colnames(dataGSE85241))]
    dataGSE85241 <- as.data.frame(t(dataGSE85241))
    label_file <- "cell_annotation.csv"
    labelsGSE85241 <- read.table(sprintf("%s/%s",  main_directory_ref, label_file),sep=",", header = TRUE)
    labelsGSE85241$Barcode <- unlist(str_split(labelsGSE85241$Barcode, 'GSE85241_'))[c(FALSE, TRUE)]
    labelsGSE85241 <- labelsGSE85241[order(labelsGSE85241$Barcode),]
    rownames(labelsGSE85241) <- labelsGSE85241$Barcode
    labelsGSE85241 <- subset(labelsGSE85241, select = -c(Barcode, Set, Platform))
    deGSE85241 <- merge(dataGSE85241, labelsGSE85241, by=0, all=TRUE)
    deGSE85241 <- na.omit(deGSE85241)
    row.names(deGSE85241) <- deGSE85241$Row.names
    labelsGSE85241 <- deGSE85241$Type
    dataGSE85241 <- subset(deGSE85241, select = -c(Row.names, Type))
    labelsGSE85241 <- as.data.frame(labelsGSE85241)
    row.names(labelsGSE85241) <- row.names(deGSE85241)
    colnames(labelsGSE85241) <- c("type")
    fprintf("Shape of dataset is [%d x %d]\n", dim(dataGSE85241)[1],dim(dataGSE85241)[2] )
    fprintf("Shape of labels df is [%d x %d]\n", dim(labelsGSE85241)[1],dim(labelsGSE85241)[2] )
    
    dataset_query = "GSE81076"
    main_directory_query <- sprintf("%s/%s", folder, dataset_query)
    file <- 'GSE81076_D2_3_7_10_17.txt.gz'
    # read count data
    dataGSE81076 <- read.table(gzfile(sprintf("%s/%s",  main_directory_query, file)),sep="\t", header = TRUE)
    dataGSE81076 <- dataGSE81076[dataGSE81076$X != "NA",] # Remove Gene/feature whose name is 
    dataGSE81076 <- na.omit(dataGSE81076)
    rownames(dataGSE81076) <- dataGSE81076$X
    dataGSE81076 <- subset(dataGSE81076, select = -c(X))
    dataGSE81076 <- dataGSE81076[,sort(colnames(dataGSE81076))]
    dataGSE81076 <- as.data.frame(t(dataGSE81076))
    label_file <- "cell_annotation.csv"
    labelsGSE81076 <- read.table(sprintf("%s/%s",  main_directory_query, label_file),sep=",", header = TRUE)
    labelsGSE81076$Barcode <- unlist(str_split(labelsGSE81076$Barcode, 'GSE81076_'))[c(FALSE, TRUE)]
    labelsGSE81076 <- labelsGSE81076[order(labelsGSE81076$Barcode),]
    rownames(labelsGSE81076) <- labelsGSE81076$Barcode
    labelsGSE81076 <- subset(labelsGSE81076, select = -c(Barcode, Set, Platform))
    deGSE81076 <- merge(dataGSE81076, labelsGSE81076, by=0, all=TRUE)
    deGSE81076 <- na.omit(deGSE81076)
    row.names(deGSE81076) <- deGSE81076$Row.names
    labelsGSE81076 <- deGSE81076$Type
    dataGSE81076 <- subset(deGSE81076, select = -c(Row.names, Type))
    labelsGSE81076 <- as.data.frame(labelsGSE81076)
    row.names(labelsGSE81076) <- row.names(deGSE81076)
    colnames(labelsGSE81076) <- c("type")
    fprintf("Shape of dataset is [%d x %d]\n", dim(dataGSE81076)[1],dim(dataGSE81076)[2] )
    fprintf("Shape of labels df is [%d x %d]\n", dim(labelsGSE81076)[1],dim(labelsGSE81076)[2] )
    
    # now keep common genes only
    common_genes <- intersect(colnames(dataGSE85241), colnames(dataGSE81076))
    dataGSE85241 <- dataGSE85241[,common_genes]
    dataGSE81076 <- dataGSE81076[,common_genes]
    CHETAH_prediction_inter(dataGSE85241, dataGSE81076, labelsGSE85241,  labelsGSE81076, folder_to_save = "./", dataset_name = "GSE85241_GSE81076")
} else if(dataset == "GSE81076_GSE85241"){
    
    folder <- 'scData'
    dataset_query = "GSE85241"
    main_directory_query <- sprintf("%s/%s", folder, dataset_query)
    file <- 'GSE85241_cellsystems_dataset_4donors_updated.csv.gz'
    # read count data
    dataGSE85241 <- read.table(gzfile(sprintf("%s/%s",  main_directory_query, file)),sep="\t", header = TRUE)
    dataGSE85241 <- na.omit(dataGSE85241)
    dataGSE85241 <- dataGSE85241[,sort(colnames(dataGSE85241))]
    dataGSE85241 <- as.data.frame(t(dataGSE85241))
    label_file <- "cell_annotation.csv"
    labelsGSE85241 <- read.table(sprintf("%s/%s",  main_directory_query, label_file),sep=",", header = TRUE)
    labelsGSE85241$Barcode <- unlist(str_split(labelsGSE85241$Barcode, 'GSE85241_'))[c(FALSE, TRUE)]
    labelsGSE85241 <- labelsGSE85241[order(labelsGSE85241$Barcode),]
    rownames(labelsGSE85241) <- labelsGSE85241$Barcode
    labelsGSE85241 <- subset(labelsGSE85241, select = -c(Barcode, Set, Platform))
    deGSE85241 <- merge(dataGSE85241, labelsGSE85241, by=0, all=TRUE)
    deGSE85241 <- na.omit(deGSE85241)
    row.names(deGSE85241) <- deGSE85241$Row.names
    labelsGSE85241 <- deGSE85241$Type
    dataGSE85241 <- subset(deGSE85241, select = -c(Row.names, Type))
    labelsGSE85241 <- as.data.frame(labelsGSE85241)
    row.names(labelsGSE85241) <- row.names(deGSE85241)
    colnames(labelsGSE85241) <- c("type")
    fprintf("Shape of dataset is [%d x %d]\n", dim(dataGSE85241)[1],dim(dataGSE85241)[2] )
    fprintf("Shape of labels df is [%d x %d]\n", dim(labelsGSE85241)[1],dim(labelsGSE85241)[2] )
    
    dataset_ref = "GSE81076"
    main_directory_ref <- sprintf("%s/%s", folder, dataset_ref)
    file <- 'GSE81076_D2_3_7_10_17.txt.gz'
    # read count data
    dataGSE81076 <- read.table(gzfile(sprintf("%s/%s",  main_directory_ref, file)),sep="\t", header = TRUE)
    dataGSE81076 <- dataGSE81076[dataGSE81076$X != "NA",] # Remove Gene/feature whose name is 
    dataGSE81076 <- na.omit(dataGSE81076)
    rownames(dataGSE81076) <- dataGSE81076$X
    dataGSE81076 <- subset(dataGSE81076, select = -c(X))
    dataGSE81076 <- dataGSE81076[,sort(colnames(dataGSE81076))]
    dataGSE81076 <- as.data.frame(t(dataGSE81076))
    label_file <- "cell_annotation.csv"
    labelsGSE81076 <- read.table(sprintf("%s/%s",  main_directory_ref, label_file),sep=",", header = TRUE)
    labelsGSE81076$Barcode <- unlist(str_split(labelsGSE81076$Barcode, 'GSE81076_'))[c(FALSE, TRUE)]
    labelsGSE81076 <- labelsGSE81076[order(labelsGSE81076$Barcode),]
    rownames(labelsGSE81076) <- labelsGSE81076$Barcode
    labelsGSE81076 <- subset(labelsGSE81076, select = -c(Barcode, Set, Platform))
    deGSE81076 <- merge(dataGSE81076, labelsGSE81076, by=0, all=TRUE)
    deGSE81076 <- na.omit(deGSE81076)
    row.names(deGSE81076) <- deGSE81076$Row.names
    labelsGSE81076 <- deGSE81076$Type
    dataGSE81076 <- subset(deGSE81076, select = -c(Row.names, Type))
    labelsGSE81076 <- as.data.frame(labelsGSE81076)
    row.names(labelsGSE81076) <- row.names(deGSE81076)
    colnames(labelsGSE81076) <- c("type")
    fprintf("Shape of dataset is [%d x %d]\n", dim(dataGSE81076)[1],dim(dataGSE81076)[2] )
    fprintf("Shape of labels df is [%d x %d]\n", dim(labelsGSE81076)[1],dim(labelsGSE81076)[2] )
    
    # now keep common genes only
    common_genes <- intersect(colnames(dataGSE85241), colnames(dataGSE81076))
    dataGSE85241 <- dataGSE85241[,common_genes]
    dataGSE81076 <- dataGSE81076[,common_genes]
    CHETAH_prediction_inter(dataGSE81076, dataGSE85241, labelsGSE81076, labelsGSE85241, folder_to_save = "./", dataset_name = "GSE81076_GSE85241")
} else if(dataset == "GSE83139_GSE84133") {
    
    dataset_ref = "GSE83139"
    main_directory_ref <- sprintf("%s/%s", folder, dataset_ref)
    file <- 'GSE83139_tbx-v-f-norm-ntv-cpms.csv.gz'
    # read count data
    dataGSE83139 <- read.table(gzfile(sprintf("%s/%s", main_directory_ref, file)),sep="\t", header = TRUE)
    dataGSE83139 <- subset(dataGSE83139, select = -c(chr, start, end, strand, transcript, part) )
    dataGSE83139 <- dataGSE83139[dataGSE83139$gene != "NA",] # Remove Gene/feature whose name is
    dataGSE83139 <- na.omit(dataGSE83139)
    rownames(dataGSE83139) <- dataGSE83139$gene
    dataGSE83139 <- subset(dataGSE83139, select = -c(gene))
    dataGSE83139 <- dataGSE83139[,sort(colnames(dataGSE83139))]
    dataGSE83139 <- as.data.frame(t(dataGSE83139))
    label_file <- "cell_annotations.csv"
    labelsGSE83139 <- read.table(sprintf("%s/%s", main_directory_ref, label_file),sep=",", header = TRUE)
    labelsGSE83139$Barcode <- unlist(str_split(labelsGSE83139$Barcode, '_'))[c(FALSE, TRUE)]
    labelsGSE83139 <- labelsGSE83139[order(labelsGSE83139$Barcode),]
    rownames(labelsGSE83139) <- labelsGSE83139$Barcode
    labelsGSE83139 <- subset(labelsGSE83139, select = -c(Barcode, Set, Platform))
    de <- merge(dataGSE83139, labelsGSE83139, by=0, all=TRUE)
    de <- na.omit(de)
    row.names(de) <- de$Row.names
    labelsGSE83139 <- de$Type
    dataGSE83139 <- subset(de, select = -c(Row.names, Type))
    labelsGSE83139 <- as.data.frame(labelsGSE83139)
    row.names(labelsGSE83139) <- row.names(dataGSE83139)
    colnames(labelsGSE83139) <- c("type")
    fprintf("Shape of dataset is [%d x %d]\n", dim(dataGSE83139)[1],dim(dataGSE83139)[2] )
    fprintf("Shape of labels df is [%d x %d]\n", dim(labelsGSE83139)[1],dim(labelsGSE83139)[2] )
    
    
    dataset_query = "GSE84133"
    main_directory_query <- sprintf("%s/%s", folder, dataset_query)
    file1 <- 'GSM2230757_human1_umifm_counts.csv.gz'
    data1 <- read.table(gzfile(sprintf("%s/%s",  main_directory_query, file1)),sep=",", header = TRUE)
    file2 <- 'GSM2230758_human2_umifm_counts.csv.gz'
    data2 <- read.table(gzfile(sprintf("%s/%s",  main_directory_query, file2)),sep=",", header = TRUE)
    file3 <- 'GSM2230759_human3_umifm_counts.csv.gz'
    data3 <- read.table(gzfile(sprintf("%s/%s",  main_directory_query, file3)),sep=",", header = TRUE)
    file4 <- 'GSM2230760_human4_umifm_counts.csv.gz'
    data4 <- read.table(gzfile(sprintf("%s/%s",  main_directory_query, file4)),sep=",", header = TRUE)
    dataGSE84133 <- rbind(data1, data2, data3, data4)
    rownames(dataGSE84133) <- dataGSE84133$X
    labelsGSE84133 <- dataGSE84133$assigned_cluster
    dataGSE84133 <- subset(dataGSE84133, select = -c(X, barcode, assigned_cluster))
    labelsGSE84133 <- as.data.frame(labelsGSE84133)
    row.names(labelsGSE84133) <- row.names(dataGSE84133)
    colnames(labelsGSE84133) <- c("type")
    fprintf("Shape of dataset is [%d x %d]\n", dim(dataGSE84133)[1],dim(dataGSE84133)[2] )
    fprintf("Shape of labels df is [%d x %d]\n", dim(labelsGSE84133)[1],dim(labelsGSE84133)[2] )
    
    # now keep common genes only
    common_genes <- intersect(colnames(dataGSE83139), colnames(dataGSE84133))
    fprintf("Common genes length are [%d]\n",length(common_genes) )
    dataGSE83139 <- dataGSE83139[,common_genes]
    dataGSE84133 <- dataGSE84133[,common_genes]
    CHETAH_prediction_inter(dataGSE83139, dataGSE84133, labelsGSE83139, labelsGSE84133, folder_to_save = "./", dataset_name = "GSE83139_GSE84133")
} else if(dataset == "GSE84133_GSE83139") {
    
    dataset_ref = "GSE84133"
    main_directory_ref <- sprintf("%s/%s", folder, dataset_ref)
    file1 <- 'GSM2230757_human1_umifm_counts.csv.gz'
    data1 <- read.table(gzfile(sprintf("%s/%s",  main_directory_ref, file1)),sep=",", header = TRUE)
    file2 <- 'GSM2230758_human2_umifm_counts.csv.gz'
    data2 <- read.table(gzfile(sprintf("%s/%s",  main_directory_ref, file2)),sep=",", header = TRUE)
    file3 <- 'GSM2230759_human3_umifm_counts.csv.gz'
    data3 <- read.table(gzfile(sprintf("%s/%s",  main_directory_ref, file3)),sep=",", header = TRUE)
    file4 <- 'GSM2230760_human4_umifm_counts.csv.gz'
    data4 <- read.table(gzfile(sprintf("%s/%s",  main_directory_ref, file4)),sep=",", header = TRUE)
    dataGSE84133 <- rbind(data1, data2, data3, data4)
    rownames(dataGSE84133) <- dataGSE84133$X
    labelsGSE84133 <- dataGSE84133$assigned_cluster
    dataGSE84133 <- subset(dataGSE84133, select = -c(X, barcode, assigned_cluster))
    labelsGSE84133 <- as.data.frame(labelsGSE84133)
    row.names(labelsGSE84133) <- row.names(dataGSE84133)
    colnames(labelsGSE84133) <- c("type")
    fprintf("Shape of dataset is [%d x %d]\n", dim(dataGSE84133)[1],dim(dataGSE84133)[2] )
    fprintf("Shape of labels df is [%d x %d]\n", dim(labelsGSE84133)[1],dim(labelsGSE84133)[2] )
    
    dataset_query = "GSE83139"
    main_directory_query <- sprintf("%s/%s", folder, dataset_query)
    file <- 'GSE83139_tbx-v-f-norm-ntv-cpms.csv.gz'
    # read count data
    dataGSE83139 <- read.table(gzfile(sprintf("%s/%s", main_directory_query, file)),sep="\t", header = TRUE)
    dataGSE83139 <- subset(dataGSE83139, select = -c(chr, start, end, strand, transcript, part) )
    dataGSE83139 <- dataGSE83139[dataGSE83139$gene != "NA",] # Remove Gene/feature whose name is
    dataGSE83139 <- na.omit(dataGSE83139)
    rownames(dataGSE83139) <- dataGSE83139$gene
    dataGSE83139 <- subset(dataGSE83139, select = -c(gene))
    dataGSE83139 <- dataGSE83139[,sort(colnames(dataGSE83139))]
    dataGSE83139 <- as.data.frame(t(dataGSE83139))
    label_file <- "cell_annotations.csv"
    labelsGSE83139 <- read.table(sprintf("%s/%s", main_directory_query, label_file),sep=",", header = TRUE)
    labelsGSE83139$Barcode <- unlist(str_split(labelsGSE83139$Barcode, '_'))[c(FALSE, TRUE)]
    labelsGSE83139 <- labelsGSE83139[order(labelsGSE83139$Barcode),]
    rownames(labelsGSE83139) <- labelsGSE83139$Barcode
    labelsGSE83139 <- subset(labelsGSE83139, select = -c(Barcode, Set, Platform))
    de <- merge(dataGSE83139, labelsGSE83139, by=0, all=TRUE)
    de <- na.omit(de)
    row.names(de) <- de$Row.names
    labelsGSE83139 <- de$Type
    dataGSE83139 <- subset(de, select = -c(Row.names, Type))
    labelsGSE83139 <- as.data.frame(labelsGSE83139)
    row.names(labelsGSE83139) <- row.names(dataGSE83139)
    colnames(labelsGSE83139) <- c("type")
    fprintf("Shape of dataset is [%d x %d]\n", dim(dataGSE83139)[1],dim(dataGSE83139)[2] )
    fprintf("Shape of labels df is [%d x %d]\n", dim(labelsGSE83139)[1],dim(labelsGSE83139)[2] )
    
    # now keep common genes only
    common_genes <- intersect(colnames(dataGSE83139), colnames(dataGSE84133))
    fprintf("Common genes length are [%d]\n",length(common_genes) )
    dataGSE83139 <- dataGSE83139[,common_genes]
    dataGSE84133 <- dataGSE84133[,common_genes]
    CHETAH_prediction_inter(dataGSE84133, dataGSE83139, labelsGSE84133, labelsGSE83139, folder_to_save = "./", dataset_name = "GSE84133_GSE83139")
} else if(dataset == "GSE84133_human_mouse") {
    
    dataset_ref = "GSE84133"
    main_directory_ref <- sprintf("%s/%s", folder, dataset_ref)
    file1 <- 'GSM2230757_human1_umifm_counts.csv.gz'
    data1 <- read.table(gzfile(sprintf("%s/%s",  main_directory_ref, file1)),sep=",", header = TRUE)
    file2 <- 'GSM2230758_human2_umifm_counts.csv.gz'
    data2 <- read.table(gzfile(sprintf("%s/%s",  main_directory_ref, file2)),sep=",", header = TRUE)
    file3 <- 'GSM2230759_human3_umifm_counts.csv.gz'
    data3 <- read.table(gzfile(sprintf("%s/%s",  main_directory_ref, file3)),sep=",", header = TRUE)
    file4 <- 'GSM2230760_human4_umifm_counts.csv.gz'
    data4 <- read.table(gzfile(sprintf("%s/%s",  main_directory_ref, file4)),sep=",", header = TRUE)
    dataGSE84133_human <- rbind(data1, data2, data3, data4)
    rownames(dataGSE84133_human) <- dataGSE84133_human$X
    labelsGSE84133_human <- dataGSE84133_human$assigned_cluster
    dataGSE84133_human <- subset(dataGSE84133_human, select = -c(X, barcode, assigned_cluster))
    labelsGSE84133_human <- as.data.frame(labelsGSE84133_human)
    row.names(labelsGSE84133_human) <- row.names(dataGSE84133_human)
    colnames(labelsGSE84133_human) <- c("type")
    fprintf("Shape of dataset is [%d x %d]\n", dim(dataGSE84133_human)[1],dim(dataGSE84133_human)[2] )
    fprintf("Shape of labels df is [%d x %d]\n", dim(labelsGSE84133_human)[1],dim(labelsGSE84133_human)[2] )
    
    file5 <- 'GSM2230761_mouse1_umifm_counts.csv.gz'
    data5 <- read.table(gzfile(sprintf("%s/%s",  main_directory_ref, file5)),sep=",", header = TRUE)
    file6 <- 'GSM2230762_mouse2_umifm_counts.csv.gz'
    data6 <- read.table(gzfile(sprintf("%s/%s",  main_directory_ref, file6)),sep=",", header = TRUE)
    dataGSE84133_mouse <- rbind(data5, data6)
    rownames(dataGSE84133_mouse) <- dataGSE84133_mouse$X
    labelsGSE84133_mouse <- dataGSE84133_mouse$assigned_cluster
    dataGSE84133_mouse <- subset(dataGSE84133_mouse, select = -c(X, barcode, assigned_cluster))
    labelsGSE84133_mouse <- as.data.frame(labelsGSE84133_mouse)
    row.names(labelsGSE84133_mouse) <- row.names(dataGSE84133_mouse)
    colnames(labelsGSE84133_mouse) <- c("type")
    
    common_genes <- intersect(colnames(dataGSE84133_human), colnames(dataGSE84133_mouse))
    fprintf("Common genes length are [%d]\n",length(common_genes) )
    dataGSE84133_human <- dataGSE84133_human[,common_genes]
    dataGSE84133_mouse <- dataGSE84133_mouse[,common_genes]
    
    fprintf("Shape of dataset is [%d x %d]\n", dim(dataGSE84133_mouse)[1],dim(dataGSE84133_mouse)[2] )
    fprintf("Shape of labels df is [%d x %d]\n", dim(labelsGSE84133_mouse)[1],dim(labelsGSE84133_mouse)[2] )
    CHETAH_prediction_inter(dataGSE84133_human, dataGSE84133_mouse, labelsGSE84133_human, labelsGSE84133_mouse, folder_to_save = "./", dataset_name = "GSE84133_human_mouse")
    
} else if(dataset == "GSE84133_mouse_human") {
    
    dataset_ref = "GSE84133"
    main_directory_ref <- sprintf("%s/%s", folder, dataset_ref)
    file1 <- 'GSM2230757_human1_umifm_counts.csv.gz'
    data1 <- read.table(gzfile(sprintf("%s/%s",  main_directory_ref, file1)),sep=",", header = TRUE)
    file2 <- 'GSM2230758_human2_umifm_counts.csv.gz'
    data2 <- read.table(gzfile(sprintf("%s/%s",  main_directory_ref, file2)),sep=",", header = TRUE)
    file3 <- 'GSM2230759_human3_umifm_counts.csv.gz'
    data3 <- read.table(gzfile(sprintf("%s/%s",  main_directory_ref, file3)),sep=",", header = TRUE)
    file4 <- 'GSM2230760_human4_umifm_counts.csv.gz'
    data4 <- read.table(gzfile(sprintf("%s/%s",  main_directory_ref, file4)),sep=",", header = TRUE)
    dataGSE84133_human <- rbind(data1, data2, data3, data4)
    rownames(dataGSE84133_human) <- dataGSE84133_human$X
    labelsGSE84133_human <- dataGSE84133_human$assigned_cluster
    dataGSE84133_human <- subset(dataGSE84133_human, select = -c(X, barcode, assigned_cluster))
    labelsGSE84133_human <- as.data.frame(labelsGSE84133_human)
    row.names(labelsGSE84133_human) <- row.names(dataGSE84133_human)
    colnames(labelsGSE84133_human) <- c("type")
    fprintf("Shape of dataset is [%d x %d]\n", dim(dataGSE84133_human)[1],dim(dataGSE84133_human)[2] )
    fprintf("Shape of labels df is [%d x %d]\n", dim(labelsGSE84133_human)[1],dim(labelsGSE84133_human)[2] )
    
    file5 <- 'GSM2230761_mouse1_umifm_counts.csv.gz'
    data5 <- read.table(gzfile(sprintf("%s/%s",  main_directory_ref, file5)),sep=",", header = TRUE)
    file6 <- 'GSM2230762_mouse2_umifm_counts.csv.gz'
    data6 <- read.table(gzfile(sprintf("%s/%s",  main_directory_ref, file6)),sep=",", header = TRUE)
    dataGSE84133_mouse <- rbind(data5, data6)
    rownames(dataGSE84133_mouse) <- dataGSE84133_mouse$X
    labelsGSE84133_mouse <- dataGSE84133_mouse$assigned_cluster
    dataGSE84133_mouse <- subset(dataGSE84133_mouse, select = -c(X, barcode, assigned_cluster))
    labelsGSE84133_mouse <- as.data.frame(labelsGSE84133_mouse)
    row.names(labelsGSE84133_mouse) <- row.names(dataGSE84133_mouse)
    colnames(labelsGSE84133_mouse) <- c("type")
    fprintf("Shape of dataset is [%d x %d]\n", dim(dataGSE84133_mouse)[1],dim(dataGSE84133_mouse)[2] )
    fprintf("Shape of labels df is [%d x %d]\n", dim(labelsGSE84133_mouse)[1],dim(labelsGSE84133_mouse)[2] )
    CHETAH_prediction_inter(dataGSE84133_mouse, dataGSE84133_human, labelsGSE84133_mouse, labelsGSE84133_human, folder_to_save = "./", dataset_name = "GSE84133_mouse_human")
} else if(dataset == "GSE124952"){
    main_directory <- sprintf("%s/%s", folder, dataset)
    file <- 'GSE124952_expression_matrix.csv.gz'
    # read count data
    dataGSE124952 <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep=",", header = TRUE)
    dataGSE124952 <- dataGSE124952[dataGSE124952$X != "NA",] # Remove Gene/feature whose name is
    dataGSE124952 <- na.omit(dataGSE124952)
    rownames(dataGSE124952) <- dataGSE124952$X
    dataGSE124952 <- subset(dataGSE124952, select = -c(X))
    dataGSE124952 <- dataGSE124952[,sort(colnames(dataGSE124952))]
    dataGSE124952 <- as.data.frame(t(dataGSE124952))
    
    label_file <- "GSE124952_meta_data.csv.gz"
    labels <- read.table(sprintf("%s/%s", main_directory, label_file),sep=",", header = TRUE)
    rownames(labels) <- labels$X
    labels <- labels[sort(rownames(labels)), ]
    rownames(labels) <- labels$X
    labels <- subset(labels, select = c(CellType))
    colnames(labels) <- c("type")
    fprintf("Shape of dataset is [%d x %d]\n", dim(dataGSE124952)[1],dim(dataGSE124952)[2] )
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    CHETAH_prediction_intra(dataGSE124952, labels, folder_to_save = "./", dataset_name = "GSE124952")
} else if(dataset == "GSE103892"){
    main_directory <- sprintf("%s/%s", folder, dataset)
    file <- 'GSE103892_Expression_Count_Matrix.txt.gz'
    # read count data
    dataGSE103892 <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep="\t", header = TRUE)
    dataGSE103892 <- dataGSE103892[dataGSE103892$Gene != "NA",] # Remove Gene/feature whose name is
    dataGSE103892 <- na.omit(dataGSE103892)
    rownames(dataGSE103892) <- dataGSE103892$Gene
    dataGSE103892 <- subset(dataGSE103892, select = -c(Gene))
    dataGSE103892 <- dataGSE103892[,sort(colnames(dataGSE103892))]
    dataGSE103892 <- as.data.frame(t(dataGSE103892))
    
    label_file <- "GSE103892_Sample_Cell_Cluster_Information.txt.gz"
    labels <- read.table(sprintf("%s/%s", main_directory, label_file),sep="\t", header = TRUE)
    rownames(labels) <- labels$sample_cellbarcode
    labels <- labels[sort(rownames(labels)), ]
    labels <- subset(labels, select = c(cell.type))
    colnames(labels) <- c("type")
    fprintf("Shape of dataset is [%d x %d]\n", dim(dataGSE103892)[1],dim(dataGSE103892)[2] )
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    CHETAH_prediction_intra(dataGSE103892, labels, folder_to_save = "./", dataset_name = "GSE103892")
} else if(dataset == 'xin') {
    main_directory <- sprintf("%s/%s", folder, 'ADClust_datasets')
    fprintf("%s\n",main_directory)
    file <- 'Xin_normalized_counts.csv'
    labelfile <- 'Xin_normalized_labels.csv'
    data <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep=",", header = TRUE)
    rownames(data) <- data$X
    data <-  subset(data, select = -c(X))
    fprintf("Shape of dataset is [%d x %d]\n", dim(data)[1],dim(data)[2] )
    labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = ",", header = TRUE)
    rownames(labels) <- labels$X
    labels <-  subset(labels, select = c(CellType) )
    colnames(labels) <- c("type")
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    stopifnot(identical(rownames(data), rownames(labels)))
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = "Xin")
} else if(dataset == 'TM') {
    main_directory <- sprintf("%s/%s", folder, 'ADClust_datasets')
    fprintf("%s\n",main_directory)
    file <- 'TM_normalized_counts.csv'
    labelfile <- 'TM_normalized_labels.csv'
    data <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep=",", header = TRUE)
    rownames(data) <- data$X
    data <-  subset(data, select = -c(X))
    fprintf("Shape of dataset is [%d x %d]\n", dim(data)[1],dim(data)[2] )
    labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = ",", header = TRUE)
    rownames(labels) <- labels$X
    labels <-  subset(labels, select = c(CellType) )
    colnames(labels) <- c("type")
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    stopifnot(identical(rownames(data), rownames(labels)))
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = "TM")
} else if(dataset == 'Tasic') {
    main_directory <- sprintf("%s/%s", folder, 'ADClust_datasets')
    fprintf("%s\n",main_directory)
    file <- 'Tasic_normalized_counts.csv'
    labelfile <- 'Tasic_normalized_labels.csv'
    data <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep=",", header = TRUE)
    rownames(data) <- data$X
    data <-  subset(data, select = -c(X))
    fprintf("Shape of dataset is [%d x %d]\n", dim(data)[1],dim(data)[2] )
    labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = ",", header = TRUE)
    rownames(labels) <- labels$X
    labels <-  subset(labels, select = c(CellType) )
    colnames(labels) <- c("type")
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    stopifnot(identical(rownames(data), rownames(labels)))
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = "Tasic")
} else if(dataset == 'Segerstolpe') {
    main_directory <- sprintf("%s/%s", folder, 'ADClust_datasets')
    fprintf("%s\n",main_directory)
    file <- 'Segerstolpe_normalized_counts.csv'
    labelfile <- 'Segerstolpe_normalized_labels.csv'
    data <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep=",", header = TRUE)
    rownames(data) <- data$X
    data <-  subset(data, select = -c(X))
    fprintf("Shape of dataset is [%d x %d]\n", dim(data)[1],dim(data)[2] )
    labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = ",", header = TRUE)
    rownames(labels) <- labels$X
    labels <-  subset(labels, select = c(CellType) )
    colnames(labels) <- c("type")
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    stopifnot(identical(rownames(data), rownames(labels)))
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = "Segerstolpe")
} else if(dataset == 'Romanov') {
    main_directory <- sprintf("%s/%s", folder, 'ADClust_datasets')
    fprintf("%s\n",main_directory)
    file <- 'Romanov_normalized_counts.csv'
    labelfile <- 'Romanov_normalized_labels.csv'
    data <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep=",", header = TRUE)
    rownames(data) <- data$X
    data <-  subset(data, select = -c(X))
    fprintf("Shape of dataset is [%d x %d]\n", dim(data)[1],dim(data)[2] )
    labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = ",", header = TRUE)
    rownames(labels) <- labels$X
    labels <-  subset(labels, select = c(CellType) )
    colnames(labels) <- c("type")
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    stopifnot(identical(rownames(data), rownames(labels)))
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = "Romanov")
} else if(dataset == 'PBMC_68k') {
    main_directory <- sprintf("%s/%s", folder, 'ADClust_datasets')
    fprintf("%s\n",main_directory)
    file <- 'PBMC_68k_normalized_counts.csv'
    labelfile <- 'PBMC_68k_normalized_labels.csv'
    data <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep=",", header = TRUE)
    rownames(data) <- data$X
    data <-  subset(data, select = -c(X))
    fprintf("Shape of dataset is [%d x %d]\n", dim(data)[1],dim(data)[2] )
    labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = ",", header = TRUE)
    rownames(labels) <- labels$X
    labels <-  subset(labels, select = c(CellType) )
    colnames(labels) <- c("type")
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    stopifnot(identical(rownames(data), rownames(labels)))
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = "PBMC_68k")
} else if(dataset == 'Mouse_retina') {
    main_directory <- sprintf("%s/%s", folder, 'ADClust_datasets')
    fprintf("%s\n",main_directory)
    file <- 'Mouse_retina_normalized_counts.csv'
    labelfile <- 'Mouse_retina_normalized_labels.csv'
    data <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep=",", header = TRUE)
    rownames(data) <- data$X
    data <-  subset(data, select = -c(X))
    fprintf("Shape of dataset is [%d x %d]\n", dim(data)[1],dim(data)[2] )
    labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = ",", header = TRUE)
    rownames(labels) <- labels$X
    labels <-  subset(labels, select = c(CellType) )
    colnames(labels) <- c("type")
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    stopifnot(identical(rownames(data), rownames(labels)))
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = "Mouse_retina")
} else if(dataset == 'Klein') {
    main_directory <- sprintf("%s/%s", folder, 'ADClust_datasets')
    fprintf("%s\n",main_directory)
    file <- 'Klein_normalized_counts.csv'
    labelfile <- 'Klein_normalized_labels.csv'
    data <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep=",", header = TRUE)
    rownames(data) <- data$X
    data <-  subset(data, select = -c(X))
    fprintf("Shape of dataset is [%d x %d]\n", dim(data)[1],dim(data)[2] )
    labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = ",", header = TRUE)
    rownames(labels) <- labels$X
    labels <-  subset(labels, select = c(CellType) )
    colnames(labels) <- c("type")
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    stopifnot(identical(rownames(data), rownames(labels)))
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = "Klein")
} else if (dataset == 'villani_mgh_covid19'){
    main_directory <- sprintf("%s/%s", folder, 'Covid_cell_atlas')
    fprintf("%s\n",main_directory)
    file <- 'villani_mgh_covid19_counts.csv.gz'
    labelfile <- 'villani_mgh_covid19_labels.csv.gz'
    data <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep=",", header = TRUE)
    rownames(data) <- data$cell
    data <-  subset(data, select = -c(cell))
    fprintf("Shape of dataset is [%d x %d]\n", dim(data)[1],dim(data)[2] )
    labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = ",", header = TRUE)
    rownames(labels) <- labels$cell
    labels <-  subset(labels, select = c(CellType) )
    colnames(labels) <- c("type")
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    stopifnot(identical(rownames(data), rownames(labels)))
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = dataset)   
} else if (dataset == 'vento_pbmc_covid19'){
    main_directory <- sprintf("%s/%s", folder, 'Covid_cell_atlas')
    fprintf("%s\n",main_directory)
    file <- 'vento_pbmc_covid19_counts.csv.gz'
    labelfile <- 'vento_pbmc_covid19_labels.csv.gz'
    data <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep=",", header = TRUE)
    rownames(data) <- data$cell
    data <-  subset(data, select = -c(cell))
    fprintf("Shape of dataset is [%d x %d]\n", dim(data)[1],dim(data)[2] )
    labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = ",", header = TRUE)
    rownames(labels) <- labels$cell
    labels <-  subset(labels, select = c(CellType) )
    colnames(labels) <- c("type")
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    stopifnot(identical(rownames(data), rownames(labels)))
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = dataset)   
} else if (dataset == 'shalek_nasal_epithelia_covid19'){
    main_directory <- sprintf("%s/%s", folder, 'Covid_cell_atlas')
    fprintf("%s\n",main_directory)
    file <- 'shalek_nasal_epithelia_covid19_counts.csv.gz'
    labelfile <- 'shalek_nasal_epithelia_covid19_labels.csv.gz'
    data <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep=",", header = TRUE)
    rownames(data) <- data$cell
    data <-  subset(data, select = -c(cell))
    fprintf("Shape of dataset is [%d x %d]\n", dim(data)[1],dim(data)[2] )
    labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = ",", header = TRUE)
    rownames(labels) <- labels$cell
    labels <-  subset(labels, select = c(CellType) )
    colnames(labels) <- c("type")
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    stopifnot(identical(rownames(data), rownames(labels)))
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = dataset)   
} else if (dataset == 'blish_pbmc_covid19'){
    main_directory <- sprintf("%s/%s", folder, 'Covid_cell_atlas')
    fprintf("%s\n",main_directory)
    file <- 'blish_pbmc_covid19_counts.csv.gz'
    labelfile <- 'blish_pbmc_covid19_labels.csv.gz'
    data <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep=",", header = TRUE)
    rownames(data) <- data$cell
    data <-  subset(data, select = -c(cell))
    fprintf("Shape of dataset is [%d x %d]\n", dim(data)[1],dim(data)[2] )
    labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = ",", header = TRUE)
    rownames(labels) <- labels$cell
    labels <-  subset(labels, select = c(CellType) )
    colnames(labels) <- c("type")
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    stopifnot(identical(rownames(data), rownames(labels)))
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = dataset)   
} else if (dataset == 'meyer_airway_covid19'){
    main_directory <- sprintf("%s/%s", folder, 'Covid_cell_atlas')
    fprintf("%s\n",main_directory)
    file <- 'meyer_airway_covid19_counts.csv.gz'
    labelfile <- 'meyer_airway_covid19_labels.csv.gz'
    data <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep=",", header = TRUE)
    rownames(data) <- data$cell
    data <-  subset(data, select = -c(cell))
    fprintf("Shape of dataset is [%d x %d]\n", dim(data)[1],dim(data)[2] )
    labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = ",", header = TRUE)
    rownames(labels) <- labels$cell
    labels <-  subset(labels, select = c(CellType) )
    colnames(labels) <- c("type")
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    stopifnot(identical(rownames(data), rownames(labels)))
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = dataset)   
} else if (dataset == 'meyer_pbmc_covid19'){
    main_directory <- sprintf("%s/%s", folder, 'Covid_cell_atlas')
    fprintf("%s\n",main_directory)
    file <- 'meyer_pbmc_covid19_counts.csv.gz'
    labelfile <- 'meyer_pbmc_covid19_labels.csv.gz'
    data <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep=",", header = TRUE)
    rownames(data) <- data$cell
    data <-  subset(data, select = -c(cell))
    fprintf("Shape of dataset is [%d x %d]\n", dim(data)[1],dim(data)[2] )
    labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = ",", header = TRUE)
    rownames(labels) <- labels$cell
    labels <-  subset(labels, select = c(CellType) )
    colnames(labels) <- c("type")
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    stopifnot(identical(rownames(data), rownames(labels)))
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = dataset)   
} else if (id >= 26 && id <= 57){
    main_directory <- sprintf("%s/%s", folder, 'Covid_cell_atlas')
    fprintf("%s\n",main_directory)
    file <- sprintf('%s_counts.csv.gz', dataset)
    labelfile <- sprintf('%s_labels.csv.gz', dataset)
    data <- read.table(gzfile(sprintf("%s/%s", main_directory, file)),sep=",", header = TRUE)
    rownames(data) <- data$cell
    data <-  subset(data, select = -c(cell))
    fprintf("Shape of dataset is [%d x %d]\n", dim(data)[1],dim(data)[2] )
    labels <- read.table(sprintf("%s/%s", main_directory, labelfile), sep = ",", header = TRUE)
    rownames(labels) <- labels$cell
    labels <-  subset(labels, select = c(CellType) )
    colnames(labels) <- c("type")
    fprintf("Shape of labels df is [%d x %d]\n", dim(labels)[1],dim(labels)[2] )
    stopifnot(identical(rownames(data), rownames(labels)))
    CHETAH_prediction_intra(data, labels, folder_to_save = "./", dataset_name = dataset)   
} 
