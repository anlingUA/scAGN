#' This functions provides the optional way to generate scGCN labels
#' @param Dat1  reference data; rows are genes and columns are cells
#' @param Dat2  query data; rows are genes and columns are cells
#' @param Lab1  label of reference data; rows are genes and columns are cells
#' @param Lab2  query data; rows are genes and columns are cells

#' @return This function returns files saved in folders "input" & "process_data"
#' @export: all files are saved in current path
#' @examples: load count.list and label.list from folder "example_data"
#' save_processed_data(count.list,label.list)

suppressMessages(library(Seurat)); suppressMessages(library(entropy))
DataTransform <- function(obj)
{
    obj <- NormalizeData(obj,verbose=F)
    obj <- FindVariableFeatures(obj,
                                selection.method = "vst",
                                nfeatures = 2000,verbose=F)
    obj <- ScaleData(obj,features=rownames(obj),verbose=FALSE)
    obj <- RunPCA(obj, features=rownames(obj), verbose = FALSE)
    return(obj)
}

GenerateGraph <- function(Dat1,Dat2,Lab1,K,check.unknown)
{
    object1 <- CreateSeuratObject(counts=Dat1,project = "1",assay = "Data1",
                                  min.cells = 0,min.features = 0,
                                  names.field = 1,names.delim = "_")
    objects <- list(object1)    
    objects1 <- lapply(objects, DataTransform)

    #'  Intra-data graph  
    d2.list <- list(objects1[[1]],objects1[[1]])
    d2.nn <- FindIntegrationAnchors(object.list =d2.list,k.anchor=K,verbose=F)    
    d2.arc=d2.nn@anchors
    d2.arc1=cbind(d2.arc[d2.arc[,4]==1,1],d2.arc[d2.arc[,4]==1,2],d2.arc[d2.arc[,4]==1,3])
    d2.grp=d2.arc1[d2.arc1[,3]>0,1:2]-1
    final <- list(intraG=d2.grp)
    fprintf("\nGraph Calculation Done...\n")
    return (final)
}


