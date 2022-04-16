suppressMessages(library(Seurat)); suppressMessages(library(entropy))
source('select_feature.R')

#' This functions takes refrence data and labels to identify variable gene features
#' @param count_list list of reference data and query data; rows are genes and columns are cells
#' @param label_list list of reference label and query label (if any), both are data frames with rownames identical with colnames of data; the first column is cell type
#' @return This function returns four elements, including the normalized data, scaled data, hvg features, and selected features
#' @export
#' @examples
#' pre_process(count_list,label_list)
pre_process <- function(count_list,label_list)
{
    fprintf("\nPre-processing data...\n")
    sel.features <- select_feature(count_list[[1]],label_list[[1]])
    count_list_new <- list(count_list[[1]][sel.features,])
    fprintf("\nPre-processing Done...\n")
    return (count_list_new)
}