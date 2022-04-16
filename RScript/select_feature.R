suppressMessages(library(Seurat)); suppressMessages(library(entropy))

#' This functions takes refrence data and labels to identify variable gene features
#' @param data reference data; rows are genes and columns are cells
#' @param label data frame with rownames identical with colnames of data; the first column is cell type
#' @param nf number of variable features
#' @return This function returns variable gene features using ANOVA
#' @export
#' @examples
#' select_feature(data=reference.data,label=reference.label)
select_feature <- function(data,label,nf=2000)
{
    fprintf("\nEntering select_feature...\n")
    M <- nrow(data); new.label <- label[,1]
    pv1 <- sapply(1:M, function(i){
        mydataframe <- data.frame(y=as.numeric(data[i,]), ig=new.label)
        fit <- aov(y ~ ig, data=mydataframe)
        summary(fit)[[1]][["Pr(>F)"]][1]})
    names(pv1) <- rownames(data)
    pv1.sig <- names(pv1)[order(pv1)[1:nf]]
    egen <- unique(pv1.sig)
    fprintf("\nFeature selection done...\n")
    return (egen)
}