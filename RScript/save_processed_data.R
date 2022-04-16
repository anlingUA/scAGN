suppressMessages(library(Seurat)); suppressMessages(library(entropy))
library(pracma)
source('pre_process.R')
source('GenerateGraph.R')

#' This functions takes raw counts and labels of reference/query set to generate scGCN training input
#' @param count.list list of reference data and query data; rows are genes and columns are cells
#' @param label.list list of reference label and query label (if any), both are data frames with rownames identical with colnames of data; the first column is cell type
#' @return This function returns files saved in folders "input" & "process_data"
#' @export: all files are saved in current path
#' @examples: load count.list and label.list from folder "example_data"
#' save_processed_data(count.list,label.list)
save_processed_data <- function(count.list,label.list,Rgraph=TRUE,check_unknown=FALSE, folder_name='input')
{
    fprintf("\nEntering saved_process_data...\n")
    count.list <- pre_process(count_list=count.list,label_list=label.list)
    #' save counts data to certain path: 'input'
    dir.create(folder_name); 
    write.csv(t(count.list[[1]]),file=sprintf('%s/Data1.csv', folder_name),quote=F,row.names=T)
   
    if (Rgraph)
    {
        #' use R generated graph
        new.dat1 <- count.list[[1]];
        new.lab1 <- label.list[[1]];
        fprintf("\nCalculating Graph...\n")
        graphs <- suppressWarnings(GenerateGraph(Dat1=new.dat1,
                                                 Lab1=new.lab1,K=5,
                                                 check.unknown=check_unknown))
        write.csv(graphs[[1]],file=sprintf('%s/intra_graph.csv', folder_name),quote=F,row.names=T)
        write.csv(new.lab1,file=sprintf('%s/label1.csv', folder_name),quote=F,row.names=T)
    }
    fprintf("Exiting saved_process_data...\n")
}
