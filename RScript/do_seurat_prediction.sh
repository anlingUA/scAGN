#!/bin/bash
#Rscript do_seuratv3.R -d 1 -f scData -s 0.75 -t 50 > Seurat_output1.log
#Rscript do_seuratv3.R -d 2 -f scData -s 0.75 -t 50 &> Seurat_output2.log
#sleep 5
#Rscript do_seuratv3.R -d 3 -f scData -s 0.75 -t 50 &> Seurat_output3.log
#sleep 5
#Rscript do_seuratv3.R -d 4 -f scData -s 0.75 -t 50 &> Seurat_output4.log
#sleep 5
#Rscript do_seuratv3.R -d 5 -f scData -s 0.75 -t 50 &> Seurat_output5.log
#sleep 5
#Rscript do_seuratv3.R -d 6 -f scData -s 0.75 -t 50 &> Seurat_output6.log
#sleep 5
#Rscript do_seuratv3.R -d 8 -f scData -s 0.75 -t 50 &> Seurat_output8.log
#sleep 5
#Rscript do_seuratv3.R -d 9 -f scData -s 0.75 -t 10 &> Seurat_output9.log
#sleep 5

#Rscript do_seuratv3.R -d 21 -f scData -s 0.75 -t 10 &> Seurat_output21.log &
#Rscript do_seuratv3.R -d 22 -f scData -s 0.75 -t 10 &> Seurat_output22.log &
#Rscript do_seuratv3.R -d 23 -f scData -s 0.75 -t 10 &> Seurat_output23.log &
#Rscript do_seuratv3.R -d 24 -f scData -s 0.75 -t 10 &> Seurat_output24.log &
#Rscript do_seuratv3.R -d 25 -f scData -s 0.75 -t 10 &> Seurat_output25.log &
#Rscript do_seuratv3.R -d 26 -f scData -s 0.75 -t 10 &> Seurat_output26.log &

#for i in {10..57}
#for i in {13..57}
for i in {18..57}
do
   Rscript do_seuratv3.R -d $i -f scData -s 0.75 -t 10 &> Seurat_output$i.log
done


