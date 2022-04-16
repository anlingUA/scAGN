#!/bin/bash
Rscript generate_data_for_CCA.R -d 1 -f scData -s 0.75 -t 50 &> output1.log &
Rscript generate_data_for_CCA.R -d 2 -f scData -s 0.75 -t 50 &> output2.log &
Rscript generate_data_for_CCA.R -d 3 -f scData -s 0.75 -t 50 &> output3.log &
Rscript generate_data_for_CCA.R -d 4 -f scData -s 0.75 -t 50 &> output4.log &
Rscript generate_data_for_CCA.R -d 5 -f scData -s 0.75 -t 50 &> output5.log &
Rscript generate_data_for_CCA.R -d 6 -f scData -s 0.75 -t 50 &> output6.log &
Rscript generate_data_for_CCA.R -d 8 -f scData -s 0.75 -t 50 &> output8.log &
Rscript generate_data_for_CCA.R -d 9 -f scData -s 0.75 -t 10 &> output9.log &
Rscript generate_data_for_CCA.R -d 21 -f scData -s 0.75 -t 10 &> output21.log &
Rscript generate_data_for_CCA.R -d 22 -f scData -s 0.75 -t 10 &> output22.log &
Rscript generate_data_for_CCA.R -d 23 -f scData -s 0.75 -t 10 &> output23.log &
Rscript generate_data_for_CCA.R -d 24 -f scData -s 0.75 -t 10 &> output24.log &
Rscript generate_data_for_CCA.R -d 25 -f scData -s 0.75 -t 10 &> output25.log &
Rscript generate_data_for_CCA.R -d 26 -f scData -s 0.75 -t 10 &> output26.log &

for i in {26..57}
do
   Rscript generate_data_for_CCA.R -d $i -f scData -s 0.75 -t 10 &> output$i.log &
done


