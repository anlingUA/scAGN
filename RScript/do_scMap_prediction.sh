#!/bin/bash
#Rscript do_scMAP.R -d 1 -f scData -s 0.75 -t 50 > scMap_output1.log
#Rscript do_scMAP.R -d 2 -f scData -s 0.75 -t 50 &> scMap_output2.log
#sleep 5
#Rscript do_scMAP.R -d 3 -f scData -s 0.75 -t 50 &> scMap_output3.log
#sleep 5
#Rscript do_scMAP.R -d 4 -f scData -s 0.75 -t 50 &> scMap_output4.log
#sleep 5
#Rscript do_scMAP.R -d 5 -f scData -s 0.75 -t 50 &> scMap_output5.log
#sleep 5
#Rscript do_scMAP.R -d 6 -f scData -s 0.75 -t 50 &> scMap_output6.log
#sleep 5
#Rscript do_scMAP.R -d 8 -f scData -s 0.75 -t 50 &> scMap_output8.log
#sleep 5
#Rscript do_scMAP.R -d 9 -f scData -s 0.75 -t 10 &> scMap_output9.log
#sleep 5

#Rscript do_scMAP.R -d 21 -f scData -s 0.75 -t 10 &> scMap_output21.log &
#Rscript do_scMAP.R -d 22 -f scData -s 0.75 -t 10 &> scMap_output22.log &
#Rscript do_scMAP.R -d 23 -f scData -s 0.75 -t 10 &> scMap_output23.log &
#Rscript do_scMAP.R -d 24 -f scData -s 0.75 -t 10 &> scMap_output24.log &
#Rscript do_scMAP.R -d 25 -f scData -s 0.75 -t 10 &> scMap_output25.log &
#Rscript do_scMAP.R -d 26 -f scData -s 0.75 -t 10 &> scMap_output26.log &

#for i in {2..20}
#for i in {21..40}
#for i in {10..57}
for i in {25..57}
do
   Rscript do_scMAP.R -d $i -f scData -s 0.75 -t 10 &> scMap_output$i.log
   sleep 5
done


