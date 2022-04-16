#!/bin/bash
#Rscript do_SingleR.R -d 1 -f scData -s 0.75 -t 50 >  singleR_output1.log
#Rscript do_SingleR.R -d 2 -f scData -s 0.75 -t 50 &>  singleR_output2.log
#sleep 5
#Rscript do_SingleR.R -d 3 -f scData -s 0.75 -t 50 &>  singleR_output3.log
#sleep 5
#Rscript do_SingleR.R -d 4 -f scData -s 0.75 -t 50 &>  singleR_output4.log
#sleep 5
#Rscript do_SingleR.R -d 5 -f scData -s 0.75 -t 50 &>  singleR_output5.log
#sleep 5
#Rscript do_SingleR.R -d 6 -f scData -s 0.75 -t 50 &>  singleR_output6.log
#sleep 5
#Rscript do_SingleR.R -d 8 -f scData -s 0.75 -t 50 &>  singleR_output8.log
#sleep 5
#Rscript do_SingleR.R -d 9 -f scData -s 0.75 -t 10 &>  singleR_output9.log
#sleep 5

#Rscript do_SingleR.R -d 21 -f scData -s 0.75 -t 10 &>  singleR_output21.log &
#Rscript do_SingleR.R -d 22 -f scData -s 0.75 -t 10 &>  singleR_output22.log &
#Rscript do_SingleR.R -d 23 -f scData -s 0.75 -t 10 &>  singleR_output23.log &
#Rscript do_SingleR.R -d 24 -f scData -s 0.75 -t 10 &>  singleR_output24.log &
#Rscript do_SingleR.R -d 25 -f scData -s 0.75 -t 10 &>  singleR_output25.log &
#Rscript do_SingleR.R -d 26 -f scData -s 0.75 -t 10 &>  singleR_output26.log &

#for i in {52..57}
#for i in {10..15}
#for i in {16..32}
#for i in {32..40}
# skipped 37
#for i in {38..40}
#for i in {41..45}
#for i in {46..50}
#for i in {51..52}
#for i in {10..57}
#for i in {34..57}
for i in {39..57}
do
   Rscript do_SingleR.R -d $i -f scData -s 0.75 -t 10 &> singleR_output$i.log
   sleep 5
done


