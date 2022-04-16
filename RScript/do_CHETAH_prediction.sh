#!/bin/bash
#Rscript do_CHETAH.R -d 1 -f scData -s 0.75 -t 50 > CHETAH_output1.log
Rscript do_CHETAH.R -d 2 -f scData -s 0.75 -t 50 > CHETAH_output2.log
sleep 5
Rscript do_CHETAH.R -d 3 -f scData -s 0.75 -t 50 > CHETAH_output3.log
sleep 5
Rscript do_CHETAH.R -d 4 -f scData -s 0.75 -t 50 > CHETAH_output4.log
sleep 5
Rscript do_CHETAH.R -d 5 -f scData -s 0.75 -t 50 > CHETAH_output5.log
sleep 5
Rscript do_CHETAH.R -d 6 -f scData -s 0.75 -t 50 > CHETAH_output6.log
sleep 5
Rscript do_CHETAH.R -d 8 -f scData -s 0.75 -t 50 > CHETAH_output8.log
sleep 5
Rscript do_CHETAH.R -d 9 -f scData -s 0.75 -t 10 > CHETAH_output9.log
sleep 5

#Rscript do_CHETAH.R -d 21 -f scData -s 0.75 -t 10 &> CHETAH_output21.log &
#Rscript do_CHETAH.R -d 53 -f scData -s 0.75 -t 10 &> CHETAH_output53.log &
#Rscript do_CHETAH.R -d 22 -f scData -s 0.75 -t 10 &> CHETAH_output22.log &
#Rscript do_CHETAH.R -d 23 -f scData -s 0.75 -t 10 &> CHETAH_output23.log &
#Rscript do_CHETAH.R -d 24 -f scData -s 0.75 -t 10 &> CHETAH_output24.log &
#Rscript do_CHETAH.R -d 25 -f scData -s 0.75 -t 10 &> CHETAH_output25.log &
#Rscript do_CHETAH.R -d 26 -f scData -s 0.75 -t 10 &> CHETAH_output26.log &

#for i in {2..20}
#for i in {21..40}
for i in {10..57}
do
   Rscript do_CHETAH.R -d $i -f scData -s 0.75 -t 10 > CHETAH_output$i.log
   sleep 5
done


