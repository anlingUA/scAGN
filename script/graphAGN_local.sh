#!/bin/bash

conda activate gcn

cd ~/gene-cn/notebooks

DatasetArray=("GSE108989" "GSE72056" "GSE85241" "GSE98638" "GSE99254" "GSE118389" "GSE84133" "GSE165092" "GSE83139" "GSE85241_GSE81076" "GSE81076_GSE85241" "GSE124952" "Xin" "TM" "Tasic" "Segerstolpe" "Romanov" "PBMC_68k" "Mouse_retina" "Klein" "villani_mgh_covid19" "vento_pbmc_covid19" "shalek_nasal_epithelia_covid19" "blish_pbmc_covid19" "meyer_airway_covid19" "meyer_pbmc_covid19" "vandamme_bronchoalveolar_lavage_fluid_covid19" "aldinger_fetal_cerebellum" "popescu19_fetal_liver" "litvinukova20_heart" "madissoon_oesophagus" "vento_placenta_10X" "martin_ileum" "park_fetal_thymus" "guo_testis" "warner_salivery_gland" "madissoon_spleen" "miller_fetal_lung" "vento_placenta_SmartSeq2" "habib_brain" "voigt_retina" "james_colon_immune" "stewart_kidney" "vallier_restricted_gallbladder" "menon_retina10xV3" "byrd_gingiva" "wang_rectum" "lukowski_retina10xV2" "Wyler_LungCellLine_Calu3" "wang_colon" "henry_prostate" "lako_cornea" "Wyler_LungCellLine_H1299" "cheng18_skin" "wang_ileum" "smillie_colon_epithelial" "macparland18_adult_liver")


for dataset in ${DatasetArray[*]}; do
    # varies the number of layers
    echo $dataset
    now=$(date +%Y_%m_%d_%H_%M_%S_%N)
    ~/gene-cn/notebooks/graph_AGN.py -i 128 -l 2 -d 0.0 -t 0.05 -r 0.1 -w 0.0000 -s 10 -G 0.8 -D $dataset -P 2000 > ~/gene-cn/slurm/GRAPH_agn_1_$dataset$now.txt 
    sleep 5
    ~/gene-cn/notebooks/graph_AGN.py -i 128 -l 3 -d 0.0 -t 0.05 -r 0.1 -w 0.0000 -s 10 -G 0.8 -D $dataset -P 2000 > ~/gene-cn/slurm/GRAPH_agn_2_$dataset$now.txt
    sleep 5
    ~/gene-cn/notebooks/graph_AGN.py -i 128 -l 4 -d 0.0 -t 0.05 -r 0.1 -w 0.0000 -s 10 -G 0.8 -D $dataset -P 2000 > ~/gene-cn/slurm/GRAPH_agn_3_$dataset$now.txt
    sleep 5
    ~/gene-cn/notebooks/graph_AGN.py -i 128 -l 5 -d 0.0 -t 0.05 -r 0.1 -w 0.0000 -s 10 -G 0.8 -D $dataset -P 2000 > ~/gene-cn/slurm/GRAPH_agn_4_$dataset$now.txt
    sleep 5

    # varies hidden units
    ~/gene-cn/notebooks/graph_AGN.py -i 32 -l 3 -d 0.0 -t 0.05 -r 0.1 -w 0.0000 -s 10 -G 0.8 -D $dataset -P 2000 > ~/gene-cn/slurm/GRAPH_agn_5_$dataset$now.txt
    sleep 5
    ~/gene-cn/notebooks/graph_AGN.py -i 64 -l 3 -d 0.0 -t 0.05 -r 0.1 -w 0.0000 -s 10 -G 0.8 -D $dataset -P 2000 > ~/gene-cn/slurm/GRAPH_agn_6_$dataset$now.txt
    sleep 5
    ~/gene-cn/notebooks/graph_AGN.py -i 128 -l 3 -d 0.0 -t 0.05 -r 0.1 -w 0.0000 -s 10 -G 0.8 -D $dataset -P 2000 > ~/gene-cn/slurm/GRAPH_agn_7_$dataset$now.txt
    sleep 5
    ~/gene-cn/notebooks/graph_AGN.py -i 256 -l 3 -d 0.0 -t 0.05 -r 0.1 -w 0.0000 -s 10 -G 0.8 -D $dataset -P 2000 > ~/gene-cn/slurm/GRAPH_agn_8_$dataset$now.txt
    sleep 5
done
