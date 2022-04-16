#!/bin/bash

DatasetArray=('Segerstolpe' 'Romanov' 'PBMC_68k' 'Mouse_retina' 'Klein' 'villani_mgh_covid19' 'vento_pbmc_covid19' 'shalek_nasal_epithelia_covid19' 'blish_pbmc_covid19' 'meyer_airway_covid19' 'meyer_pbmc_covid19' 'vandamme_bronchoalveolar_lavage_fluid_covid19' 'aldinger_fetal_cerebellum' 'popescu19_fetal_liver'  'litvinukova20_heart' 'madissoon_oesophagus' 'vento_placenta_10X' 'martin_ileum' 'park_fetal_thymus' 'guo_testis' 'warner_salivery_gland' 'madissoon_spleen' 'miller_fetal_lung' 'vento_placenta_SmartSeq2' 'habib_brain' 'voigt_retina' 'james_colon_immune' 'stewart_kidney' 'vallier_restricted_gallbladder' 'menon_retina10xV3' 'byrd_gingiva' 'wang_rectum' 'lukowski_retina10xV2' 'Wyler_LungCellLine_Calu3' 'wang_colon' 'henry_prostate' 'lako_cornea' 'Wyler_LungCellLine_H1299' 'cheng18_skin' 'wang_ileum' 'smillie_colon_epithelial' 'macparland18_adult_liver')


for dataset in ${DatasetArray[*]}; do
	mkdir $dataset

done
