#!/usr/bin/env python

# Initial Date: April 2021
# Author: Rahul Bhadani
# Copyright (c)  Rahul Bhadani
# All rights reserved.

__author__ = 'Rahul Bhadani'
__email__  = 'rahulbhadani@email.arizona.edu'

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import seaborn as sea
import sklearn.metrics
from sklearn.metrics import confusion_matrix

datasets_to_include = ['GSE108989', 'GSE72056', 'GSE85241', 'GSE98638', 'GSE99254',
 'GSE118389', 'GSE84133', 'GSE165092', 'GSE83139', 'GSE85241_GSE81076',
'GSE81076_GSE85241', 'GSE124952', 'Xin', 'TM', 'Tasic',
'Segerstolpe', 'Romanov', 'PBMC_68k', 'Mouse_retina', 'Klein',
'villani_mgh_covid19', 'vento_pbmc_covid19', 'shalek_nasal_epithelia_covid19', 'blish_pbmc_covid19', 'meyer_airway_covid19',
'meyer_pbmc_covid19', 'vandamme_bronchoalveolar_lavage_fluid_covid19', 'aldinger_fetal_cerebellum', 'popescu19_fetal_liver',  'litvinukova20_heart',
'madissoon_oesophagus', 'vento_placenta_10X', 'martin_ileum', 'park_fetal_thymus', 'guo_testis',
'warner_salivery_gland', 'madissoon_spleen', 'miller_fetal_lung', 'vento_placenta_SmartSeq2','habib_brain',
'voigt_retina', 'james_colon_immune', 'stewart_kidney', 'vallier_restricted_gallbladder','menon_retina10xV3',
'byrd_gingiva', 'wang_rectum', 'lukowski_retina10xV2', 'Wyler_LungCellLine_Calu3','wang_colon',
'henry_prostate', 'lako_cornea', 'Wyler_LungCellLine_H1299', 'cheng18_skin','wang_ileum',
'smillie_colon_epithelial', 'macparland18_adult_liver',
'GSE83139_GSE84133', 'GSE84133_GSE83139', 'GSE103892']

# datasets_to_include = ['GSE108989', 'GSE72056', 'GSE85241', 'GSE98638', 'GSE99254',
#  'GSE118389', 'GSE84133', 'GSE165092', 'GSE83139', 'GSE85241_GSE81076',
# 'GSE81076_GSE85241', 'GSE124952', 'Xin', 'TM', 'Tasic', 'Romanov', 'PBMC_68k', 'Mouse_retina', 'Klein',
# 'villani_mgh_covid19', 'vento_pbmc_covid19', 'shalek_nasal_epithelia_covid19', 'blish_pbmc_covid19', 'meyer_airway_covid19',
# 'meyer_pbmc_covid19', 'vandamme_bronchoalveolar_lavage_fluid_covid19', 'aldinger_fetal_cerebellum', 'popescu19_fetal_liver',  'litvinukova20_heart',
# 'madissoon_oesophagus', 'vento_placenta_10X', 'martin_ileum', 'park_fetal_thymus', 'guo_testis',
# 'warner_salivery_gland', 'madissoon_spleen', 'miller_fetal_lung', 'vento_placenta_SmartSeq2','habib_brain',
# 'voigt_retina', 'james_colon_immune', 'stewart_kidney', 'vallier_restricted_gallbladder','menon_retina10xV3',
# 'byrd_gingiva', 'wang_rectum', 'lukowski_retina10xV2', 'Wyler_LungCellLine_Calu3','wang_colon',
# 'henry_prostate', 'lako_cornea', 'Wyler_LungCellLine_H1299', 'cheng18_skin','wang_ileum',
# 'smillie_colon_epithelial', 'macparland18_adult_liver',
# 'GSE83139_GSE84133', 'GSE84133_GSE83139', 'GSE103892']


print(datasets_to_include)

## Calculate AGN's data


AGN_output_folder = "../data"
df_metric = pd.DataFrame(columns=['dataset', 'prediction_accuracy', 'precision_score', 'recall_score', 'f1_score'])
runtime_df = pd.DataFrame()



df_pred_accuracy = pd.read_csv('Predication_Comparison2.csv')


df_pred_accuracy = df_pred_accuracy.rename(columns={'Unnamed: 0': 'Dataset'})
df_pred_accuracy = df_pred_accuracy.set_index('Dataset')

print(df_pred_accuracy)


# df_pred_accuracy.drop(columns=['Dataset'], inplace=True)
data1 = ['TM', 'Tasic', 'Mouse_retina', 'PBMC_68k']
df1 = df_pred_accuracy.filter(items = data1, axis=0)
print(df1)
sea.set_context('paper')
sea.set_style("dark")
fig, ax = plt.subplots(figsize=(7, 4))
ax.set_ylim([0, 1.35])
df1[['scAGN_prediction', 'seurat_prediction', 'singleR_prediction', 'scmap-cluster_prediction', 'scmap-cell_prediction', 'chetah_prediction']].plot(kind="bar", ax = ax, rot=45)
x_offset = -0.03
y_offset = 0.06
for p in ax.patches:
    #ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    b = p.get_bbox()
    val = "{:.2f}".format(b.y1 + b.y0)        
    ax.annotate(val, ((b.x0 + b.x1)/2 + x_offset, b.y1 + y_offset), rotation=90)

plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('Prediction Accuracy: Healthy Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_Prediction_Accuracy_Healthy.png', bbox_inches='tight')
plt.savefig('Comparison_Prediction_Accuracy_Healthy.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
ax.set_ylim([0, 1.35])
df1[['scAGN_precision_score', 'Seurat_precision_score', 'SingleR_precision_score', 'Scmap-cluster_precision_score', 'Scmap-cell_precision_score', 'Chetah_precision_score']].plot(kind="bar", ax = ax, rot=45)
x_offset = -0.03
y_offset = 0.06
for p in ax.patches:
    #ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    b = p.get_bbox()
    val = "{:.2f}".format(b.y1 + b.y0)        
    ax.annotate(val, ((b.x0 + b.x1)/2 + x_offset, b.y1 + y_offset), rotation=90)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('Precision Score: Healthy Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_Precision_Score_Healthy.png', bbox_inches='tight')
plt.savefig('Comparison_Precision_Score_Healthy.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
ax.set_ylim([0, 1.35])
df1[['scAGN_recall_score', 'Seurat_recall_score', 'SingleR_recall_score', 'Scmap-cluster_recall_score', 'Scmap-cell_recall_score', 'Chetah_recall_score']].plot(kind="bar", ax = ax, rot=45)
x_offset = -0.03
y_offset = 0.06
for p in ax.patches:
    #ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    b = p.get_bbox()
    val = "{:.2f}".format(b.y1 + b.y0)        
    ax.annotate(val, ((b.x0 + b.x1)/2 + x_offset, b.y1 + y_offset), rotation=90)

plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('Recall Score: Healthy Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_Recall_Score_Healthy.png', bbox_inches='tight')
plt.savefig('Comparison_Recall_Score_Healthy.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
ax.set_ylim([0, 1.35])
df1[['scAGN_f1_score', 'Seurat_f1_score', 'SingleR_f1_score', 'Scmap-cluster_f1_score', 'Scmap-cell_f1_score', 'Chetah_f1_score']].plot(kind="bar", ax = ax, rot=45)
x_offset = -0.03
y_offset = 0.06
for p in ax.patches:
    #ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    b = p.get_bbox()
    val = "{:.2f}".format(b.y1 + b.y0)        
    ax.annotate(val, ((b.x0 + b.x1)/2 + x_offset, b.y1 + y_offset), rotation=90)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('F1 Score: Healthy Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_F1_Score_Healthy.png', bbox_inches='tight')
plt.savefig('Comparison_F1_Score_Healthy.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
ax.set_ylim([0, 1.35])
df1[['scAGN_mcc_score', 'Seurat_mcc_score', 'SingleR_mcc_score', 'Scmap-cluster_mcc_score', 'Scmap-cell_mcc_score', 'Chetah_mcc_score']].plot(kind="bar", ax = ax , rot = 45)
x_offset = -0.03
y_offset = 0.06
for p in ax.patches:
    #ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    b = p.get_bbox()
    val = "{:.2f}".format(b.y1 + b.y0)        
    ax.annotate(val, ((b.x0 + b.x1)/2 + x_offset, b.y1 + y_offset), rotation=90)

plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title("Matthew's Correlation Coefficient", fontsize=18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize=14)
plt.savefig('Comparison_mcc_Score_Healthy.png', bbox_inches='tight')
plt.savefig('Comparison_mcc_Score_Healthy.pdf', bbox_inches='tight')

data2 = ['GSE108989', 'GSE118389', 'GSE72056', 'GSE98638', 'GSE99254']
df2 = df_pred_accuracy.filter(items = data2, axis=0)
print(df2)
sea.set_context('paper')
sea.set_style("dark")
fig, ax = plt.subplots(figsize=(7, 4))
ax.set_ylim([0, 1.35])
df2[['scAGN_prediction', 'seurat_prediction', 'singleR_prediction', 'scmap-cluster_prediction', 'scmap-cell_prediction', 'chetah_prediction']].plot(kind="bar", ax = ax, rot=45)
x_offset = -0.03
y_offset = 0.06
for p in ax.patches:
    #ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    b = p.get_bbox()
    val = "{:.2f}".format(b.y1 + b.y0)        
    ax.annotate(val, ((b.x0 + b.x1)/2 + x_offset, b.y1 + y_offset), rotation=90)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('Prediction Accuracy: Cancer Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_Prediction_Accuracy_Cancer.png', bbox_inches='tight')
plt.savefig('Comparison_Prediction_Accuracy_Cancer.pdf', bbox_inches='tight')


fig, ax = plt.subplots(figsize=(7, 4))
ax.set_ylim([0, 1.35])
df2[['scAGN_precision_score', 'Seurat_precision_score', 'SingleR_precision_score', 'Scmap-cluster_precision_score', 'Scmap-cell_precision_score', 'Chetah_precision_score']].plot(kind="bar", ax = ax, rot=45)
x_offset = -0.03
y_offset = 0.06
for p in ax.patches:
    #ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    b = p.get_bbox()
    val = "{:.2f}".format(b.y1 + b.y0)        
    ax.annotate(val, ((b.x0 + b.x1)/2 + x_offset, b.y1 + y_offset), rotation=90)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('Precision Score: Cancer Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_Precision_Score_Cancer.png', bbox_inches='tight')
plt.savefig('Comparison_Precision_Score_Cancer.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
ax.set_ylim([0, 1.35])
df2[['scAGN_recall_score', 'Seurat_recall_score', 'SingleR_recall_score', 'Scmap-cluster_recall_score', 'Scmap-cell_recall_score', 'Chetah_recall_score']].plot(kind="bar", ax = ax, rot=45)
x_offset = -0.03
y_offset = 0.06
for p in ax.patches:
    #ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    b = p.get_bbox()
    val = "{:.2f}".format(b.y1 + b.y0)        
    ax.annotate(val, ((b.x0 + b.x1)/2 + x_offset, b.y1 + y_offset), rotation=90)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('Recall Score: Cancer Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_Recall_Score_Cancer.png', bbox_inches='tight')
plt.savefig('Comparison_Recall_Score_Cancer.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
ax.set_ylim([0, 1.35])
df2[['scAGN_f1_score', 'Seurat_f1_score', 'SingleR_f1_score', 'Scmap-cluster_f1_score', 'Scmap-cell_f1_score', 'Chetah_f1_score']].plot(kind="bar", ax = ax, rot=45)
x_offset = -0.03
y_offset = 0.06
for p in ax.patches:
    #ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    b = p.get_bbox()
    val = "{:.2f}".format(b.y1 + b.y0)        
    ax.annotate(val, ((b.x0 + b.x1)/2 + x_offset, b.y1 + y_offset), rotation=90)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('F1 Score: Cancer Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_F1_Score_Cancer.png', bbox_inches='tight')
plt.savefig('Comparison_F1_Score_Cancer.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
ax.set_ylim([0, 1.35])
df2[['scAGN_mcc_score', 'Seurat_mcc_score', 'SingleR_mcc_score', 'Scmap-cluster_mcc_score', 'Scmap-cell_mcc_score', 'Chetah_mcc_score']].plot(kind="bar", ax = ax , rot = 45)
x_offset = -0.03
y_offset = 0.06
for p in ax.patches:
    #ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    b = p.get_bbox()
    val = "{:.2f}".format(b.y1 + b.y0)        
    ax.annotate(val, ((b.x0 + b.x1)/2 + x_offset, b.y1 + y_offset), rotation=90)

plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title("Matthew's Correlation Coefficient", fontsize=18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize=14)
plt.savefig('Comparison_mcc_Score_Cancer.png', bbox_inches='tight')
plt.savefig('Comparison_mcc_Score_Cancer.pdf', bbox_inches='tight')

data3 = ['blish_pbmc_covid19', 'shalek_nasal_epithelia_covid19', 'villani_mgh_covid19']
df3 = df_pred_accuracy.filter(items = data3, axis=0)
print(df3)
sea.set_context('paper')
sea.set_style("dark")
fig, ax = plt.subplots(figsize=(7, 4))
ax.set_ylim([0, 1.35])
df3[['scAGN_prediction', 'seurat_prediction', 'singleR_prediction', 'scmap-cluster_prediction', 'scmap-cell_prediction', 'chetah_prediction']].plot(kind="bar", ax = ax, rot=45)
x_offset = -0.03
y_offset = 0.06
for p in ax.patches:
    #ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    b = p.get_bbox()
    val = "{:.2f}".format(b.y1 + b.y0)        
    ax.annotate(val, ((b.x0 + b.x1)/2 + x_offset, b.y1 + y_offset), rotation=90)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('Prediction Accuracy: COVID19 Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_Prediction_Accuracy_COVID19.png', bbox_inches='tight')
plt.savefig('Comparison_Prediction_Accuracy_COVID19.pdf', bbox_inches='tight')


fig, ax = plt.subplots(figsize=(7, 4))
ax.set_ylim([0, 1.35])
df3[['scAGN_precision_score', 'Seurat_precision_score', 'SingleR_precision_score', 'Scmap-cluster_precision_score', 'Scmap-cell_precision_score', 'Chetah_precision_score']].plot(kind="bar", ax = ax, rot=45)
x_offset = -0.03
y_offset = 0.06
for p in ax.patches:
    #ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    b = p.get_bbox()
    val = "{:.2f}".format(b.y1 + b.y0)        
    ax.annotate(val, ((b.x0 + b.x1)/2 + x_offset, b.y1 + y_offset), rotation=90)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('Precision Score: COVID19 Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_Precision_Score_COVID19.png', bbox_inches='tight')
plt.savefig('Comparison_Precision_Score_COVID19.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
ax.set_ylim([0, 1.35])
df3[['scAGN_recall_score', 'Seurat_recall_score', 'SingleR_recall_score', 'Scmap-cluster_recall_score', 'Scmap-cell_recall_score', 'Chetah_recall_score']].plot(kind="bar", ax = ax, rot=45)
x_offset = -0.03
y_offset = 0.06
for p in ax.patches:
    #ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    b = p.get_bbox()
    val = "{:.2f}".format(b.y1 + b.y0)        
    ax.annotate(val, ((b.x0 + b.x1)/2 + x_offset, b.y1 + y_offset), rotation=90)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('Recall Score: COVID19 Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_Recall_Score_COVID19.png', bbox_inches='tight')
plt.savefig('Comparison_Recall_Score_COVID19.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
ax.set_ylim([0, 1.35])
df3[['scAGN_f1_score', 'Seurat_f1_score', 'SingleR_f1_score', 'Scmap-cluster_f1_score', 'Scmap-cell_f1_score', 'Chetah_f1_score']].plot(kind="bar", ax = ax, rot=45)
x_offset = -0.03
y_offset = 0.06
for p in ax.patches:
    #ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    b = p.get_bbox()
    val = "{:.2f}".format(b.y1 + b.y0)        
    ax.annotate(val, ((b.x0 + b.x1)/2 + x_offset, b.y1 + y_offset), rotation=90)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('F1 Score: COVID19 Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_F1_Score_COVID19.png', bbox_inches='tight')
plt.savefig('Comparison_F1_Score_COVID19.pdf', bbox_inches='tight')


fig, ax = plt.subplots(figsize=(7, 4))
ax.set_ylim([0, 1.35])
df3[['scAGN_mcc_score', 'Seurat_mcc_score', 'SingleR_mcc_score', 'Scmap-cluster_mcc_score', 'Scmap-cell_mcc_score', 'Chetah_mcc_score']].plot(kind="bar", ax = ax , rot = 45)
x_offset = -0.03
y_offset = 0.06
for p in ax.patches:
    #ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    b = p.get_bbox()
    val = "{:.2f}".format(b.y1 + b.y0)        
    ax.annotate(val, ((b.x0 + b.x1)/2 + x_offset, b.y1 + y_offset), rotation=90)

plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title("Matthew's Correlation Coefficient", fontsize=18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize=14)
plt.savefig('Comparison_mcc_Score_COVID19.png', bbox_inches='tight')
plt.savefig('Comparison_mcc_Score_COVID19.pdf', bbox_inches='tight')