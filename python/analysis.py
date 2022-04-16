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

for data in datasets_to_include:
    print("======== {} =========".format(data))
    folder = AGN_output_folder+ "/" + data
    AGN_prediction_files = glob.glob(folder + "/" + "*prediction_metric.csv") 
    
    df_main = pd.DataFrame(columns=['predication_accuracy', 'precision_score', 'recall_score', 'f1_score',
       'hidden_units', 'n_layers', 'dropout', 'loss_threshold',
       'learning_rate', 'weight_decay', 'step_size', 'gamma', 'patience',
       'dataset'])

    for AGN_pred_file in AGN_prediction_files:
        splits = AGN_pred_file.split('_')
        hidden_units = int(splits[splits.index('i') + 1])
        n_layers = int(splits[splits.index('l') + 1])
        dropout = float(splits[splits.index('d') + 1])
        loss_threshold = float(splits[splits.index('t') + 1])
        learning_rate = float(splits[splits.index('r') + 1])
        weight_decay = float(splits[splits.index('w') + 1])
        step_size = int(splits[splits.index('s') + 1])
        gamma = float(splits[splits.index('G') + 1])
        patience = int(splits[splits.index('P') + 1])
        
        df = pd.read_csv(AGN_pred_file, index_col = 0)
        df = df.T
        df.index.name = None
        df = df.reset_index(drop=True)
        df['hidden_units'] = hidden_units
        df['n_layers'] = n_layers
        df['dropout'] = dropout
        df['loss_threshold'] = loss_threshold
        df['learning_rate'] = learning_rate
        df['weight_decay'] = weight_decay
        df['step_size'] = step_size
        df['gamma'] = gamma
        df['patience'] = patience
        df['dataset'] = data
        df_main = pd.concat([df_main, df])
   
    df_main2 = df_main.reset_index(drop=True)
    #print(df_main)

    if df_main.shape[0]:
        sea.set_context('paper')
        sea.set_style("dark")
        fig, ax = plt.subplots(figsize=(3, 3))
        sea.lineplot(x = 'hidden_units', y = 'predication_accuracy', data = df_main,  ax = ax, markers = True)
        sea.scatterplot(x = 'hidden_units', y = 'predication_accuracy', data = df_main,  ax = ax, markers = '.', size = 'n_layers')
        ax.set_xlabel('# Hidden Units')
        ax.set_ylabel('Prediction Accuracy')
        ax.set_title('scAGN: {}'.format(data))
        plt.tight_layout()
        plt.savefig('AGN_'+data + '_predictionVsHiddenUnits.png', bbox_inches='tight')
        plt.savefig('AGN_'+data + '_predictionVsHiddenUnits.pdf', bbox_inches='tight')
        
        sea.set_context('paper')
        sea.set_style("dark")
        fig, ax = plt.subplots(figsize=(3, 3))
        sea.lineplot(x = 'hidden_units', y = 'precision_score', data = df_main,  ax = ax, markers = True)
        sea.scatterplot(x = 'hidden_units', y = 'precision_score', data = df_main,  ax = ax, markers = '.', size = 'n_layers')
        ax.set_xlabel('# Hidden Units')
        ax.set_ylabel('Precision Score')
        ax.set_title('scAGN: {}'.format(data))
        plt.tight_layout()
        plt.savefig('AGN_'+data + '_precision_scoreVsHiddenUnits.png', bbox_inches='tight')
        plt.savefig('AGN_'+data + '_precision_scoreVsHiddenUnits.pdf', bbox_inches='tight')
        
        sea.set_context('paper')
        sea.set_style("dark")
        fig, ax = plt.subplots(figsize=(3, 3))
        sea.lineplot(x = 'hidden_units', y = 'recall_score', data = df_main,  ax = ax, markers = True)
        sea.scatterplot(x = 'hidden_units', y = 'recall_score', data = df_main,  ax = ax, markers = '.', size = 'n_layers')
        ax.set_xlabel('# Hidden Units')
        ax.set_ylabel('Recall Score')
        ax.set_title('scAGN: {}'.format(data))
        plt.tight_layout()
        plt.savefig('AGN_'+data + '_recall_scoreVsHiddenUnits.png', bbox_inches='tight')
        plt.savefig('AGN_'+data + '_recall_scoreVsHiddenUnits.pdf', bbox_inches='tight')
        
        sea.set_context('paper')
        sea.set_style("dark")
        fig, ax = plt.subplots(figsize=(3, 3))
        sea.lineplot(x = 'hidden_units', y = 'f1_score', data = df_main,  ax = ax, markers = True)
        sea.scatterplot(x = 'hidden_units', y = 'f1_score', data = df_main,  ax = ax, markers = '.', size = 'n_layers')
        ax.set_xlabel('# Hidden Units')
        ax.set_ylabel('F1 Score')
        ax.set_title('scAGN: {}'.format(data))
        plt.tight_layout()
        plt.savefig('AGN_'+data + '_f1_scoreVsHiddenUnits.png', bbox_inches='tight')
        plt.savefig('AGN_'+data + '_f1_scoreVsHiddenUnits.pdf', bbox_inches='tight')
            
        print("****")
        id_file = df_main2['predication_accuracy'].idxmax()

        runtime = pd.read_csv('{}_AGN_elapsed_time.txt'.format(data), index_col = 0, header = None)
        time_taken = runtime.iloc[id_file][2]
        runtime_df.at[data, 'scAGN'] = time_taken
         

        print(id_file)
        true_labels_file = AGN_prediction_files[id_file][0:-21] + 'true_labels.npy'
        predicted_labels_file = AGN_prediction_files[id_file][0:-21] + 'predicted_labels.npy'
        true_labels = np.load(true_labels_file)
        predicted_labels = np.load(predicted_labels_file)
        cf_matrix = confusion_matrix(true_labels, predicted_labels)
        normed_c = (cf_matrix.T / cf_matrix.astype(float).sum(axis=1)).T
        label_mapping = pd.read_csv(folder + '/label_mapping.csv')
        fig, ax = plt.subplots(1)
        g = sea.heatmap(normed_c, annot=False,  cmap='Blues')
        g.invert_yaxis()
        ax.set_xticks(np.arange(0, len(label_mapping[label_mapping.columns[1]].values.tolist()) ))
        ax.set_yticks(np.arange(0, len(label_mapping[label_mapping.columns[1]].values.tolist()) ))
        ax.xaxis.set_ticklabels(label_mapping[label_mapping.columns[1]].values.tolist(), rotation = 90)
        ax.yaxis.set_ticklabels(label_mapping[label_mapping.columns[1]].values.tolist(), rotation = 0)
        ax.set_title(data)
        plt.savefig(folder + '/Confusion_matrix.png', bbox_inches='tight')
        plt.savefig(folder + '/Confusion_matrix.pdf', bbox_inches='tight')
        
        df_m = pd.DataFrame(columns=['dataset', 'prediction_accuracy', 'precision_score', 'recall_score', 'f1_score'])
        df_m.loc[0] = [data, np.max(df_main['predication_accuracy']), np.max(df_main['precision_score']), np.max(df_main['recall_score']), np.max(df_main['f1_score'])]
        print(df_m)
        df_metric = pd.concat([df_metric, df_m])


df_metric = df_metric.set_index('dataset')
df_metric.index.name = None
#print(df_metric)
# now we need to calculate metrics for other methods.

folders_for_baseline = "../RScript/scData"
folders_for_prediction_acc = "../RScript"
df_predicted_seu = pd.read_csv(folders_for_prediction_acc + '/seurat_metric.txt', index_col = 0, header = None)
df_predicted_seu.columns = ["seurat_prediction"]
#print(df_predicted_seu)

df_predicted_singleR = pd.read_csv(folders_for_prediction_acc + '/SingleR_metric.txt', index_col = 0, header = None)
df_predicted_singleR.columns = ["singleR_prediction"]
#print(df_predicted_singleR)

df_predicted_scmap_cluster = pd.read_csv(folders_for_prediction_acc + '/scMAP-cluster_metric.txt', index_col = 0, header = None)
df_predicted_scmap_cluster.columns = ["scmap-cluster_prediction"]
#print(df_predicted_scmap_cluster)

df_predicted_scmap_cell = pd.read_csv(folders_for_prediction_acc + '/scMAP-cell_metric.txt', index_col = 0, header = None)
df_predicted_scmap_cell.columns = ["scmap-cell_prediction"]
#print(df_predicted_scmap_cell)

df_predicted_chetah = pd.read_csv(folders_for_prediction_acc + '/Chetah_metric.txt', index_col = 0, header = None)
df_predicted_chetah.columns = ["chetah_prediction"]
#print(df_predicted_chetah)

df_pred_accuracy = pd.concat([df_metric, df_predicted_seu, df_predicted_singleR, df_predicted_scmap_cluster, df_predicted_scmap_cell, df_predicted_chetah], axis = 1)
df_pred_accuracy.rename(columns={'prediction_accuracy':'scAGN_prediction', 'precision_score': 'scAGN_precision_score', 'recall_score': 'scAGN_recall_score', 'precision_score': 'scAGN_precision_score', 'f1_score': 'scAGN_f1_score'}, inplace=True)
#df_pred_accuracy.drop(columns=['precision_score', 'recall_score','f1_score'], inplace=True)

for data in datasets_to_include:
    print(data)
    print("----------")
    folder = folders_for_baseline + "/" + data
    #SeuratV3
    if os.path.exists(folder + '/true_labels_SeuratV3.csv'): 
        df_true_labels_seu = pd.read_csv(folder + '/true_labels_SeuratV3.csv', index_col = 0)
        df_pred_labels_seu = pd.read_csv(folder + '/predicted_labels_SeuratV3.csv', index_col = 0)
        print(df_true_labels_seu.columns)
        print(df_pred_labels_seu.columns)
        precision_score = sklearn.metrics.precision_score(df_true_labels_seu['query.object$type'].values, df_pred_labels_seu['predicted.id'].values, average='weighted')
        recall_score = sklearn.metrics.recall_score(df_true_labels_seu['query.object$type'].values, df_pred_labels_seu['predicted.id'].values, average='weighted')
        f1_score = sklearn.metrics.f1_score(df_true_labels_seu['query.object$type'].values, df_pred_labels_seu['predicted.id'].values, average='weighted')
        df_pred_accuracy.at[data,'Seurat_precision_score'] = precision_score
        df_pred_accuracy.at[data,'Seurat_recall_score'] = recall_score
        df_pred_accuracy.at[data,'Seurat_f1_score'] = f1_score
 
    #SingleR
    if os.path.exists(folder + '/true_labels_SingleR.csv'): 
        df_true_labels_singleR = pd.read_csv(folder + '/true_labels_SingleR.csv', index_col = 0)
        df_pred_labels_singleR = pd.read_csv(folder + '/predicted_labels_SingleR.csv', index_col = 0)
        print(df_pred_labels_singleR.columns)
        print(df_pred_labels_singleR.columns)
        precision_score = sklearn.metrics.precision_score(df_true_labels_singleR['type'].values, df_pred_labels_singleR['singler$labels'].values, average='weighted')
        recall_score = sklearn.metrics.recall_score(df_true_labels_singleR['type'].values, df_pred_labels_singleR['singler$labels'].values, average='weighted')
        f1_score = sklearn.metrics.f1_score(df_true_labels_singleR['type'].values, df_pred_labels_singleR['singler$labels'].values, average='weighted')
        df_pred_accuracy.at[data,'SingleR_precision_score'] = precision_score
        df_pred_accuracy.at[data,'SingleR_recall_score'] = recall_score
        df_pred_accuracy.at[data,'SingleR_f1_score'] = f1_score
 
    #scmap-Cluster
    if os.path.exists(folder + '/true_labels_scMAP-cluster.csv'): 
        df_true_labels_scmap_cluster = pd.read_csv(folder + '/true_labels_scMAP-cluster.csv', index_col = 0)
        df_pred_labels_scmap_cluster = pd.read_csv(folder + '/predicted_labels_scMAP-cluster.csv', index_col = 0)
        df_labs = df_true_labels_scmap_cluster
        df_labs['predicted_type'] = df_pred_labels_scmap_cluster[df_pred_labels_scmap_cluster.columns[0]].values
        df_labs = df_labs[df_labs['predicted_type'] != 'unassigned']
        df_labs.to_csv('testcheck.csv')
        df_labs = pd.read_csv('testcheck.csv',index_col = 0)
        precision_score = sklearn.metrics.precision_score(df_labs['type'].values, df_labs['predicted_type'].values, average='weighted')
        recall_score = sklearn.metrics.recall_score(df_labs['type'].values, df_labs['predicted_type'].values, average='weighted')
        f1_score = sklearn.metrics.f1_score(df_labs['type'].values, df_labs['predicted_type'].values, average='weighted')

        df_pred_accuracy.at[data,'Scmap-cluster_precision_score'] = precision_score
        df_pred_accuracy.at[data,'Scmap-cluster_recall_score'] = recall_score
        df_pred_accuracy.at[data,'Scmap-cluster_f1_score'] = f1_score
 
    #scmap-Cell
    if os.path.exists(folder + '/true_labels_scMAP-cell.csv'): 
        df_true_labels_scmap_cell = pd.read_csv(folder + '/true_labels_scMAP-cell.csv', index_col = 0)
        df_pred_labels_scmap_cell = pd.read_csv(folder + '/predicted_labels_scMAP-cell.csv', index_col = 0)
        # df_true_labels_scmap_cell.loc[(df_true_labels_scmap_cell['type'] == 'unassigned'),'type'] = -1
        # df_pred_labels_scmap_cell.loc[(df_pred_labels_scmap_cell[df_pred_labels_scmap_cell.columns[0]] == 'unassigned'),df_pred_labels_scmap_cell.columns[0]] = -1
        df_labs = df_true_labels_scmap_cell
        df_labs['predicted_type'] = df_pred_labels_scmap_cell[df_pred_labels_scmap_cell.columns[0]].values
        df_labs = df_labs[df_labs['predicted_type'] != 'unassigned']
        df_labs.to_csv('testcheck.csv')
        df_labs = pd.read_csv('testcheck.csv',index_col = 0)

        precision_score = sklearn.metrics.precision_score(df_labs['type'].values, df_labs['predicted_type'].values, average='weighted')
        recall_score = sklearn.metrics.recall_score(df_labs['type'].values, df_labs['predicted_type'].values, average='weighted')
        f1_score = sklearn.metrics.f1_score(df_labs['type'].values, df_labs['predicted_type'].values, average='weighted')
        df_pred_accuracy.at[data,'Scmap-cell_precision_score'] = precision_score
        df_pred_accuracy.at[data,'Scmap-cell_recall_score'] = recall_score
        df_pred_accuracy.at[data,'Scmap-cell_f1_score'] = f1_score
 
    #Chetah
    if os.path.exists(folder + '/true_labels_Chetah.csv'): 
        df_true_labels_chetah = pd.read_csv(folder + '/true_labels_Chetah.csv', index_col = 0, dtype={'type': object})
        df_pred_labels_chetah = pd.read_csv(folder + '/predicted_labels_Chetah.csv', index_col = 0, dtype={'test$celltype_CHETAH': object})
        print(df_true_labels_chetah.columns)
        print(df_pred_labels_chetah.columns)
        print(df_true_labels_chetah['type'].values)
        print(df_pred_labels_chetah['test$celltype_CHETAH'].values)
        precision_score = sklearn.metrics.precision_score(df_true_labels_chetah['type'].values, df_pred_labels_chetah['test$celltype_CHETAH'].values, average='weighted')
        recall_score = sklearn.metrics.recall_score(df_true_labels_chetah['type'].values, df_pred_labels_chetah['test$celltype_CHETAH'].values, average='weighted')
        f1_score = sklearn.metrics.f1_score(df_true_labels_chetah['type'].values, df_pred_labels_chetah['test$celltype_CHETAH'].values, average='weighted')
        df_pred_accuracy.at[data,'Chetah_precision_score'] = precision_score
        df_pred_accuracy.at[data,'Chetah_recall_score'] = recall_score
        df_pred_accuracy.at[data,'Chetah_f1_score'] = f1_score
 

print(df_pred_accuracy)


df_pred_accuracy.to_csv('Predication_Comparison.csv')


dftime_seurat = pd.read_csv('../RScript/seurat_time_taken.txt', index_col = 0, header = None)
dftime_seurat.columns = ['Seurat']

dftime_singleR = pd.read_csv('../RScript/SingleR_time_taken.txt', index_col = 0, header = None)
dftime_singleR.columns = ['SingleR']

dftime_scmap_cluster = pd.read_csv('../RScript/scMP-cluster_time_taken.txt', index_col = 0, header = None)
dftime_scmap_cluster.columns = ['Scmap-cluster']

dftime_scmap_cell = pd.read_csv('../RScript/scMAP-cell_time_taken.txt', index_col = 0, header = None)
dftime_scmap_cell.columns = ['Scmap-cell']

dftime_chetah = pd.read_csv('../RScript/Chetah_time_taken.txt', index_col = 0, header = None)
dftime_chetah.columns = ['Chetah']

dftime = pd.concat([runtime_df, dftime_seurat, dftime_singleR, dftime_scmap_cluster, dftime_scmap_cell, dftime_chetah], axis = 1)
print(dftime)
dftime.to_csv('alltime.csv')
df11 = dftime[['TM', 'Tasic', 'Mouse_retina', 'PBMC_68k', 'GSE108989', 'GSE118389', 'GSE72056', 'GSE98638', 'GSE99254', 'blish_pbmc_covid19', 'shalek_nasal_epithelia_covid19', 'villani_mgh_covid19']]

fig, ax = plt.subplots(figsize=(7, 4))
df11.plot(kind="bar", ax = ax)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.ylabel('Runtime  [s]')
plt.yticks(fontsize=16) # for yticks
ax.set_yscale('log')
plt.title('Runtime of classifiers', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('runtime.png', bbox_inches='tight')
plt.savefig('runtime.pdf', bbox_inches='tight')



sea.set_context('paper')
sea.set_style("dark")
fig, ax = plt.subplots(figsize=(14, 14))
df_pred_accuracy['Dataset'] = df_pred_accuracy.index.tolist()
#sea.lineplot(x = 'hidden_units', y = 'f1_score', data = df_main,  ax = ax, markers = True)
sea.scatterplot(x = 'seurat_prediction', y = 'scAGN_prediction', data = df_pred_accuracy,  ax = ax, markers = '.', hue = 'Dataset')
ax.set_xlabel('Seurat')
ax.set_ylabel('scAGN')
ax.set_title('Prediction Accuracy Comparison')
ax.axis('equal')
x = np.linspace(*ax.get_xlim())
ax.plot(x, x)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.savefig('AGNvsSeurat.png', bbox_inches='tight')
plt.savefig('AGNvsSeurat.pdf', bbox_inches='tight')

df_pred_accuracy.drop(columns=['Dataset'], inplace=True)
data1 = ['TM', 'Tasic', 'Mouse_retina', 'PBMC_68k']
df1 = df_pred_accuracy.filter(items = data1, axis=0)
print(df1)
sea.set_context('paper')
sea.set_style("dark")
fig, ax = plt.subplots(figsize=(7, 4))
df1[['scAGN_prediction', 'seurat_prediction', 'singleR_prediction', 'scmap-cluster_prediction', 'scmap-cell_prediction', 'chetah_prediction']].plot(kind="bar", ax = ax)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('Prediction Accuracy: Healthy Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_Prediction_Accuracy_Healthy.png', bbox_inches='tight')
plt.savefig('Comparison_Prediction_Accuracy_Healthy.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
df1[['scAGN_precision_score', 'Seurat_precision_score', 'SingleR_precision_score', 'Scmap-cluster_precision_score', 'Scmap-cell_precision_score', 'Chetah_precision_score']].plot(kind="bar", ax = ax)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('Precision Score: Healthy Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_Precision_Score_Healthy.png', bbox_inches='tight')
plt.savefig('Comparison_Precision_Score_Healthy.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
df1[['scAGN_recall_score', 'Seurat_recall_score', 'SingleR_recall_score', 'Scmap-cluster_recall_score', 'Scmap-cell_recall_score', 'Chetah_recall_score']].plot(kind="bar", ax = ax)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('Recall Score: Healthy Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_Recall_Score_Healthy.png', bbox_inches='tight')
plt.savefig('Comparison_Recall_Score_Healthy.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
df1[['scAGN_f1_score', 'Seurat_f1_score', 'SingleR_f1_score', 'Scmap-cluster_f1_score', 'Scmap-cell_f1_score', 'Chetah_f1_score']].plot(kind="bar", ax = ax)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('F1 Score: Healthy Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_F1_Score_Healthy.png', bbox_inches='tight')
plt.savefig('Comparison_F1_Score_Healthy.pdf', bbox_inches='tight')

data2 = ['GSE108989', 'GSE118389', 'GSE72056', 'GSE98638', 'GSE99254']
df2 = df_pred_accuracy.filter(items = data2, axis=0)
print(df2)
sea.set_context('paper')
sea.set_style("dark")
fig, ax = plt.subplots(figsize=(7, 4))
df2[['scAGN_prediction', 'seurat_prediction', 'singleR_prediction', 'scmap-cluster_prediction', 'scmap-cell_prediction', 'chetah_prediction']].plot(kind="bar", ax = ax)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('Prediction Accuracy: Cancer Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_Prediction_Accuracy_Cancer.png', bbox_inches='tight')
plt.savefig('Comparison_Prediction_Accuracy_Cancer.pdf', bbox_inches='tight')


fig, ax = plt.subplots(figsize=(7, 4))
df2[['scAGN_precision_score', 'Seurat_precision_score', 'SingleR_precision_score', 'Scmap-cluster_precision_score', 'Scmap-cell_precision_score', 'Chetah_precision_score']].plot(kind="bar", ax = ax)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('Precision Score: Cancer Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_Precision_Score_Cancer.png', bbox_inches='tight')
plt.savefig('Comparison_Precision_Score_Cancer.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
df2[['scAGN_recall_score', 'Seurat_recall_score', 'SingleR_recall_score', 'Scmap-cluster_recall_score', 'Scmap-cell_recall_score', 'Chetah_recall_score']].plot(kind="bar", ax = ax)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('Recall Score: Cancer Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_Recall_Score_Cancer.png', bbox_inches='tight')
plt.savefig('Comparison_Recall_Score_Cancer.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
df2[['scAGN_f1_score', 'Seurat_f1_score', 'SingleR_f1_score', 'Scmap-cluster_f1_score', 'Scmap-cell_f1_score', 'Chetah_f1_score']].plot(kind="bar", ax = ax)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('F1 Score: Cancer Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_F1_Score_Cancer.png', bbox_inches='tight')
plt.savefig('Comparison_F1_Score_Cancer.pdf', bbox_inches='tight')

data3 = ['blish_pbmc_covid19', 'shalek_nasal_epithelia_covid19', 'villani_mgh_covid19']
df3 = df_pred_accuracy.filter(items = data3, axis=0)
print(df3)
sea.set_context('paper')
sea.set_style("dark")
fig, ax = plt.subplots(figsize=(7, 4))
df3[['scAGN_prediction', 'seurat_prediction', 'singleR_prediction', 'scmap-cluster_prediction', 'scmap-cell_prediction', 'chetah_prediction']].plot(kind="bar", ax = ax)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('Prediction Accuracy: COVID19 Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_Prediction_Accuracy_COVID19.png', bbox_inches='tight')
plt.savefig('Comparison_Prediction_Accuracy_COVID19.pdf', bbox_inches='tight')


fig, ax = plt.subplots(figsize=(7, 4))
df3[['scAGN_precision_score', 'Seurat_precision_score', 'SingleR_precision_score', 'Scmap-cluster_precision_score', 'Scmap-cell_precision_score', 'Chetah_precision_score']].plot(kind="bar", ax = ax)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('Precision Score: COVID19 Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_Precision_Score_COVID19.png', bbox_inches='tight')
plt.savefig('Comparison_Precision_Score_COVID19.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
df3[['scAGN_recall_score', 'Seurat_recall_score', 'SingleR_recall_score', 'Scmap-cluster_recall_score', 'Scmap-cell_recall_score', 'Chetah_recall_score']].plot(kind="bar", ax = ax)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('Recall Score: COVID19 Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_Recall_Score_COVID19.png', bbox_inches='tight')
plt.savefig('Comparison_Recall_Score_COVID19.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(7, 4))
df3[['scAGN_f1_score', 'Seurat_f1_score', 'SingleR_f1_score', 'Scmap-cluster_f1_score', 'Scmap-cell_f1_score', 'Chetah_f1_score']].plot(kind="bar", ax = ax)
plt.tight_layout()
plt.xticks(fontsize=16)  # for xticks
plt.yticks(fontsize=16) # for yticks
plt.title('F1 Score: COVID19 Cells', fontsize = 18)
plt.legend(['scAGN', 'Seurat', 'SingleR', 'scmap-cluster', 'scmap-cell', 'Chetah'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 14)
plt.savefig('Comparison_F1_Score_COVID19.png', bbox_inches='tight')
plt.savefig('Comparison_F1_Score_COVID19.pdf', bbox_inches='tight')
