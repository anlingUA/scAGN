#!/usr/bin/env python

# Initial Date: October 2021
# Author: Rahul Bhadani
# Copyright (c)  Rahul Bhadani
# All rights reserved.

__author__ = 'Rahul Bhadani'
__email__  = 'rahulbhadani@email.arizona.edu'

import torch
from torch.functional import norm
from functions import configure_logworker
import pandas as pd
import numpy as np
import torch_geometric
from torch_geometric.data import Data
from sklearn.preprocessing import StandardScaler
from prepare_data import generate_dataset
from torch.optim.lr_scheduler import StepLR
import sklearn.metrics
import time

from AGNNet import AGNNet

import sys, getopt
from functions import dict_reverser
import torch.optim.lr_scheduler as Scheduler
from EarlyStopping import EarlyStopping
import seaborn as sns
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt

def main(argv):
    # What kind of different parameters to the script?
    # 1. number of hidden inputs
    # 2. number of layers
    # 3. dropout value
    # 4. learning rate
    # 5. model weight decay
    nhidden = 16
    nlayers = 5
    dropout = 0.1
    learning_rate = 0.0001
    weight_decay = 5e-4
    loss_threshold = 0.1
    normalize_data = False
    inductive_training = False
    graphmethod = "umap"
    dataset = 'GSE84133'
    console_output = False
    gamma = 0.8
    step_size = 10
    patience = 100

    try:
        opts, args = getopt.getopt(argv,"Ihnpi:d:r:w:t:g:l:D:G:s:P:",["nhidden=","dropout=","learning_rate=", "weight_decay=" , "loss_threshold=", "graph=", "nlayers=", "dataset=","gamma=", "step_size=", 'patience='])
        if len(opts) == 0:
            print('Check options by typing:\n{} -h'.format(__file__))
            sys.exit()

    except getopt.GetoptError:
        print('Check options by typing:\n{} -h'.format(__file__))
        sys.exit(2)

    print("OPTS: {}".format(opts))
    for opt, arg in opts:
        if(opt == '-h'):
            print('\n{} [OPTIONS]'.format(__file__))
            print('\t -h, --help\t\t Get help')
            print('\t -p, --print\t\t\t Print output')
            print('\t -n, --normalize\t Normalize the data')
            print('\t -i, --nhidden\t\t Number of hidden units')
            print('\t -l, --nlayers\t\t Number of hidden layers')
            print('\t -d, --dropout\t\t Dropout value')
            print('\t -r, --learning_rate\t Learning rate')
            print('\t -w, --weight_decay\t Weight Decay rate')
            print('\t -t, --loss_threshold\t Loss Threshold')
            print('\t -D, --dataset\t\t Dataset to process')
            print('\t -g, --graph\t\t Graph construction Method. Available options: "knn", "radius", "umap", "pearson", "spearman", "kendall"')
            print('\t -G, --gamma\t\t Multiplication factor to learning rate to reduce learning with epochs')
            print('\t -s, --step_size\t\t After every s stepsize, learning rate reduces')
            print('\t -P, --patience\t\t Patience, i.e. number of steps to wait for early stopping of training')
            sys.exit()
        elif(opt == '-I'):
            inductive_training = True
        elif(opt in ("-p", "--print")):
            console_output = True
        elif(opt == '-n'):
            normalize_data = True
        elif(opt in ("-i", "--nhidden")):
            nhidden = eval(arg)
        elif(opt in ("-d", "--dropout")):
            dropout = eval(arg)
        elif(opt in ("-l", "--nlayers")):
            nlayers = eval(arg)
        elif(opt in ("-r", "--learning_rate")):
            learning_rate = eval(arg)
        elif(opt in ("-w", "--weight_decay")):
            weight_decay = eval(arg)
        elif(opt in ("-D", "--dataset")):
            dataset = arg
        elif(opt in ("-t", "--loss_threshold")):
            loss_threshold = eval(arg)
        elif(opt in ("-g", "--graph")):
            graphmethod = arg
        elif(opt in ("-G", "--gamma")):
            gamma = eval(arg)
        elif(opt in ("-s", "--step_size")):
            step_size = eval(arg)
        elif(opt in ("-P", "--patience")):
            patience = eval(arg)

    print("Step Size of Learning Rate Reduction: {}".format(step_size))
    print("Gamma of Learning Rate Reduction: {}".format(gamma))


    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    if inductive_training:
        print("Doing Inductive Training with graph method: {}".format(graphmethod))

        folder = "../RScript/{}_ref_{}/".format(dataset, graphmethod)
        edge_file = 'intra_graph.csv'
        label1_file = 'label1.csv'
        edge_df = pd.read_csv(folder + edge_file)
        edge_df = edge_df.drop(columns=['Unnamed: 0', 'index'])
        labels =  pd.read_csv(folder + label1_file)

        feature_df = pd.read_csv(folder + 'Data1.csv')
        feature_df = feature_df.set_index(keys='Unnamed: 0')
        labels = labels.set_index(keys='Unnamed: 0')
        feature_df.index = feature_df.index.rename('cellid')
        labels.index = labels.index.rename('cellid')
        Y = labels['type'].tolist()

        if (normalize_data):
            print("We will normalize feature matrix before feeding to AGN Net.")
            feature_df = pd.DataFrame(StandardScaler().fit_transform(feature_df), columns=feature_df.columns, index=feature_df.index)
        else:
            print("We will not normalize feature matrix of reference before feeding to AGN Net.")
        # Construct graph for training
        X = feature_df.values.tolist()
        y = torch.tensor(Y, dtype=torch.long)
        edge_index_tensor = torch.tensor(edge_df.values.tolist(), dtype=torch.long)
        x = torch.tensor(X, dtype=torch.float)
        reference_data = Data(x=x, edge_index=edge_index_tensor.t().contiguous(), y=y)


        ########################################################################################
        folder = "../RScript/{}_query_{}/".format(dataset, graphmethod)
        edge_file = 'intra_graph.csv'
        label1_file = 'label1.csv'
        edge_df = pd.read_csv(folder + edge_file)
        edge_df = edge_df.drop(columns=['Unnamed: 0', 'index'])
        labels =  pd.read_csv(folder + label1_file)


        feature_df = pd.read_csv(folder + 'Data1.csv')
        feature_df = feature_df.set_index(keys='Unnamed: 0')
        labels = labels.set_index(keys='Unnamed: 0')
        feature_df.index = feature_df.index.rename('cellid')
        labels.index = labels.index.rename('cellid')
        Y = labels['type'].tolist()


        if (normalize_data):
            print("We will normalize feature matrix before feeding to AGN Net.")
            feature_df = pd.DataFrame(StandardScaler().fit_transform(feature_df), columns=feature_df.columns, index=feature_df.index)
        else:
            print("We will not normalize feature matrix of query data before feeding to AGN Net.")
        # Construct graph for training
        X = feature_df.values.tolist()
        y = torch.tensor(Y, dtype=torch.long)
        edge_index_tensor = torch.tensor(edge_df.values.tolist(), dtype=torch.long)
        x = torch.tensor(X, dtype=torch.float)
        query_data = Data(x=x, edge_index=edge_index_tensor.t().contiguous(), y=y)

        query_data = query_data.to(device)
        reference_data = reference_data.to(device)

        print("Query data dimensions: {}".format(query_data))
        print("Reference data dimensions: {}".format(reference_data))

        model = AGNNet(reference_data, nhidden=nhidden, dropout=dropout, description='Training AGNNet on {} with {} Graph'.format(dataset, graphmethod), lr =learning_rate, \
            loss_threshold = loss_threshold, weight_decay = weight_decay, query_data = query_data, nlayers = nlayers).to(device)

        _LOGGER = configure_logworker(logfile =  model.logfile)
        _LOGGER.info(model.train())

        _LOGGER.info(model.config)
        _LOGGER.info("Model learning rate is {}".format(model.lr))

        if (normalize_data):
            _LOGGER.info("We will normalize feature matrix before feeding to AGN Net.")
        else:
            _LOGGER.info("We will not normalize feature matrix before feeding to AGN Net.")

        _LOGGER.info("We will do inductive learning.")

        model.optimizer = torch.optim.Adam(model.parameters(), lr=model.lr, weight_decay=model.weight_decay)
        # model.scheduler  = Scheduler.ReduceLROnPlateau(model.optimizer, 'min', patience = 200, threshold=0.0000001, factor = 0.5)
        model.inductive_loss_training()


    else:
        print("Transductive Training")

        DataDir = '../data/{}/'.format(dataset)
        data = generate_dataset(DataDir=DataDir)

        data = data.to(device)

        description ='AGNNet_graph_i_{}_l_{}_d_{}_t_{}_r_{}_w_{}_s_{}_G_{}_D_{}_P_{}_'.format(nhidden, nlayers, dropout, loss_threshold, learning_rate, weight_decay, step_size, gamma,  dataset, patience)
        description = description.replace(" ", "")

        model = AGNNet(data, nhidden=nhidden, dropout=dropout, description=description, nlayers = nlayers, \
            lr =learning_rate, weight_decay = weight_decay, loss_threshold = loss_threshold).to(device)

        _LOGGER = configure_logworker(logfile =  model.logfile)
        _LOGGER.info(model.train())

        _LOGGER.info(model.config)
        _LOGGER.info("Model learning rate is {}".format(model.lr))

        if (normalize_data):
            _LOGGER.info("We will normalize feature matrix before feeding to AGN Net.")
        else:
            _LOGGER.info("We will not normalize feature matrix before feeding to AGN Net.")


        model.optimizer = torch.optim.Adam(model.parameters(), lr=model.lr, weight_decay=model.weight_decay)
        model.earlystopping = EarlyStopping(patience=patience, verbose=True, trace_func=_LOGGER.info)

        model.scheduler = StepLR(model.optimizer, step_size=step_size, gamma=gamma, verbose=False)
        start = time.time()
        model.loss_training(console_output = console_output)
        end = time.time()

        sourceFile = open('./{}_AGN_elapsed_time.txt'.format(dataset), 'a')
        print('{},{},{}'.format(dataset, model.description, end - start), file = sourceFile)
        sourceFile.close()

        pred = model.evaluate(data, mask='pred_mask')
        prediction_probabilities = np.exp(pred.cpu().detach().numpy())
        pred_numpy  = pred.cpu().detach().numpy()
        prediction_ranking = pd.DataFrame(pred_numpy).rank(axis = 1, ascending = False)
        predicted = pred.max(1)[1].cpu().detach().numpy()
        true_labels = data.y[data['pred_mask']].cpu().detach().numpy()

        predication_accuracy = pred.max(1)[1].eq(data.y[data['pred_mask']]).sum().item() / data['pred_mask'].sum().item()
        precision_score = sklearn.metrics.precision_score(true_labels, predicted, average='weighted')
        recall_score = sklearn.metrics.recall_score(true_labels, predicted, average='weighted')
        f1_score= sklearn.metrics.f1_score(true_labels, predicted, average='weighted')

        label_mapping = pd.read_csv(DataDir + '/label_mapping.csv')
        print(label_mapping[label_mapping.columns[1]].values.tolist())

        cl_report = sklearn.metrics.classification_report(true_labels, predicted, target_names=label_mapping['label_names'].astype('str').tolist())
        sourceFile = open(DataDir + model.description + 'classification_report.csv', 'w')
        print(cl_report, file = sourceFile)
        sourceFile.close()

        _LOGGER.info("Prediction Accuracy is {}".format(predication_accuracy))

        cf_matrix = confusion_matrix(true_labels, predicted)
        normed_c = (cf_matrix.T / cf_matrix.astype(float).sum(axis=1)).T

        

        # save some data
        df_metric = pd.DataFrame([predication_accuracy, precision_score, recall_score, f1_score], columns=['metric'], index = ['predication_accuracy', 'precision_score', 'recall_score', 'f1_score'])

        df_metric.to_csv(DataDir + model.description + 'prediction_metric.csv')
        np.save(DataDir + model.description + 'prediction_probabilities.npy', prediction_probabilities)
        np.save(DataDir + model.description + 'predicted_labels.npy', predicted)
        np.save(DataDir + model.description + 'true_labels.npy', true_labels)

        sns.set_context('talk')
        fig, ax = plt.subplots(1)
        for i in range(0, len(label_mapping)):
            fpr, tpr, _ = sklearn.metrics.roc_curve(np.int32(data.y[data['pred_mask']].cpu().detach().numpy() == i), prediction_probabilities[:,i])
            plt.plot(fpr, tpr, label = label_mapping['label_names'].loc[i])
        plt.legend(bbox_to_anchor=(-0.1,1))
        plt.savefig(DataDir + model.description + model.dt + '_ROC_Curve.png', bbox_inches='tight')
        plt.savefig(DataDir + model.description + model.dt + '_ROC_Curve.pdf', bbox_inches='tight')

        fig, ax = plt.subplots(1)
        g = sns.heatmap(normed_c, annot=False,  cmap='Blues')
        g.invert_yaxis()
        ax.set_xticks(np.arange(0, len(label_mapping[label_mapping.columns[1]].values.tolist()) ))
        ax.set_yticks(np.arange(0, len(label_mapping[label_mapping.columns[1]].values.tolist()) ))
        ax.xaxis.set_ticklabels(label_mapping[label_mapping.columns[1]].values.tolist(), rotation = 90)
        ax.yaxis.set_ticklabels(label_mapping[label_mapping.columns[1]].values.tolist(), rotation = 0)
        ax.set_title(dataset)
        plt.savefig(DataDir + model.description + model.dt + '_confusion_matrix.png', bbox_inches='tight')
        plt.savefig(DataDir + model.description + model.dt + '_confusion_matrix.pdf', bbox_inches='tight')
        

if __name__ == "__main__":
   main(sys.argv[1:])
