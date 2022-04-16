#!/usr/bin/env python

# Initial Date: October 2021
# Author: Rahul Bhadani
# Copyright (c)  Rahul Bhadani
# All rights reserved.

__author__ = 'Rahul Bhadani'
__email__  = 'rahulbhadani@email.arizona.edu'

from functions import configure_logworker
import sys, time, subprocess, os, datetime, glob, signal
import ntpath
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from anndata import AnnData
import scanpy as sc
import torch
from torch_geometric.data import Data
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import AGNNConv
import numpy as np
import pandas as pd
import torch_geometric
import csv
import torch.optim.lr_scheduler as Scheduler
import yaml
import socket
sock_name = socket.gethostname()

from EarlyStopping import EarlyStopping

class AGNNet(nn.Module):
    """
    `AGNNet` A neural network with AGN Conv Layers.
    
    Parameters
    ------------
    data: `torch_geometric.datasets` type dataset
     
    
    """
    def __init__(self, data, nhidden = 16, dropout = 0.1, description = 'AGNNet', **kwargs):
        super(AGNNet, self).__init__()

        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        dt_object = datetime.datetime.fromtimestamp(time.time())
        dt = dt_object.strftime('%Y_%m_%d__%H_%M_%S_%f')
        self.dt = dt
        self.device = device
        self.data = data.to(device)
        self.query_data = kwargs.get("query_data", None)
        if self.query_data is not None:
            self.query_data = self.query_data.to(device)
        
        self.description = description
        if 'arizona' in sock_name:
            self.dir = kwargs.get("dir","/xdisk/anling/rahulbhadani/gene-cn-output/")
        else:
            self.dir = kwargs.get("dir", "../../gene-cn-output/")

        self.logfile = self.dir +  self.description.strip().replace(" ","_") + "_AGNNet_log_" + self.dt + ".log"
        self.modelckdir = self.dir + self.description.strip().replace(" ","_") + "_AGNNet_modelck_" + self.dt
        self.checkpoint_tracker = 1 #checkpoint tracker to save model to alternate file in every epoch to avoid incomplete writing upon termination     
        self.configfile = self.dir +  self.description.strip().replace(" ","_") + "_AGNNet_config_" + self.dt + ".yaml"
        
        self._LOGGER = configure_logworker(logfile = self.logfile)
        self.num_classes = len(np.unique(self.data.y.cpu()))
        self._LOGGER.info("Number of classes is {}".format(self.num_classes))
        
        self._LOGGER.info("Reference Data feature: {}, dim: {}".format(self.data.x, self.data.x.shape))
        self._LOGGER.info("Reference Data edge matrix: {}, dim: {}".format(self.data.edge_index, self.data.edge_index.shape))
        self._LOGGER.info("Reference Data labels: {}, dim: {}".format(self.data.y, self.data.y.shape))

        if self.query_data is not None:
            self._LOGGER.info("Query Data feature: {}, dim: {}".format(self.query_data.x, self.query_data.x.shape))
            self._LOGGER.info("Query Data edge matrix: {}, dim: {}".format(self.query_data.edge_index, self.query_data.edge_index.shape))
            self._LOGGER.info("Query Data labels: {}, dim: {}".format(self.query_data.y, self.query_data.y.shape))
            
        self.nhidden = nhidden
        self.nlayers = kwargs.get("nlayers", 2)
        self.convs = nn.ModuleList()

        self.linear1 = nn.Linear(self.data.num_features, nhidden)
        self.conv2 = AGNNConv(requires_grad=False)
        self.conv3 = AGNNConv(requires_grad=True)
        self.convs.append(self.conv2)
        self.convs.append(self.conv3)
        if self.nlayers > 4:
            for i in range(0, self.nlayers-3):
                self.convs.append(AGNNConv(requires_grad=True))
        self.linear4 = nn.Linear(nhidden, self.num_classes) 

        self.dropout = dropout
        self.lr = kwargs.get("lr", 0.00001)
        self.weight_decay = kwargs.get("weight_decay", 5e-4)
        self.optimizer = None
        self.scheduler = None
        self.earlystopping = None
        self.loss_threshold = kwargs.get("loss_threshold", 0.01)
        self.config = {}
        self.config['nhidden'] = self.nhidden
        self.config['nlayers'] = self.nlayers
        self.config['dropout'] = self.dropout
        self.config['description'] = self.description
        self.config['learning_rate'] = self.lr
        self.config['weight_decay'] = self.weight_decay
        self.config['loss_threshold'] = self.loss_threshold
        self.config['model_class'] = self.__class__

        print("Description: {}".format(self.description))
        print("Config: {}".format(self.config))
        with open(self.configfile,'w') as file:
            configdoc = yaml.dump(self.config, file)

    def forward(self, DATA):
        x, edge_index = DATA.x, DATA.edge_index
        x = F.dropout(x, p=self.dropout, training=self.training)
        x = F.relu(self.linear1(x))
        #x = self.conv2(x, edge_index)
        #x = self.conv3(x, edge_index)
        for conv in self.convs:
            x = conv(x, edge_index)

        x = F.dropout(x, p=self.dropout, training=self.training)
        x = self.linear4(x)
        return F.log_softmax(x, dim=1)

    def train_step(self):
        self.train()
        self.optimizer.zero_grad()
        loss = F.nll_loss( self(self.data)[self.data.train_mask], self.data.y[self.data.train_mask] )
        L = loss.item()
        loss.backward()
        self.optimizer.step()
        return L

    def inductive_train_step(self):
        self.train()
        self.optimizer.zero_grad()
        out = self(self.data)
        loss = F.nll_loss(out, self.data.y)
        L = loss.item()
        loss.backward()
        self.optimizer.step()
        return L

    def test_step(self):
        self.eval()
        logits, accs = self(self.data), []

        for _, mask in self.data('train_mask', 'val_mask', 'test_mask'):
            pred = logits[mask].max(1)[1]
            acc = pred.eq(self.data.y[mask]).sum().item() / mask.sum().item()
            accs.append(acc)
        return accs

    def evaluate(self, DATA, mask=None):
        logits = self(DATA)
        pred = None
        if mask is None:
            pred = logits
        else:
            pred = logits[DATA[mask]]

        return pred

    def inductive_test_step(self):
        self.eval()

        if self.query_data is not None:

            out = self(self.query_data)
            predicted = out.argmax(dim=-1)
            correct_prediction =  predicted.eq(self.query_data.y)
            query_accs = correct_prediction.sum().item()/len(self.query_data.y)
        else:
            query_accs = 0.0

        out = self(self.data)
        predicted = out.argmax(dim=-1)
        correct_prediction =  predicted.eq(self.data.y)
        train_accs = correct_prediction.sum().item()/len(self.data.y)

        return train_accs, query_accs
    
    def epoch_training(self, n_epochs = 20000, console_output = False):     
        best_val_acc = test_acc = 0.0
        filename = self.dir + self.description.strip().replace(" ","_") + "_AGNNet" + self.dt + ".csv"
        self._LOGGER.info("Training Information will be saved to {}".format(filename))
        filehandler = open(filename, 'a')
        csvwriter = csv.writer(filehandler)
        csvwriter.writerow(['Loss', 'TrainAccuracy', 'ValidationAccuracy', 'TestAccuracy'])

        for epoch in range(0, n_epochs):
            loss = self.train_step()
            train_acc, val_acc, tmp_test_acc = self.test_step()
            if val_acc > best_val_acc:
                best_val_acc = val_acc
                test_acc = tmp_test_acc
            log = 'Loss: {:.4f}, Train: {:.4f}, Val: {:.4f}, Test: {:.4f}'
            if(console_output):
                self._LOGGER.info(log.format(loss, train_acc, best_val_acc, test_acc))
            csvwriter.writerow([loss, train_acc, best_val_acc, test_acc])

    def loss_training(self, console_output = False):
        
        best_val_acc = test_acc = 0.0

        filename = self.dir + self.description.strip().replace(" ","_") + "_AGNNet_" + self.dt + ".csv"
        filehandler = open(filename, 'a')
        csvwriter = csv.writer(filehandler)
        csvwriter.writerow(['Loss', 'TrainAccuracy', 'ValidationAccuracy', 'TestAccuracy'])

        epoch = 0

        while(True):
            loss = self.train_step()
            train_acc, val_acc, tmp_test_acc = self.test_step()
            if val_acc > best_val_acc:
                best_val_acc = val_acc
                test_acc = tmp_test_acc
            log = 'Loss: {:.4f}, Train: {:.4f}, Val: {:.4f}, Test: {:.4f}'

            if(epoch%100 == 0):
                checkpoint = {}
                checkpoint['epoch'] = epoch
                checkpoint['val_acc'] = best_val_acc
                checkpoint['state_dict'] = self.state_dict()
                checkpoint['optimizer'] = self.optimizer.state_dict()
                self.savecheckpt(checkpoint)

            if(console_output):
                print(log.format(loss, train_acc, best_val_acc, test_acc))
            self._LOGGER.info(log.format(loss, train_acc, best_val_acc, test_acc))
            csvwriter.writerow([loss, train_acc, best_val_acc, test_acc])
            if(loss <= self.loss_threshold):
                checkpoint = {}
                checkpoint['epoch'] = epoch
                checkpoint['val_acc'] = best_val_acc
                checkpoint['state_dict'] = self.state_dict()
                checkpoint['optimizer'] = self.optimizer.state_dict()
                self.savecheckpt(checkpoint)
                break

            if self.earlystopping is not None:
                self.earlystopping(best_val_acc, self)
                print("early_stopping.early_stop: {}".format(self.earlystopping.early_stop))
                if self.earlystopping.early_stop:
                    self._LOGGER.info("Early stopping")
                    break

            epoch = epoch + 1
            if self.scheduler is not None:
                 self.scheduler.step()

    def inductive_loss_training(self, console_output = False):
        
        best_val_acc = test_acc = 0.0

        filename = self.description.strip().replace(" ","_") + "_inductiveAGNNet_" + self.dt + ".csv"
        filehandler = open(filename, 'a')
        csvwriter = csv.writer(filehandler)
        csvwriter.writerow(['Loss', 'TrainAccuracy', 'QueryAccuracy'])

        epoch = 0

        while(True):
            loss = self.inductive_train_step()
            train_acc, query_acc = self.inductive_test_step()
            if query_acc > best_val_acc:
                best_val_acc = query_acc
            log = 'Loss: {:.4f}, Train: {:.4f}, Query: {:.4f}'

            # if self.scheduler is not None:
            #     self.scheduler.step(-1*best_val_acc)

            if(epoch%100 == 0):
                self._LOGGER.info("Epoch= {}. Saving Model Checkpoint".format(epoch))
                checkpoint = {}
                checkpoint['epoch'] = epoch
                checkpoint['val_acc'] = best_val_acc
                checkpoint['state_dict'] = self.state_dict()
                checkpoint['optimizer'] = self.optimizer.state_dict()
                self.savecheckpt(checkpoint)

            if(console_output):
                self._LOGGER.info(log.format(loss, train_acc, best_val_acc))
            csvwriter.writerow([loss, train_acc, best_val_acc])
            
            if(loss <= self.loss_threshold):
                checkpoint = {}
                checkpoint['epoch'] = epoch
                checkpoint['val_acc'] = best_val_acc
                checkpoint['state_dict'] = self.state_dict()
                checkpoint['optimizer'] = self.optimizer.state_dict()
                self.savecheckpt(checkpoint)
                break

            epoch = epoch + 1
            

    def savecheckpt(self, checkpoint):
        torch.save(checkpoint, '{}_{}'.format(self.modelckdir, self.checkpoint_tracker) )
        self.checkpoint_tracker= (self.checkpoint_tracker + 1)%3

    def loadcheckpt(self, oldck_path):
        checkpoint = torch.load(oldck_path)
        self.load_state_dict(checkpoint['state_dict'])
        self.optimizer.load_state_dict(checkpoint['optimizer'])
        return checkpoint['val_acc'].item(), checkpoint['epoch']
