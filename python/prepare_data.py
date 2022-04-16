import os
import torch
import numpy as np
import random
import pandas as pd
import time as tm
from operator import index, itemgetter
from sklearn.model_selection import train_test_split
import pickle as pkl
from collections import defaultdict
from scipy import sparse as sp
import torch_geometric
from torch_geometric.data import Data

# from graph import *

def get_value(diction, specific):
    for key, val in diction.items():
        if val == specific:
            return (key)

def graph(matrix):
    adj = defaultdict(list)  # default value of int is 0
    for i, row in enumerate(matrix):
        for j, adjacent in enumerate(row):
            if adjacent:
                adj[i].append(j)
        if adj[i].__len__ == 0:
            adj[i] = []
    return adj

def sample_mask(idx, l):
    """Create mask."""
    mask = np.zeros(l)
    mask[idx] = 1
    return np.array(mask, dtype=np.bool)

def normalize_adj(adj):
    """Symmetrically normalize adjacency matrix."""
    adj = sp.coo_matrix(adj)
    rowsum = np.array(adj.sum(1))
    d_inv_sqrt = np.power(rowsum, -0.5).flatten()
    d_inv_sqrt[np.isinf(d_inv_sqrt)] = 0.
    d_mat_inv_sqrt = sp.diags(d_inv_sqrt)
    return adj.dot(d_mat_inv_sqrt).transpose().dot(d_mat_inv_sqrt).tocoo()

def sparse_to_tuple(sparse_mx):
    """Convert sparse matrix to tuple representation."""
    def to_tuple(mx):
        if not sp.isspmatrix_coo(mx):
            mx = mx.tocoo()
        coords = np.vstack((mx.row, mx.col)).transpose()
        values = mx.data
        shape = mx.shape
        return coords, values, shape

    if isinstance(sparse_mx, list):
        for i in range(len(sparse_mx)):
            sparse_mx[i] = to_tuple(sparse_mx[i])
    else:
        sparse_mx = to_tuple(sparse_mx)

    return sparse_mx

#' data preperation
def generate_dataset(DataDir, **kwargs):
    """
    Dataset preparation before training
    What does this function do:
    1. Split the dataset into equal parts per types so that there is no unbalance in the training
    2. 
    """
    print('DataDir = {}'.format(DataDir))
    # if Rgraph==False:
    #     graph_construct(outputdir='process_data')

    datafile1 = kwargs.get("Data1", "Data1.csv") # Reference dataset
    datafile2 = kwargs.get("Data2", "Data2.csv") # Query dataset
    labelfile1 = kwargs.get("Label1", "Label1.csv") # Labels for reference
    labelfile2 = kwargs.get("Label2", "Label2.csv") # Labels for query Data
 
    # DataDir = '/home/starfire/VersionControl/SingleCellGraph/scGCN/input'

    datafile1 = "Data1.csv"
    datafile2 = "Data2.csv"
    labelfile1 = "Label1.csv"
    labelfile2 = "Label2.csv"
   
    DataPath1 = '{}/{}'.format(DataDir, datafile1)
    DataPath2 = '{}/{}'.format(DataDir,datafile2)
    LabelsPath1 = '{}/{}'.format(DataDir, labelfile1)
    LabelsPath2 = '{}/{}'.format(DataDir, labelfile2)

    #' read the data
    data_df1 = pd.read_csv(DataPath1, index_col=0, sep=',')
    data_df2 = pd.read_csv(DataPath2, index_col=0, sep=',')
    print(LabelsPath1)
    label_df1 = pd.read_csv(LabelsPath1, header=0, index_col=False, sep='\n')
    label_df2 = pd.read_csv(LabelsPath2, header=0, index_col=False, sep='\n')

    data_df1 = data_df1.reset_index(drop=True)  #.transpose()
    data_df2 = data_df2.reset_index(drop=True)  #.transpose()
    label_df1.columns = ['type']
    label_df2.columns = ['type']

    types = np.unique(label_df1['type']).tolist()

    random.seed(123)
    p_data = []
    p_label = []

    # First split dataset according types (or unique cell labels)
    # and then randomize the data and then add to a list of dataframes
    for i in types:
        indices_of_type = label_df1[label_df1['type'] == i].index
        label_of_type = label_df1[label_df1['type'] == i]
        data_per_type_df = data_df1.iloc[indices_of_type]
        num_to_select = len(data_per_type_df)
        random_items = random.sample(range(0, len(indices_of_type)), num_to_select) # Just randomize the indices
        sub_data = data_per_type_df.iloc[random_items]
        sub_label = label_of_type.iloc[random_items]
        p_data.append(sub_data)
        p_label.append(sub_label)
    # split data to training, test, valdiaton sets

    data_train = []
    data_test = []
    data_validation = []
    label_train = []
    label_test = []
    label_val = []

    for i in range(0, len(p_data)):
        train_df, test_df, train_label, test_label = train_test_split(p_data[i], p_label[i], test_size=0.1, random_state=1)        
        train_df, validation_df, teml_train, validation_label = train_test_split(train_df, train_label, test_size=0.1, random_state=1)
        data_train.append(train_df)
        label_train.append(teml_train)
        data_test.append(test_df)
        label_test.append(test_label)
        data_validation.append(validation_df)
        label_val.append(validation_label)

    # each of training, test and validation dataset has equal amount of each type
    data_train_df = pd.concat(data_train)
    data_test_df = pd.concat(data_test)
    data_validation_df = pd.concat(data_validation)
    label_train_df = pd.concat(label_train)
    label_test_df = pd.concat(label_test)
    label_validation_df = pd.concat(label_val)

    ###
    ## At this point, we have feature and labels dataset of Dataset1 as split into three parts: training, testing and validation.
    ## 
    ## Now, we combine query data set to above training dataset.
    ###

    combined_df = pd.concat([data_train_df, data_df2])
    combined_label = pd.concat([label_train_df, label_df2])

    index_guide = np.concatenate(
        (label_train_df.index, label_df2.index * (-1) - 1, label_validation_df.index,
         label_test_df.index))

    #' 3) set the unlabeled data in training set

    #' @param N; the number of labeled samples in training set
    M = len(data_train_df)

    #' 4) get the feature object by combining training, test, valiation sets    
    # Row-normalize feature matrix and convert to tuple representation
    features = pd.concat([combined_df,data_validation_df, data_test_df])
    r_inv = 1/np.array(features.sum(1))
    r_inv[np.isinf(r_inv)] = 0.
    r_mat_inv = sp.diags(r_inv)
    features_ = r_mat_inv.dot(features)
    features = pd.DataFrame(features_, columns=features.columns, index=features.index)

     #' 5) Given cell type, generate three sets of labels with the same dimension
    labels_train = combined_label['type'].values # np.array(combined_label).flatten()
    labels_test = label_test_df['type'].values # np.array(label_test_df).flatten()
    labels_validation = label_validation_df['type'].values # np.array(label_validation_df).flatten()

    Labels = pd.concat([combined_label, label_validation_df, label_test_df])

    label_mapping = {}
    print(types)
    for line in range(0, len(types)):
        key = types[line]
        label_mapping[key] = int(line)

    labelmap_df = pd.DataFrame(label_mapping.items(), columns = ['label_names', 'id'])
    labelmap_df.to_csv(DataDir+'/label_mapping.csv')

    print(label_mapping)
    Label1 = Labels.replace(label_mapping)
    indices = np.array(Label1.values, dtype='int').tolist()

    indice = np.array(Label1.values, dtype='int').flatten().tolist() #[item for sublist in indices for item in sublist]

    idx_train = range(M)
    idx_pred = range(M, len(labels_train))
    idx_val = range(len(labels_train), len(labels_train) + len(labels_validation))
    idx_test = range(
        len(labels_train) + len(labels_validation),
        len(labels_train) + len(labels_validation) + len(labels_test))

    train_mask = sample_mask(idx_train, features.shape[0])
    pred_mask = sample_mask(idx_pred, features.shape[0])
    val_mask = sample_mask(idx_val, features.shape[0])
    test_mask = sample_mask(idx_test, features.shape[0])


    print('load data succesfully....')

    id_graph1 = pd.read_csv('{}/inter_graph.csv'.format(DataDir),
                            index_col=0,
                            sep=',')
    id_graph2 = pd.read_csv('{}/intra_graph.csv'.format(DataDir),
                            sep=',',
                            index_col=0)

    fake1 = np.array([-1] * len(data_df2.index))
    index1 = np.concatenate((data_train_df.index, fake1, data_validation_df.index,
                             data_test_df.index)).flatten()
    #' (feature_data.index==index1).all()
    fake2 = np.array([-1] * len(data_train_df))
    fake3 = np.array([-1] * (len(data_validation_df) + len(data_test_df)))
    find1 = np.concatenate((fake2, np.array(data_df2.index), fake3)).flatten()



    id_grp1 = np.array([
        np.concatenate((np.where(find1 == id_graph2.iloc[i, 1])[0],
                        np.where(find1 == id_graph2.iloc[i, 0])[0]))
        for i in range(len(id_graph2))
    ])

    id_grp2 = np.array([
        np.concatenate((np.where(find1 == id_graph2.iloc[i, 0])[0],
                        np.where(find1 == id_graph2.iloc[i, 1])[0]))
        for i in range(len(id_graph2))
    ])

    # Basically id_grp2 and id_grp1 are reverse edges of each other. Combining
    # them, we can have undirected or bidirected graph


    #' ---------------------------------------------
    #'  inter-graph
    #' ---------------------------------------------
    id_gp1 = np.array([
        np.concatenate((np.where(find1 == id_graph1.iloc[i, 1])[0],
                        np.where(index1 == id_graph1.iloc[i, 0])[0]))
        for i in range(len(id_graph1))
    ])

    id_gp2 = np.array([
        np.concatenate((np.where(index1 == id_graph1.iloc[i, 0])[0],
                        np.where(find1 == id_graph1.iloc[i, 1])[0]))
        for i in range(len(id_graph1))
    ])
    

    id_gp1 = pd.DataFrame(id_gp1)
    id_gp1.columns = ['V1', 'V2']

    id_gp2 = pd.DataFrame(id_gp2)
    id_gp2.columns = ['V1', 'V2']

    id_grp1 = pd.DataFrame(id_grp1)
    id_grp1.columns = ['V1', 'V2']

    id_grp2 = pd.DataFrame(id_grp2)
    id_grp2.columns = ['V1', 'V2']


    graph_set = pd.concat([id_gp1, id_gp2, id_grp1, id_grp2 ])

    X = features.values.tolist()
    x = torch.tensor(X, dtype=torch.float)
    Y = Label1['type'].tolist()
    y = torch.tensor(Y, dtype=torch.long)
    edge_index_tensor = torch.tensor(graph_set.values.tolist(), dtype=torch.long)
    data = Data(x=x, edge_index=edge_index_tensor.t().contiguous(), y=y)

    data.train_mask = torch.tensor(train_mask.tolist())
    data.test_mask = torch.tensor(test_mask.tolist())
    data.val_mask = torch.tensor(val_mask.tolist())
    data.pred_mask = torch.tensor(pred_mask.tolist())

    return data
    
    # return features, graph_set, Label1, train_mask, pred_mask, val_mask, test_mask



