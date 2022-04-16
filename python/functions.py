#!/usr/bin/env python

from inspect import stack
import logging
import sys
import igraph as ig
from sklearn.neighbors import kneighbors_graph
from sklearn.neighbors import radius_neighbors_graph
import umap
import pandas as pd
from sklearn.cross_decomposition import CCA
from sklearn.preprocessing import StandardScaler, normalize
import scipy as sp
from sklearn.preprocessing import StandardScaler, normalize
from itertools import chain
import numpy as np

def configure_logworker(level = logging.INFO, logfile = None):
    """
    Logging Configuration
    """

    logworker = logging.getLogger('')
    logworker.setLevel(level)
    logworker.propagate = False
    logworker.handlers = []

    log_format = '[%(asctime)s] (%(name)s) %(levelname)s: %(message)s'
    log_date_format = '%Y_%m_%d_%H_%M_%S'
    formatter = logging.Formatter(log_format, log_date_format)


    if logfile is None:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(formatter)
        logworker.addHandler(handler)
    else:
        handler = logging.FileHandler(logfile)
        handler.setFormatter(formatter)
        logworker.addHandler(handler)

    return logworker

def unit_normalize(data, norm="l2", bySample=True):
	"""
	Default norm used is l2-norm. Other options: "l1", and "max"
	If bySample==True, then we independently normalize each sample. If bySample==False, then we independently normalize each feature
	"""
	assert (norm in ["l1","l2","max"]), "Norm argument has to be either one of 'max', 'l1', or 'l2'."

	if bySample==True:
		axis=1
	else:
		axis=0

	return normalize(data, norm=norm, axis=axis) 

def zscore_standardize(data):
	scaler=StandardScaler()
	scaledData=scaler.fit_transform(data)
	return scaledData

def get_spatial_distance_matrix(data, metric="eucledian"):
	Cdata= sp.spatial.distance.cdist(data,data,metric=metric)
	return Cdata/Cdata.max()

def geodesic_distances(X, num_neighbors, mode="distance", metric="minkowski"):

    assert (mode in ["connectivity", "distance"]), "Norm argument has to be either one of 'connectivity', or 'distance'. "
    if mode=="connectivity":
        include_self=True
    else:
        include_self=False
    knn = kneighbors_graph(X, num_neighbors, n_jobs=-1, mode=mode, metric=metric, include_self=include_self)
    connected_components = sp.csgraph.connected_components(knn, directed=False)[0]
    dist = sp.csgraph.dijkstra(knn, directed=False)
    connected_element = []

    ## for local connectively
    if connected_components is not 1:
        inf_matrix = []
        
        for i in range(len(X)):
            inf_matrix.append(list(chain.from_iterable(np.argwhere(np.isinf(dist[i])))))

        for i in range(len(X)):
            if i==0:
                connected_element.append([0])
            else:
                for j in range(len(connected_element)+1):
                    if j == len(connected_element):
                        connected_element.append([])
                        connected_element[j].append(i)
                        break
                    if inf_matrix[i] == inf_matrix[connected_element[j][0]]:
                        connected_element[j].append(i)
                        break

        components_dist = []
        x_index = []
        y_index = []
        components_dist.append(np.inf)
        x_index.append(-1)
        y_index.append(-1)
        for i in range(connected_components):
            for j in range(i):    
                for num1 in connected_element[i]:
                    for num2 in connected_element[j]:
                        if np.linalg.norm(X[num1]-X[num2])<components_dist[len(components_dist)-1]:
                            components_dist[len(components_dist)-1]=np.linalg.norm(X[num1]-X[num2])
                            x_index[len(x_index)-1] = num1
                            y_index[len(y_index)-1] = num2
                components_dist.append(np.inf)
                x_index.append(-1)
                y_index.append(-1)

        components_dist = components_dist[:-1]
        x_index = x_index[:-1]
        y_index = y_index[:-1]

        sort_index = np.argsort(components_dist)
        components_dist = np.array(components_dist)[sort_index]
        x_index = np.array(x_index)[sort_index]
        y_index = np.array(y_index)[sort_index]

        for i in range(len(x_index)):
            knn = knn.todense()
            knn = np.array(knn)
            knn[x_index[i]][y_index[i]] = components_dist[i]
            knn[y_index[i]][x_index[i]] = components_dist[i]
            knn = sp.csr_matrix(knn)
            connected_components = sp.csgraph.connected_components(knn, directed=False)[0]
            dist = sp.csgraph.dijkstra(knn, directed=False)
            if connected_components == 1:
                break

    return dist/dist.max()

def construct_knngraph(values, n_neighbors = 30, metric = 'euclidean', stacked = False):
    """
    `construct_knn` constructs k-nearest neighbor graph.

    Parameters
    -----------
    values: `numpy.array`
        Array with gene-features and cell/samples as columns
    
    n_neighbors: `int`
        The value `k` in nearest neighbor
    
    metric: 
    """

    # cell_ids = df.index
    knn_scRNAseq = kneighbors_graph(values, n_neighbors, metric = metric,  mode = 'connectivity', include_self=True).toarray()

    if stacked:
        #knn_scRNAseq = pd.DataFrame(knn_scRNAseq, columns = cell_ids, index = cell_ids)
        knn_scRNAseq = pd.DataFrame(knn_scRNAseq)
        knn_scRNAseq = knn_scRNAseq.stack().reset_index()
        knn_scRNAseq = knn_scRNAseq.rename(columns = {'level_0': 'V1', 'level_1': 'V2', 0: 'connectivity'})
        knn_scRNAseq = knn_scRNAseq.loc[knn_scRNAseq['connectivity'] != 0]
        knn_scRNAseq = knn_scRNAseq.drop(columns='connectivity')
        knn_scRNAseq = knn_scRNAseq.reset_index()
    # if file_to_save is not None:
    #     knn_scRNAseq.to_csv(file_to_save)

    return knn_scRNAseq

def construct_radiusgraph(values, radius = 1.5, stacked = False):
    """
    `construct_radiusgraph` constructs k-nearest neighbor graph.

    Parameters
    -----------
    df: `pandas.DataFrame`
        Dataframe with gene-features and cell/samples as columns
    
    n_neighbors: `int`
        The value `k` in nearest neighbor
    
    metric: 
    """

    #cell_ids = df.index
    radius_scRNAseq = radius_neighbors_graph(values, radius=radius, mode='connectivity',  include_self=True)    
    if stacked:
        radius_scRNAseq = pd.DataFrame(radius_scRNAseq)
        radius_scRNAseq = radius_scRNAseq.stack().reset_index()
        radius_scRNAseq = radius_scRNAseq.rename(columns = {'level_0': 'V1', 'level_1': 'V2', 0: 'connectivity'})
        radius_scRNAseq = radius_scRNAseq.loc[radius_scRNAseq['connectivity'] != 0]
        radius_scRNAseq = radius_scRNAseq.drop(columns='connectivity')
        radius_scRNAseq = radius_scRNAseq.reset_index()

    # radius_scRNAseq.to_csv(file_to_save)
    return radius_scRNAseq

def construct_umapgraph(values, stacked = False):
    """
    `construct_umapgraph` constructs umap neighbor graph.

    Parameters
    -----------
    df: `pandas.DataFrame`
        Dataframe with gene-features and cell/samples as columns
    
    n_neighbors: `int`
        The value `k` in nearest neighbor
    
    metric: 
    """
    mapper = umap.UMAP().fit(values)
    adjacency_matrix = mapper.graph_.toarray()

    if stacked:
        adjacency_matrix = pd.DataFrame(adjacency_matrix)
        adjacency_matrix = adjacency_matrix.stack().reset_index()
        adjacency_matrix = adjacency_matrix.rename(columns = {'level_0': 'V1', 'level_1': 'V2', 0: 'connectivity'})
        adjacency_matrix = adjacency_matrix.loc[adjacency_matrix['connectivity'] > 0.8]
        adjacency_matrix['connectivity'] = 1.0
        adjacency_matrix = adjacency_matrix.drop(columns='connectivity')
        adjacency_matrix = adjacency_matrix.reset_index()
    # adjacency_matrix.to_csv(file_to_save)
    return adjacency_matrix

def construct_pearsongraph(values, stacked = False):
    """
    """

    df = pd.DataFrame(values)

    adjacency_matrix = df.T.corr(method='pearson')

    if stacked:
        adjacency_matrix = adjacency_matrix.stack().reset_index()
        adjacency_matrix = adjacency_matrix.rename(columns = {'level_0': 'V1', 'level_1': 'V2', 0: 'connectivity'})
        adjacency_matrix = adjacency_matrix.loc[adjacency_matrix['connectivity'] > 0.8]
        adjacency_matrix['connectivity'] = 1.0
        adjacency_matrix = adjacency_matrix.drop(columns='connectivity')
        adjacency_matrix = adjacency_matrix.reset_index()
    # adjacency_matrix.to_csv(file_to_save)

    if not stacked:
        return adjacency_matrix.values
    return adjacency_matrix

def construct_spearmangraph(values, stacked = False):
    """
    """
    df = pd.DataFrame(values)

    adjacency_matrix = df.T.corr(method='spearman')

    if stacked:
        adjacency_matrix = adjacency_matrix.stack().reset_index()
        adjacency_matrix = adjacency_matrix.rename(columns = {'level_0': 'V1', 'level_1': 'V2', 0: 'connectivity'})
        adjacency_matrix = adjacency_matrix.loc[adjacency_matrix['connectivity'] > 0.8]
        adjacency_matrix['connectivity'] = 1.0
        adjacency_matrix = adjacency_matrix.drop(columns='connectivity')
        adjacency_matrix = adjacency_matrix.reset_index()
    # adjacency_matrix.to_csv(file_to_save)

    if not stacked:
        return adjacency_matrix.values
    return adjacency_matrix

def construct_kendallgraph(values, stacked = False):
    """
    """
    df = pd.DataFrame(values)

    adjacency_matrix = df.T.corr(method='kendall')

    if stacked:
        adjacency_matrix = adjacency_matrix.stack().reset_index()
        adjacency_matrix = adjacency_matrix.rename(columns = {'level_0': 'V1', 'level_1': 'V2', 0: 'connectivity'})
        adjacency_matrix = adjacency_matrix.loc[adjacency_matrix['connectivity'] > 0.8]
        adjacency_matrix['connectivity'] = 1.0
        adjacency_matrix = adjacency_matrix.drop(columns='connectivity')
        adjacency_matrix = adjacency_matrix.reset_index()
        # adjacency_matrix.to_csv(file_to_save)

    if not stacked:
        return adjacency_matrix.values
    return adjacency_matrix

def GW_adjmatrix(values, stacked=False):
    pass


def construct_eigenmapgraph(df, file_to_save):
    pass

def construct_tsnegraph(df, file_to_save):
    pass

def construct_ccagraph(df, file_to_save):
    pass

def make_unique(mylist):
    """
    Make all entries of the list unique
    """
    pass
    
def dict_reverser(d):
    seen = set()
    return {v: k for k, v in d.items() if v not in seen or seen.add(v)}

def splitby(entry, delimiter='_', return_indices=1):
    """
    Split a string by delimiter and return `return_indices` entry
    """
    return entry.split(delimiter)[return_indices]