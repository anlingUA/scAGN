# scAGN: Attention-based Graph Neural Network for Label Propagation in Single-Cell Omics

## Requirements
- Ubuntu 18.04 LTS or higher
- Python 3.7.5
- R version 4.1.3

## Installation
We will use Python 3.7.5. First create a python virtual environment.

```
conda create -n gcn python=3.7.5
```

### Pytorch Installation

```
conda install pytorch torchvision torchaudio cudatoolkit=11.3 -c pytorch
```

## Pytorch-Geometric Installation

```
conda install pyg -c pyg -c conda-forge
pip install python-igraph
pip install umap-learn
pip install umap-learn[plot]
```

# Preparation before running the code
1. Clone this repository in your home director

```
cd ~
git clone https://github.com/phycomlab/scAGN
# Create a folder to save intermediate outputs
mkdir gene-cn-output  

# cd ~/scAGN/data
# create folder for saving graph representation
# ./create_dir.sh
```

2. Download dataset. Dataset needs to be download to ~/scAGN/RScript/scData

```
cd  ~/scAGN/RScript/scData
wget <add-url>
./create_dir.sh

```

2. Generate CCA graph

```
cd scAGN/RScript/
./graph_generation.sh
```

3. Run scAGN on all datasets

```
cd ~
~/scAGN/python/graph_AGN.py -i 128 -l 5 -d 0.0 -t 0.05 -r 0.1 -w 0.0000 -s 10 -G 0.8 -D GSE108989 -P 2000 > ~/scAGN/script/GRAPH_agn_4_GSE108989$now.txt
```

4. Run baseline methods

```
cd ~/scAGN/RScript
./do_CHETAH_prediction.sh
./do_scMap_prediction.sh
./do_seurat_prediction.sh
./do_singleR_prediction.sh
```

5. Generate plots used in the paper

```
cd ~/scAGN/python
python analysis.py
```

## Licensing

    License: MIT License 
    Copyright Rahul Bhadani, Lingling An
    Initial Date: Apr 16, 2022
    Permission is hereby granted, free of charge, to any person obtaining 
    a copy of this software and associated documentation files 
    (the "Software"), to deal in the Software without restriction, including
    without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to 
    permit persons to whom the Software is furnished to do so, subject 
    to the following conditions:

    The above copyright notice and this permission notice shall be 
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF 
    ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED 
    TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT 
    SHALL THE AUTHORS, COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
    DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
    CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF 
    OR IN CONNECTION WITH THE SOFTWARE OR THE USE 
    OR OTHER DEALINGS IN THE SOFTWARE.
