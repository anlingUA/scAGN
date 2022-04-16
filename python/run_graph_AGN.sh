#!/bin/bash
./graph_AGN108989.py -i 16 -d 0.1 -t 0.01 &
./graph_AGN108989.py -i 8 -d 0.1 -t 0.01 &
./graph_AGN108989.py -n -i 16 -d 0.5 -t 0.001 -r 0.00001 &
./graph_AGN108989.py -n -i 16 -d 0.5 -t 0.001 -r 0.000001 &
./graph_AGN108989.py -n -i 256 -d 0.5 -t 0.001 -r 0.00001 -w 0 &
./graph_AGN108989.py -n -i 512 -d 0.5 -t 0.001 -r 0.000001 -w 0 &
./graph_AGN108989.py -n -I -i 512 -d 0.5 -t 0.001 -r 0.000001 -w 0 &

./graph_AGN108989_umap.py -n -I -i 512 -d 0.5 -t 0.001 -r 0.000001 -w 0.000001 &
./graph_AGN108989_umap.py -n -I -i 512 -d 0.5 -t 0.001 -r 0.001 -w 0.000001 &

./graph_AGN108989_umap.py -I -i 512 -d 0.5 -t 0.001 -r 0.000001 -w 0.000001 -g knn &
./graph_AGN108989_umap.py -I -i 512 -d 0.5 -t 0.001 -r 0.000001 -w 0.000001 -g umap &

./graph_AGN108989_umap.py -I -i 512 -d 0.5 -t 0.001 -r 0.000001 -w 0.000001 -g knn &

./graph_AGN.py -I -i 512 -d 0.5 -t 0.001 -r 0.000001 -w 0.000001 -g knn -D GSE118389 &
./graph_AGN.py -I -i 512 -d 0.5 -t 0.001 -r 0.0000001 -w 0.000001 -g knn -D GSE99254 &

# Works well with sample data
 ./graph_AGN.py -i 128 -l 3 -d 0.0 -t 0.05 -r 0.1 -w 0.0000 -D GSE84133 -P 100000 -p


# Tried last time with GSE108989
#./graph_GCN.py -i [64,32,64] -l 4 -d 0.0 -r 0.0001 -w 0.0 -t 0.24 -P 100000 -p

./graph_AGN.py -i 128 -l 3 -d 0.0 -t 0.05 -r 0.1 -w 0.0000 -s 10 -G 0.8 -D GSE108989 -P 100000 -p

./graph_AGN.py -i 256 -l 3 -d 0.0 -t 0.05 -r 0.1 -w 0.0000 -s 10 -G 0.75 -D GSE108989 -P 100000 -p

./graph_AGN.py -i 256 -l 3 -d 0.0 -t 0.05 -r 0.1 -w 0.0000 -s 10 -G 0.55 -D GSE108989 -P 100000 -p

./graph_AGN.py -i 256 -l 3 -d 0.0 -t 0.05 -r 0.5 -w 0.0000 -s 10 -G 0.5 -D GSE108989 -P 100000 -p

./graph_AGN.py -i 128 -l 4 -d 0.0 -t 0.05 -r 0.1 -w 0.0000 -s 10 -G 0.8 -D GSE108989 -P 100000 -p


# works well with GSE72056
./graph_AGN.py -i 128 -l 4 -d 0.0 -t 0.05 -r 0.1 -w 0.0000 -s 10 -G 0.8 -D GSE72056 -P 100000 -p
./graph_AGN.py -i 128 -l 3 -d 0.0 -t 0.01 -r 0.1 -w 0.0000 -s 10 -G 0.99 -D GSE72056 -P 100000 -p # best result so far
./graph_AGN.py -i 128 -l 3 -d 0.10 -t 0.001 -r 0.1 -w 0.0000 -s 10 -G 0.999 -D GSE72056 -P 100000 -p # best result so far

# works well with GSE85241
./graph_AGN.py -i 128 -l 3 -d 0.50 -t 0.05 -r 0.1 -w 0.0000 -s 10 -G 0.999  -D GSE85241 -P 100000 -p

# works well with GSE98638
 ./graph_AGN.py -i 256 -l 4 -d 0.50 -t 0.001 -r 0.1 -w 0.000001 -s 10 -G 0.99 -D GSE98638 -P 100000 -p # Run this on HPC with samepleType as label

# works well with GSE118389
 ./graph_AGN.py -i 128 -l 4 -d 0.50 -t 0.05 -r 0.1 -w 0.0001 -s 10 -G 0.8 -D GSE118389 -P 100000 -p
