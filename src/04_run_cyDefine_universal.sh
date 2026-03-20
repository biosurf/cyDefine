#!/bin/bash

datasets=("POISED")
echo "Running cyDefine"
for dataset in "${datasets[@]}"; do
    echo "Processing dataset: $dataset"
    Rscript R/universal_transfer_cyDefine.R \
        --name "$dataset" \
        --seed 145 \
        --threads 39

done
