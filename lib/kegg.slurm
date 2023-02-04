#!/bin/bash
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --mem=32768  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G 1  # Number of GPUs
#SBATCH -t 06:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = job ID

chmod +x miner.py

mineloc=$(dirname "$(realpath $0)")/"miner.py"
directory="/home/robcli_umass_edu/Post_annotation/testing"
src="unity"

cd $directory
arr=()

for file in $directory/*; do

    base="${file##*/}"
    arr+=("${base}")
    
    export f="$base"
    export s="$src"
    
    mineloc=$(dirname "$(realpath $0)")/"miner.py"
    
    mkdir "${base%.tsv}"
    mv "${src}_${base%.tsv}_pathways.tsv" $directory/${base%.tsv}
    mv "${src}_${base%.tsv}_counted.tsv" $directory/${base%.tsv}

done
