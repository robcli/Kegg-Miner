#!/bin/bash

chmod +x miner.py

mineloc=$(dirname "$(realpath $0)")/"miner.py"

read -p "Enter directory path: " directory
read -p "Input file source (unity or kbase): " src

cd $directory || exit
arr=()

for file in "$directory"/*; do

    base="${file##*/}"
    arr+=("${base}")
    
    export f="$base"
    export s="$src"
    
    $mineloc
    
    mkdir "${base%.tsv}"
    mv "${src}_${base%.tsv}_pathways.tsv" "$directory"/"${base%.tsv}"
    mv "${src}_${base%.tsv}_counted.tsv" "$directory"/"${base%.tsv}"

done
