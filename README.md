# Kegg-Miner

## Description

Appends KEGG pathway descriptors to prokka tsv files. 

## How to Use
Execute as:

python Convert.py filename

Ex: python Convert.py A1_bin.012_prokka.tsv

Example prokka files to use are in the data folder. 

## Output
A tsv of only prokka entries with a KEGG id with pathway description.

A png of a plot of KEGG pathway appears by frequency. 

![Example plot](https://github.com/robcli/Kegg-Miner/blob/main/A1_bin.012_prokka_plot.png)

## Misc
Required dependencies include Requests, Seaborn, and Pandas. Full information in requirements.txt
