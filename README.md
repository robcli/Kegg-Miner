# Kegg-Miner

## Description

Appends KEGG pathway descriptors to genome annotation. 

## Recommended for multiple files 

Download the lib folder (don't separate the files). To execute in bash, run the .sh file. Input the correct directory with stored prokka/RAST tsvs. Make sure the path starts with a / e.g. /home/user/testing. Then input either "unity" or "kbase". To run as a sbatch file, use the slurm file. Edit the directory and src variables as needed. These methods allow the sequential processing of an entire folder of prokka/RAST tsvs at once. Only works on directories.

## For single use purposes

Run convert.py. Takes RAST or prokka result tsv files. Files can originate from KBase, or from RASTtk or prokka ran manually. To run, have the prokka/RAST output tsv in the same folder as the script. The function to run is convert(filename, type) where type is either "unity" or "kbase". Ex. convert("prokka_phyla.tsv", "unity") if tsv orginated from Unity or convert("prokka_phyla.tsv", "kbase") if the tsv was from kbase. 

## Output files

Pathways is a tsv file with all the products with EC numbers and their associated pathways from KEGG (with exceptions: ECs with no pathways associated).
Counted is a tsv file with the number of times each of the pathways detected from the products appeared.

## Misc
Required dependencies include BeautifulSoup and Requests. Full information in requirements.txt
