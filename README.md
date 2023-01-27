# Kegg-Miner

Takes RAST or prokka result tsv files. Files can originate from KBase, or from RASTtk or prokka ran manually. To run, have the prokka/RAST output tsv in the same folder as the script. The function to run is convert(filename, type) where type is either "unity" or "kbase". Ex. convert("prokka_phyla.tsv", "unity") if tsv orginated from Unity or convert("prokka_phyla.tsv", "kbase") if the tsv was from kbase. 

Pathways is a tsv file with all the products with EC numbers and their associated pathways from KEGG (with exceptions: ECs with no pathways associated).

Counted is a tsv file with the number of times each of the pathways detected from the products appeared.

Required dependencies include BeautifulSoup and Requests. Full information in requirements.txt
