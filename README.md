# Kegg-Miner

Takes RAST or prokka result tsv files. Files can originate from KBase, or from RASTtk or prokka ran manually. 

Pathways is a tsv file with all the products with EC numbers and their associated pathways from KEGG (with exceptions: ECs with no pathways associated).

Counted is a tsv file with the number of times each of the pathways detected from the products appeared.

Required dependencies include BeautifulSoup and Requests. Full information in requirements.txt
