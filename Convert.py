import requests
import pandas as pd
import concurrent.futures
import os
import math
import re
import time
import seaborn as sns
from matplotlib import pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("file", type=str,
                    help="prokka file to process in tsv format")
args = parser.parse_args()

class KEGG_obj:
    def __init__(self, pathways, ec):
        self.pathways = pathways
        self.ec = ec
    def get_pathways(self):
        return self.pathways
    def get_ec(self):
        return self.ec

class Data:
    def __init__(self, df, ecs):
        self.df = df
        self.ecs = ecs
        self.kegg_df = None
        self.kegg_plot = None

def read_file(file_name):
    df = pd.read_csv(file_name, sep="\t").applymap(str)
    ecs = list(df[df["EC_number"] != "nan"]["EC_number"])
    return Data(df, ecs)

def search_ec(ec):
    API_dump = requests.get("https://rest.kegg.jp/get/"+ec).content.decode("utf-8")
    pathways = re.findall(r"ec[0-9]{5}.*\n", API_dump)
    pathways_cleaned = list(map(lambda x: x.strip('\n')[9:], pathways))
    pathways_concat = ";".join(pathways_cleaned)
    if pathways_concat == "":
        return KEGG_obj(None, ec)
    return KEGG_obj(pathways_concat, ec)

def mine_kegg(ecs):
    threads = 10
    kegg_map = dict()
    count = 0
    tally = 0
    queue = []
    zero = time.time()

    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        for ec in ecs:
            count += 1
            tally += 1
            queue.append(executor.submit(search_ec, ec))
            if count == 10:
                for future in queue:
                    result = future.result()
                    if result.get_pathways() == None:
                        continue
                    kegg_map[result.get_ec()] = result.get_pathways()
                count = 0 
                print(f"#################################### {tally}/{len(ecs)} | {time.time()-zero} seconds", flush=True)
                time.sleep(1)
    print(f"#################################### {tally}/{len(ecs)} | {time.time()-zero} seconds", flush=True)
    return kegg_map

def append_kegg_id(data, kegg_map):
    data.kegg_df = data.df[data.df["EC_number"].isin(pd.Series(kegg_map.keys()))]
    data.kegg_df = data.kegg_df.assign(pathways = data.kegg_df["EC_number"].map(kegg_map))

def plot_kegg(data):
    fig, ax = plt.subplots(layout="constrained", figsize=(16,12))
    list_of_pathways = list(map(lambda x: x.split(";"), data.kegg_df["pathways"]))
    counts = dict()
    list_of_pathways = sum(list_of_pathways, [])
    for pathway in set(list_of_pathways):
        counts[pathway] = [list_of_pathways.count(pathway)]

    df = pd.DataFrame.from_dict(counts)
    df = df.transpose().reset_index().set_axis(['pathway', 'count'], axis=1).sort_values("count", ascending=False)
    data.kegg_plot = sns.barplot(df, y="pathway", x="count")
    return fig

def main():
    data = read_file(args.file)
    kegg_map = mine_kegg(data.ecs)
    append_kegg_id(data, kegg_map)
    plot_kegg(data).savefig(args.file[:-4]+"_plot.png")
    data.kegg_df.to_csv(args.file[:-4]+"_kegg_id.tsv", sep="\t")
    

main()
