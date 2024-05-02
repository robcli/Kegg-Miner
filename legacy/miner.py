#!/usr/bin/env python


import requests
from collections import defaultdict
import csv
from bs4 import BeautifulSoup
import concurrent.futures
from concurrent.futures import as_completed
import os

#Converting annotation to uniform throughput

def kbase_protein(file):
    i = [2, 8]
    r=[]
    with open(file, "r", encoding="UTF-8") as read:
        tsv_file = csv.reader(read, delimiter = '\t')
        for line in tsv_file:
            out = [line[x] for x in i]
            out.insert(1, extract_EC(out[1]))
            out[2] = extract_function(out[2])
            if("hypothetical protein" not in out[2] and out[2] != ""):
                r.append(out)
    return r

def unity_protein(file):
    i = [1, 4, 6]
    r=[]
    with open(file, "r", encoding="UTF-8") as read:
        tsv_file = csv.reader(read, delimiter = '\t')
        for line in tsv_file:
            out = [line[x] for x in i]
            if("hypothetical protein" not in out[2]):
                r.append(out)
        r.pop(0)
    return r

def extract_EC(str):
    if ("eC_number=" in str):
        index = str.index("eC_number=")+10
        end = str[index:].index(";")
        return str[index:index+end].replace("%2", ",")
    elif ("(EC" in str):
        index = str.index("(EC")+4
        end = str[index:].index(")")
        return str[index:index+end]
    else:
        return ""

def extract_function(str):
    index = str.index("product=")+8
    present = ";" in str[index:]
    end = str[index:].index(";") if present else -1
    return str[index:index+end].replace("%2", ",") if end != -1 else str[index:].replace("%2", ",")

def protein_to_tsv(file, program):
    oname = program + "_" + file[:-3] + "tsv"
    with open (oname, 'w', newline ='') as csvfile:
        w = csv.writer(csvfile, delimiter = '\t')
        rows = unity_protein(file) if program == "unity" else kbase_protein(file)
        rows.insert(0, ["Type", "EC", "Protein"])
        w.writerows(rows)
    return oname


##Scraping methods

def pathway(line):
    ec = line[1]
    URL = "https://www.genome.jp/dbget-bin/www_bget?ec:"+ec
    page = requests.get(URL)
    soup = BeautifulSoup(page.content, "html.parser")
    try:
        pathway = soup.find_all("td", {"class":"td20 defd"})[5]
    except:
        return line + ["Obsolete"]
    if pathway.select_one('a[href*="/pathway"]') is None:
        return line + ["No metabolic pathway"]
    else:
        str = ""
        rows = pathway.findAll(lambda tag: tag.name=='table')
        for elem in rows:
            str += elem.findAll(lambda tag: tag.name =='td')[1].get_text() + ";"
    return line + [str]

def ec_search(file):
    with open(file, 'r', encoding="UTF-8") as read:
        threads = 30
        files_to_search = []
        futures = []
        rr = []
        tsv_file = csv.reader(read, delimiter = "\t")
        for line in tsv_file:
            ec = line[1]
            if ec != None and ec != "" and ec != "EC" and "-" not in ec:
                files_to_search.append(line)
        update = 0
        print("Beginning to mine for " + file, flush=True)

        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            for i, line in enumerate(files_to_search):
                files_to_search[i] = executor.submit(pathway, line)
                

            for i, future in enumerate(files_to_search):
                files_to_search[i] = future.result()
                update += 1
                if(update % 20 == 0):
                    print(str(update) + "/" + str(len(files_to_search)), flush=True)
            return files_to_search
    
def file_overwrite(file):
    overwrite = ec_search(file)
    overwrite.insert(0, ["Type", "EC", "Protein", "Pathway"])
    with open(file[:-4]+"_pathways.tsv", "w", newline = '') as csvfile:
        write = csv.writer(csvfile, delimiter = '\t')
        write.writerows(overwrite)
    return file[:-4]+"_pathways.tsv"

## Produce count file

def count(file):
    with open(file, encoding="UTF-8") as read:
        tsv_file = csv.reader(read, delimiter = '\t')
        acc = []
        for line in tsv_file:
            acc.append(pull_pathways(line[3]))
        dd = defaultdict(int)
        for line in acc:
            for path in line:
                dd[path] += 1
        return dd 

def pull_pathways(str):
    alpha = ""
    r = []
    for letter in str:
        if(letter == ";"):
             r.append(alpha)
             alpha = ""
        else:
            alpha += letter
    if("Metabolic pathways" in r): 
        r.remove("Metabolic pathways")       
    return r 

def count_to_tsv(file):
    p = list(count(file).items())
    oname = file[:-12] + "counted.tsv"
    with open (oname, 'w', newline ='') as csvfile:
        w = csv.writer(csvfile, delimiter = '\t')
        w.writerow(["Pathway", "Count"])
        w.writerows(p)
    return oname

## Combined method

def convert(file, program):
    name = protein_to_tsv(file, program)
    name2 = file_overwrite(name)
    count_to_tsv(name2)
    print("Deleting intermediate file " + name, flush=True)
    os.remove(name)

convert(os.environ['f'], os.environ['s'])
