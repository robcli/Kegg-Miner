import requests
from collections import defaultdict
import csv
from bs4 import BeautifulSoup
import concurrent.futures
from concurrent.futures import as_completed
import os

#Converting annotation to uniform throughput
## automatically removes hypothetical protein readings or empty readings
## %2 are converted to commas per encoding
## all file are passed by their name

## File[contig, program, ftype, start, end, score, strand, properties] => [ftype, EC, product]
def kbase_protein(file):
    i = [2, 8]      ##2 is ftype, 8 is product 
    r=[]
    with open(file, "r", encoding="UTF-8") as read:
        tsv_file = csv.reader(read, delimiter = '\t')
        for line in tsv_file:
            out = [line[x] for x in i]      ##out = [ftype, product]
            out.insert(1, extract_EC(out[1]))       ##out = [ftype, EC, product]
            out[2] = extract_function(out[2])       ## runs scrapper
            if("hypothetical protein" not in out[2] and out[2] != ""):
                r.append(out)
    return r

## File[locus_tag, ftype, length_bp, gene, EC_number, COG, product] => [ftype, EC, product]
def unity_protein(file):
    i = [1, 4, 6]       ##1 is ftype, 4 is EC, 6 is name
    r=[]
    with open(file, "r", encoding="UTF-8") as read:
        tsv_file = csv.reader(read, delimiter = '\t')
        for line in tsv_file:
            out = [line[x] for x in i]      ##[ftype, EC, product]
            if("hypothetical protein" not in out[2]):
                r.append(out)
        r.pop(0)        ##removes headers
    return r

## string => string
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

## string => string
def extract_function(str):
    index = str.index("product=")+8
    present = ";" in str[index:]
    end = str[index:].index(";") if present else -1
    return str[index:index+end].replace("%2", ",") if end != -1 else str[index:].replace("%2", ",")

## File, string => string
## creates new tsv file with file[ftype, EC, product]
def protein_to_tsv(file, program):
    oname = program + "_" + file[:-3] + "tsv"
    with open (oname, 'w', newline ='') as csvfile:
        w = csv.writer(csvfile, delimiter = '\t')
        rows = unity_protein(file) if program == "unity" else kbase_protein(file)
        rows.insert(0, ["Type", "EC", "product"])
        w.writerows(rows)
    return oname


##Scraping methods

## [ftype, EC, product] => [ftype, EC, product, pathways]
## finds the correct html box with the pathways
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

## file[ftype, EC, product] => [ftype, EC, product, pathways]
## Runs the scrapping process on all products in the file
def ec_search(file):
    with open(file, 'r', encoding="UTF-8") as read:
        threads = 15
        files_to_search = []
        futures = []
        rr = []
        tsv_file = csv.reader(read, delimiter = "\t")
        for line in tsv_file:       ##culls only entries with EC numbers
            ec = line[1]
            if ec != None and ec != "" and ec != "EC" and "-" not in ec:
                files_to_search.append(line)
        update = 0

        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            for line in files_to_search:
                future = executor.submit(pathway, line)
                futures.append(future)
            

            for future in futures:
                rr.append(future.result())
                update += 1
                if(update % 20 == 0):
                    print(update)
            return rr
    
## file[ftype, EC, product] => file[ftype, EC, product, pathways]
## creates pathways tsv
def file_overwrite(file):
    overwrite = ec_search(file)
    overwrite.insert(0, ["Type", "EC", "Product", "Pathway"])
    with open(file[:-4]+"_pathways.tsv", "w", newline = '') as csvfile:
        write = csv.writer(csvfile, delimiter = '\t')
        write.writerows(overwrite)
    return file[:-4]+"_pathways.tsv"

## Produce count file

## file[ftype, EC, product, pathways] => file[pathway, count]
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

## string => string[]
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

## file => string
## Creates count tsv
def count_to_tsv(file):
    p = list(count(file).items())
    oname = file[:-12] + "counted.tsv"
    with open (oname, 'w', newline ='') as csvfile:
        w = csv.writer(csvfile, delimiter = '\t')
        w.writerow(["Pathway", "Count"])
        w.writerows(p)
    return oname

## Combined method

## file, string => None
## runs both the pathways and count tsvs.
def convert(file, program):
    name = protein_to_tsv(file, program)
    name2 = file_overwrite(name)
    count_to_tsv(name2)
    print("Deleting intermediate file " + name )
    os.remove(name)

## run it like convert("prokka_Rickettsiales.tsv", "unity")