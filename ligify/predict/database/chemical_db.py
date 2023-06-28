import requests
import json
import re
from pprint import pprint
import time
from math import ceil
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs




def create_database():

# Create a json file with all rhea reactions and associated chemical IDs for each

    url= "https://www.rhea-db.org/rhea?"
    parameter = {
    "query":'uniprot:*',
    "columns":"rhea-id,chebi-id",
    "format":'tsv',
    }
    response = requests.get(url,params=parameter)
    if response.ok:
        rheaDB = {}
        e = response.text.split("\n")[1:-1]
        for rhea in e:
            rhea = rhea.split("\t")
            rhea_id = rhea[0]
            rheaDB[rhea_id] = []
            chemicals = rhea[1].split(";")
            for chem in chemicals:
                rheaDB[rhea[0]].append(chem.strip("CHEBI:"))

    print("rhea IDs collected")

# Format the chemical database to link from rhea -> chem ID, to chem ID -> rhea

    chemDB = {}

    for rhea in rheaDB:
        for chem in rheaDB[rhea]:
            if chem not in chemDB:
                chemDB[chem] = [rhea]
            elif chem in chemDB:
                chemDB[chem].append(rhea)

# filter chemical database to drug-like molecules that likely can cross the cell membrane

    num_batches = ceil(len(chemDB)/1000)
    chem_json = []
    
    for i in range(0, num_batches):
        
        # extract the group of 1000 relevant chemical IDs
        counter = 0
        start = i*1000
        end = (i+1)*1000
        chebis = ""
        for chebi in chemDB:
            if counter >= start and counter < end:
                chebis += "CHEBI:"+str(chebi)+","
            counter += 1

        chebis = chebis[0:-1]

        # get chemical data from an API request to mychem
        params = {'ids':chebis, 'fields':\
            'pubchem.smiles.isomeric,\
            pubchem.smiles.canonical,\
            chebi.name'}

        res = requests.post('http://mychem.info/v1/chem', params)
        con = res.json()
        
        a = 0
        f = 0

        for chem in con:

            chebid = chem['query'].strip("CHEBI:")
            try:
                name = chem["chebi"]["name"]
            except:
                try:
                    name = chem["chebi"][0]["name"]
                except:
                    name = "unknown"

            try:
                smiles = chem["pubchem"]["smiles"]["isomeric"]
            except:
                try:
                    smiles = chem["pubchem"]["smiles"]["canonical"]
                except:
                    f += 1
                    smiles = None

            if smiles != None:
                # add to new dictionary
                chem_dict = { 
                    "name": name,
                    "smiles": smiles }
                chem_json.append(chem_dict)
                a += 1

        print("added: "+ str(a))
        print("no smiles: "+ str(f))


    with open("all_chemicals.json", "w+") as o:
        out = json.dumps(chem_json)
        o.write(out)
        print('created RHEA chemical database')




def blast_chemical(smiles, max_alt_chems):

    def tanimoto_calc(smi1, smi2):
        mol1 = Chem.MolFromSmiles(smi1)
        mol2 = Chem.MolFromSmiles(smi2)
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)
        s = round(DataStructs.TanimotoSimilarity(fp1,fp2),3)
        return s

    with open("ligify/predict/database/all_chemicals.json", "r") as data:
        db = json.load(data)
        
    scores = []

    for i in db:
        score = tanimoto_calc(smiles, i["smiles"])
        if score != 1.0:
            scores.append({"Similarity":score, "Name": i["name"], "SMILES": i["smiles"]})

    scores = sorted(scores, key=lambda x: x["Similarity"], reverse=True)[0:max_alt_chems]

    return scores





if __name__ == "__main__":

    create_database()

    #blast_chem("C1=CC(=C(C=C1CCN)O)O")