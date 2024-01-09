import requests
import json
import re
from pprint import pprint
import time
from math import ceil
import sys

from libchebipy import ChebiEntity

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs




def create_database():



    def create_chemical_map():

        # Fetch all CHEBI IDs in RHEA
        url= "https://www.rhea-db.org/rhea?"
        parameter = {
        "query":'uniprot:*',
        "columns":"chebi-id",
        "format":'tsv',
        }
        response = requests.get(url,params=parameter)
        if response.ok:
            chemDB = []
            e = response.text.split("\n")[1:-1]
            for chem in e:
                chemicals = chem.split(";")
                for chem in chemicals:
                    chemDB.append(chem)

        # Filter out redundant CHEBI IDs
        unique = []
        for chem in chemDB:
            if chem not in unique:
                unique.append(chem)
        rhea_chebis = unique

        # Map CHEBI IDs to chemical name and smiles
        db = []
        for chebi in rhea_chebis:
            smiles = ChebiEntity(chebi).get_smiles()
            name = ChebiEntity(chebi).get_name()
            entry = {
                "chebi": chebi,
                "name": name,
                "smiles": smiles,
            }
            db.append(entry)
        with open("chemical_map.json", "w+") as out:
            out.write(json.dumps(db))
        print("Created the chemical map")
        

    create_chemical_map()

    




def blast_chemical(smiles, max_alt_chems):

    def tanimoto_calc(smi1, smi2):
        mol1 = Chem.MolFromSmiles(smi1)
        mol2 = Chem.MolFromSmiles(smi2)
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)
        s = round(DataStructs.TanimotoSimilarity(fp1,fp2),3)
        return s

    with open("ligify/predict/database/chemical_map.json", "r") as data:
    #with open("chemical_map.json", "r") as data:
    #with open("all_chemicals.json", "r") as data:
        db = json.load(data)
        
    scores = []

    for i in db:
        try:
            score = tanimoto_calc(smiles, i["smiles"])
            scores.append({"Similarity":score, "Name": i["name"], "SMILES": i["smiles"]})
        except: 
            pass

    scores = sorted(scores, key=lambda x: x["Similarity"], reverse=True)[0:max_alt_chems]

    return scores





if __name__ == "__main__":

    #create_database()

    pprint(blast_chemical("OC[C@H]1OC(O)[C@H](O)[C@@H]1O", 10))