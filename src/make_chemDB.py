import requests
import json
from pprint import pprint
from libchebipy import ChebiEntity



def getAllChEBIs(outFile):

        # fetch all ChEBI IDs from Rhea
    response = requests.get("https://www.rhea-db.org/rhea?/query=&columns=chebi-id&format=tsv")
    data = response.text

    data = data.split("\n")[1:-1]

    out = []
    for i in data:
        temp = i.split(";")
        for j in temp:
            out.append(j)
        
        # filter only unique IDs
    out = list(set(out))
    out = [{"CHEBI":i} for i in out]
    out = json.dumps(out)

        # cache file
    with open(outFile, "w+") as f:
        f.write(out)
        print("Chemical IDs from Rhea cached.")




def add_NameAndSMILES(chemicals):
    with open(chemicals, "r+") as f:
        data = json.load(f)

            # append name and SMILES for each entry
        for i in range(0, len(data)):
            data[i]["name"] = ChebiEntity(data[i]["CHEBI"]).get_name()
            data[i]["smiles"] = ChebiEntity(data[i]["CHEBI"]).get_smiles()
            data[i]["inchi_key"] = ChebiEntity(data[i]["CHEBI"]).get_inchi_key()
            
            # filter out all entries with "smiles" = None
        data = [i for i in data if i["smiles"] != None]

            # save updated file
        f.seek(0)
        f.write(json.dumps(data))
        f.truncate()
        print("added chemical name and smiles info for "+str(len(data))+" entries")



def make_chemDB(out):
    getAllChEBIs(out)
    add_NameAndSMILES(out)


if __name__ == "__main__":

        # file where I'm storing info on all chemicals included in Rhea
    chemicals = "data/all_rhea_chemicals.json"

    make_chemDB(chemicals)


