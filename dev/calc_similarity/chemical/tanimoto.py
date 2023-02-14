from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import json
import requests


def tanimoto_calc(smi1, smi2):
    try:
        mol1 = Chem.MolFromSmiles(smi1)
        mol2 = Chem.MolFromSmiles(smi2)
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)
        s = round(DataStructs.TanimotoSimilarity(fp1,fp2),3)
        return s
    except:
        return 0



def search_ligand_database(db: str, target: str):

    with open(db, "r") as f:
        data = json.load(f)

    def get_score(chem):
        return chem.get('score')

    scores = []
    for i in data:
        tani = tanimoto_calc(i["smiles"], target)
        scores.append({"name": i["name"], "score": tani, "smiles": i["smiles"]})
    
    scores.sort(key=get_score, reverse=True)
    print(scores[0:10])



    # Fetch data from GroovDB via an AWS lambda function.
        # I'll need to create a new lambda function to get the data we need.
# def get_smiles():
    
#     response = requests.get('https://4lsuwlkqoe.execute-api.us-east-2.amazonaws.com/search',
#         headers={ 'Accept': 'application/json', 'Authorization': 'groovDBAllowPostToDB_%FYaeF36!8'})

#     data = response.text
#     print(data)



if __name__ == "__main__":

        # isoprene
    target = "CC(=C)C=C"

        # databases
    groov = "groovLigands.json"
    rhea = "../../src/data/all_rhea_chemicals.json"

    data = search_ligand_database(rhea, target)
    print(data)