from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import json
import requests


def tanimoto_calc(smi1, smi2):
    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)
    s = round(DataStructs.TanimotoSimilarity(fp1,fp2),3)
    return s


def compare_chemicals(target: str):
        # isoprene
    # target = "CC(=C)C=C"
    with open("groovLigands.json", "r") as f:
        data = json.load(f)

    def get_score(chem):
        return chem.get('score')

    scores = []
    for i in data:
        tani = tanimoto_calc(i["smiles"], target)
        scores.append({"name": i["name"], "score": tani})
    
    scores.sort(key=get_score, reverse=True)
    print(scores)



    # Fetch data from GroovDB via an AWS lambda function.
        # I'll need to create a new lambda function to get the data we need.
def get_smiles():
    
    response = requests.get('https://4lsuwlkqoe.execute-api.us-east-2.amazonaws.com/search',
        headers={ 'Accept': 'application/json', 'Authorization': 'groovDBAllowPostToDB_%FYaeF36!8'})

    data = response.text
    print(data)




data = get_smiles()
print(data)