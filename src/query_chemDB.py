from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import json


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


def query_chemDB(db: str, target: str):

    with open(db, "r") as f:
        data = json.load(f)

    def get_score(chem):
        return chem.get('score')

    scores = [ 
            {"name": i["name"],     
            "score": tanimoto_calc(i["smiles"], target),    
            "smiles": i["smiles"], 
            "inchi_key": i["inchi_key"]
            } for i in data]
    
    scores.sort(key=get_score, reverse=True)

    return scores



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
    rhea = "data/all_rhea_chemicals.json"

    data = query_chemDB(rhea, target)
    print(data[0:5])