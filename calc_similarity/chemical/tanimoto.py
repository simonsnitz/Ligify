from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import json


def tanimoto_calc(smi1, smi2):
    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)
    s = round(DataStructs.TanimotoSimilarity(fp1,fp2),3)
    return s

    # isoprene
target = "CC(=C)C=C"

with open("groovLigands.json", "r") as f:
    data = json.load(f)

scores = []
for i in data:
    tani = tanimoto_calc(i["smiles"], target)
    scores.append({"name": i["name"], "score": tani})

def get_score(chem):
    return chem.get('score')

scores.sort(key=get_score, reverse=True)
print(scores)

# print(tanimoto_calc(berberine, thp))