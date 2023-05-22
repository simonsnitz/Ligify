import json
import re
import requests
import pandas as pd
from pprint import pprint

from ligify.predict.accID2operon import acc2operon
#from accID2operon import acc2operon

def append_operons(data, chemical_name):

    ligand = data


    if len(ligand["rxn_data"]) == 0:
        print("No enzymes found for "+str(chemical_name))
        return

    if "context" in ligand["rxn_data"][0]["proteins"][0]:
        print("operon data for "+chemical_name+" already cached")
        pass

    else:
        counter = 0
        for rxn in range(0,len(ligand["rxn_data"])):
            for i in range(0, len(ligand["rxn_data"][rxn]["proteins"])):
                protein = ligand["rxn_data"][rxn]["proteins"][i]
                ncbi_id = protein["enzyme"]["ncbi_id"]
                if ncbi_id != None:
                    if counter <= 25:
                        context = acc2operon(ncbi_id)
                        protein["context"] = context
                        print("fetched context")
                        counter += 1
        
        return ligand


#TODO:
# Get all chemicals associated with the regulator-bearing operons (append a chemicals:[])
# Get all literature associated with each protein in the operon (append a doi:[])




def protein2chemicals(accession: str):
    url = "https://rest.uniprot.org/uniprotkb/search?query="+accession+"&format=json"

    response = requests.get(url)
    if response.ok:
        protein = json.loads(response.text)

        ligand_names = []
        ligand_ids = []

        if len(protein['results']) > 0:
            if "comments" in protein["results"][0]:
                protein = protein["results"][0]["comments"]

                protein_data = {}

                    # look for description
                Function = [i["texts"][0]["value"] for i in protein if i["commentType"] == "FUNCTION"]
                if len(Function) != 0:
                    protein_data["function"] = Function[0]
        

                    # look for catalytic activity
                CATALYSIS = [i["reaction"]["name"] for i in protein if i["commentType"] == "CATALYTIC ACTIVITY"]
                if len(CATALYSIS) != 0:
                    protein_data["catalysis"] = CATALYSIS[0]
                    

                    # look for catalytic activity
                try:
                    RXN = [i["reaction"]["reactionCrossReferences"] for i in protein if i["commentType"] == "CATALYTIC ACTIVITY"]
                    if len(RXN) != 0:
                        LIGANDS = [i["id"] for i in RXN[0] if i["database"] == "ChEBI"]
                        if len(LIGANDS) != 0:
                            protein_data["ligands"] = LIGANDS
                except:
                    pass

                    # look for induction
                INDUCTION = [i["texts"][0]["value"] for i in protein if i["commentType"] == "INDUCTION"]
                if len(INDUCTION) != 0:
                    protein_data["induction"] = INDUCTION[0]

                            # look for induction
                PATHWAY = [i["texts"][0]["value"] for i in protein if i["commentType"] == "PATHWAY"]
                if len(PATHWAY) != 0:
                    protein_data["pathway"] = PATHWAY[0]


                # add something to append all this metadata


                
                # try:
                #     pprint(protein_data["catalysis"])
                # except:
                #     pass
                return(protein_data)
    else:
        response.raise_for_status()


    






def pull_regulators(data, chemical_name):


    regulator = re.compile(r"regulator|repressor|activator")

    reg_data = []

    for rxn in data["rxn_data"]:
        for protein in rxn["proteins"]:
            ligand_names = []
            if "context" in protein.keys():
                if protein["context"] != "EMPTY":
                    operon = protein["context"]["operon"]
                    for gene in operon:
                        if "description" in gene.keys():
                            if regulator.search(gene["description"]):

                                entry = {   "refseq": gene["accession"],
                                            "annotation": gene["description"],
                                            "protein": protein,
                                            "equation": rxn["equation"],
                                            "rhea_id": rxn["rhea_id"],
                                            }

                                for gene in operon:
                                    protein_data = protein2chemicals(gene["accession"])
                                    if isinstance(protein_data, dict):
                                        if "catalysis" in protein_data.keys():
                                            ligand_names += protein_data["catalysis"].split(" ")
                                unique_ligands = list(set(ligand_names))
                                not_ligands = ["H2O", "+", "-", "=", "A", "AH2", "H(+)", "NADPH", "NADH", "NADP(+)", "NAD(+)", str(chemical_name).lower()]
                                unique_ligands = [ i for i in unique_ligands if i not in not_ligands]

                                entry['alt_ligands'] = unique_ligands

                                reg_data.append(entry)

                
    return reg_data


if __name__ == "__main__":

    with open("temp/all.json", "r") as f:
        all_chemicals = json.load(f)

        for chemical in all_chemicals:
            print(chemical)
            # data = json.load(f)

            # regs = pull_regulators(data)

            # #print(regs)
            # d = pd.DataFrame(regs)
            # print(d)