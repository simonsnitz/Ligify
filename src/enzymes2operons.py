import json
from src.accID2operon import acc2operon
import re
from pprint import pprint

def append_operons(chemical_name: str):

    with open("archives/"+chemical_name+".json", "r+") as data:
        ligand = json.load(data)

        if len(ligand["rxn_data"]) == 0:
            print("No enzymes found for "+str(chemical_name))
            return

        if "context" in ligand["rxn_data"][0]["proteins"][0]:
            print("operon data for "+chemical_name+" already cached")
            pass

        else:
            for rxn in range(0,len(ligand["rxn_data"])):
                for i in range(0, len(ligand["rxn_data"][rxn]["proteins"])):
                    protein = ligand["rxn_data"][rxn]["proteins"][i]
                    ncbi_id = protein["enzyme"]["ncbi_id"]
                    if ncbi_id != None:
                        context = acc2operon(ncbi_id)
                        protein["context"] = context
                        print("fetched context")
            
            data.seek(0)
            data.write(json.dumps(ligand))
            data.truncate()
            print("operon data for "+chemical_name+" cached in archives")




def pull_regulators(chemical_name: str):

    with open("archives/"+chemical_name+".json", "r") as data:
        regulator = re.compile(r"regulator|repressor|activator")
        ligand = json.load(data)

        for rxn in ligand["rxn_data"]:
            for protein in rxn["proteins"]:

                if "context" in protein.keys():
                    if protein["context"] != "EMPTY":
                        operon = protein["context"]["operon"]
                        for gene in operon:
                            if "description" in gene.keys():
                                if regulator.search(gene["description"]):
                                    print(gene["accession"])
                                    #print(protein["organism"])
