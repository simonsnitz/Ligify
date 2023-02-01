import json
from src.accID2operon import acc2operon
import re
import requests
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


    






def pull_regulators(chemical_name: str):

    with open("archives/"+chemical_name+".json", "r") as data:
        regulator = re.compile(r"regulator|repressor|activator")
        ligand = json.load(data)

        for rxn in ligand["rxn_data"]:
            for protein in rxn["proteins"]:
                ligand_names = []
                if "context" in protein.keys():
                    if protein["context"] != "EMPTY":
                        operon = protein["context"]["operon"]
                        for gene in operon:
                            if "description" in gene.keys():
                                if regulator.search(gene["description"]):
                                    print(gene["accession"])
                                    for gene in operon:
                                        protein_data = protein2chemicals(gene["accession"])
                                        if isinstance(protein_data, dict):
                                            if "catalysis" in protein_data.keys():
                                                ligand_names.append(protein_data["catalysis"].split(" "))
                                    print("\n")
                print(ligand_names)
