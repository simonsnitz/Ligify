import json
import re
import requests



#TODO:
# Speed up protein2chemicals using Uniprot's SPARQL API?


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


    
def fetch_reg_data(accession: str):
    url = "https://rest.uniprot.org/uniprotkb/search?query="+accession+"&format=json"




def pull_regulators(protein, rxn):


    regulator = re.compile(r"regulator|repressor|activator")

    reg_data = []

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

                        ### This is where a Uniprot API query goes to fetch more info on the regulator.


                        # Fetch possible alternative inducer molecules associated with the operon
                        for gene in operon:
                            protein_data = protein2chemicals(gene["accession"])
                            if isinstance(protein_data, dict):
                                if "catalysis" in protein_data.keys():
                                    ligand_names += protein_data["catalysis"].split(" ")
                        unique_ligands = list(set(ligand_names))
                        not_ligands = ["H2O", "+", "-", "=", "A", "AH2", "H(+)", "NADPH", "NADH", "NADP(+)", "NAD(+)"]
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