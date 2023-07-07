import pandas as pd
import requests
import sys
import json
sys.path.append("..") # Adds higher directory to python modules path.

from chemical2enzymes import fetch_reactions, fetch_genes, filter_genes
from enzymes2operons import pull_regulators
from accID2operon import acc2operon
from rank import calculate_rank
from pubchem import get_inchikey


def fetch_reg_protein_seq(accession: str):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accession+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        seq = "".join(i for i in response.text.split("\n")[1:])
        return seq
    else:
        print("Bad eFetch request "+ str(response.status_code))
        return None




def benchmark():

    with open("Ligify_benchmarking_dataset.xlsx", "rb+") as f:

        df = pd.read_excel(f)

        refseqs = df.loc[:,"refseq"].values

        # if protein sequence doesn't exist, get it.
        for i in range(0, len(refseqs)):
            
            if pd.isna(df.loc[i,"protein_seq"]):
                protein_seq = fetch_reg_protein_seq(df.loc[i, "refseq"])
                df.loc[i, "protein_seq"] = protein_seq
                print("fetched protein seq for "+str(df.loc[i, "refseq"]))

        df.to_excel("Ligify_benchmarking_dataset.xlsx")

        # Collect ligify predictions for each ligand
        smiles = df.loc[:,"smiles"].values
        for i in range(0, len(smiles)):

            if pd.isna(df.loc[i, "ligify_predictions"]):
                prediction = fetch_data(df.loc[i, "smiles"])
                prediction = json.dumps(prediction)
                print(prediction)
                df.loc[i, "ligify_predictions"] = prediction
                df.to_excel("Ligify_benchmarking_dataset.xlsx")
                print("appended a prediction for "+str(smiles[i]))





def fetch_data(smiles):

    InChiKey = get_inchikey(smiles, "smiles")

    metrics = {}

    # FETCH REACTIONS

    print("fetching reactions for "+str(smiles))

    reactions = fetch_reactions(InChiKey = InChiKey, max_reactions = 20)
    total_rxns = len(reactions["rxn_data"])
    metrics["RHEA Reactions"] = total_rxns

    if total_rxns > 0:

        # FETCH GENES

        print("fetching genes for "+str(smiles))

        for i in reactions["rxn_data"]:
            associated_proteins = fetch_genes(i["rhea_id"], True, 20)
            i["proteins"] = associated_proteins

        metrics["Total genes"] = sum([len(i["proteins"]) for i in reactions["rxn_data"]])

        # Filter homologous genes
        reactions = filter_genes(reactions, lineage_filter_name = "Family")
        metrics["Filtered genes"] = sum([len(i["proteins"]) for i in reactions["rxn_data"]])

        # FETCH OPERONS

        print("fetching operons for "+str(smiles))

        if len(reactions["rxn_data"]) == 0:
            print("No enzymes found")
            return {"reg_list": None, "metrics": metrics}
        else:
            operon_counter = 0

            # Get the total number of valid protein entries for operon fetching
            # This can be shortened using list comprehension
            total_genes = 0
            for rxn in range(0,len(reactions["rxn_data"])):
                for i in range(0, len(reactions["rxn_data"][rxn]["proteins"])):
                    protein = reactions["rxn_data"][rxn]["proteins"][i]
                    if protein["enzyme"]["ncbi_id"] !=  None:
                        total_genes +=1

            for rxn in range(0,len(reactions["rxn_data"])):
                for i in range(0, len(reactions["rxn_data"][rxn]["proteins"])):
                    protein = reactions["rxn_data"][rxn]["proteins"][i]
                    refseq_id = protein["enzyme"]["ncbi_id"]
                    if refseq_id != None:

                        # Limit number of operons evaluated to avoid program taking too long to complete.
                        if operon_counter <= 20:
                            protein["context"] = acc2operon(refseq_id)
                            operon_counter += 1

            metrics["Total operons"] = operon_counter


            # FETCH REGULATORS

            print("fetching regulators for "+str(smiles))

            if reactions == None:
                return {"reg_list": None, "metrics": metrics}
            
            else:

                # This is where all of the display data is created
                regulators = []
                for rxn in reactions["rxn_data"]:
                    for protein in rxn["proteins"]:
                        regs = pull_regulators(protein, rxn)
                        for r in regs:
                            regulators.append(r)

                metrics["Total regulators"] = len(regulators)

                # Filter out duplicate regulators
                refseq_ids = []
                filtered_regulators = []
                for i in regulators:
                    if i["refseq"] not in refseq_ids:
                        filtered_regulators.append(i)
                        refseq_ids.append(i["refseq"])


                # Create a rank for each regulator
                for r in filtered_regulators:
                    rank = calculate_rank(r)
                    r["rank"] = rank


                if filtered_regulators == None or len(filtered_regulators) == 0:
                    return {"reg_list": None, "metrics": metrics}
                else:

                    reg_list = []
                    for i in filtered_regulators:
                        reg = {}
                        reg["refseq"] = i["refseq"]
                        reg["protein_seq"] = i["reg_protein_seq"]
                        reg["rank"] = i["rank"]
                        reg_list.append(reg)

                    out = {"reg_list": reg_list, "metrics": metrics}

                    # Filter this to return only metrics and regulator refseq IDs.
                    return out
    else:
        return {"reg_list": None, "metrics": metrics}



if __name__ == "__main__":

    benchmark()

    #data = fetch_data("C=CC(=O)[O-]")
    #print(data)