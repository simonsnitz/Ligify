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
from compare import blast_align


def fetch_reg_protein_seq(accession: str):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accession+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        seq = "".join(i for i in response.text.split("\n")[1:])
        return seq
    else:
        print("Bad eFetch request "+ str(response.status_code))
        return None




def benchmark(file_name, rigerous=False):

    with open(file_name, "rb+") as f:

        df = pd.read_excel(f)

        refseqs = df.loc[:,"refseq"].values

        # Get protein sequences, if they don't already exist.
        for i in range(0, len(refseqs)):
            
            if pd.isna(df.loc[i,"protein_seq"]):
                protein_seq = fetch_reg_protein_seq(df.loc[i, "refseq"])
                df.loc[i, "protein_seq"] = protein_seq
                print("fetched protein seq for "+str(df.loc[i, "refseq"]))

        df.to_excel(file_name)

        # Collect ligify predictions for each ligand, if they don't alredy exist.
        smiles = df.loc[:,"smiles"].values
        for i in range(0, len(smiles)):

            if pd.isna(df.loc[i, "ligify_predictions"]):
                genus = df.loc[i, "Genus"]
                prediction = fetch_data(df.loc[i, "smiles"], rigerous, genus)
                prediction = json.dumps(prediction)
                #print(prediction)
                df.loc[i, "ligify_predictions"] = prediction
                df.to_excel(file_name)


                if pd.isna(df.loc[i, "match"]):
                                true_reg_seq = df.loc[i, "protein_seq"]
                                data = prediction
                                #data = df.loc[i, "ligify_predictions"]

                                reg_list = json.loads(data)["reg_list"]
                                alignments = []
                                if reg_list != None:
                                    for reg in reg_list:
                                        predicted_reg_seq = reg["protein_seq"]
                                        if len(predicted_reg_seq) > 0:
                                            align = blast_align(predicted_reg_seq, true_reg_seq)
                                            align["refseq"] = reg["refseq"]
                                            align["rank"] = reg["rank"]["rank"]
                                            alignments.append(align)
                                    # pull out the best match
                                    min_e_score = min([i["e value"] for i in alignments])
                                    top_prediction = [i for i in alignments if i["e value"] == min_e_score][0]
                                    df.loc[i, "align_details"] = json.dumps(top_prediction)
                                    # assign pass/fail
                                    if top_prediction["identity"] >80 and top_prediction["coverage"] >90:
                                        df.loc[i, "match"] = True
                                        print("Matches!")
                                    else:
                                        df.loc[i, "match"] = False
                                        print("Does not match :(")
                                    print("appended comparison for "+str(refseqs[i]))


                df.to_excel(file_name)
                print("appended a prediction for "+str(smiles[i]))

        # Find the best true regulator-prediction match.
        # for i in range(0, len(refseqs)):

        #     if pd.isna(df.loc[i, "match"]):
        #         true_reg_seq = df.loc[i, "protein_seq"]
        #         data = df.loc[i, "ligify_predictions"]

        #         reg_list = json.loads(data)["reg_list"]
        #         alignments = []
        #         if reg_list != None:
        #             for reg in reg_list:
        #                 predicted_reg_seq = reg["protein_seq"]
        #                 if len(predicted_reg_seq) > 0:
        #                     align = blast_align(predicted_reg_seq, true_reg_seq)
        #                     align["refseq"] = reg["refseq"]
        #                     align["rank"] = reg["rank"]["rank"]
        #                     alignments.append(align)
        #             # pull out the best match
        #             min_e_score = min([i["e value"] for i in alignments])
        #             top_prediction = [i for i in alignments if i["e value"] == min_e_score][0]
        #             df.loc[i, "align_details"] = json.dumps(top_prediction)
        #             # assign pass/fail
        #             if top_prediction["identity"] >80 and top_prediction["coverage"] >90:
        #                 df.loc[i, "match"] = True
        #             else:
        #                 df.loc[i, "match"] = False
        #             print("appended comparison for "+str(refseqs[i]))

        #             df.to_excel(file_name)



def filter_by_genus(output, genus):

    rxns = output["rxn_data"]

        #DOUBLE LIST COMPREHENSION!!!
    # num_proteins = len([protein for rxn in rxns for protein in rxn["proteins"]])
    # output["metadata"]["number_reviewed_enzymes"] = num_proteins


    #     # have to map name to number because names are more human readable, but numbers are how
    #         # lineages are retrieved programmatically via the Uniprot API.


    # filter out highly similar proteins
    filtered_rxns = []
    for rxn in rxns:
        filtered_proteins = []
        for protein in rxn["proteins"]:
                # Sometimes the family name isn't provided in the lineage returned.
            try:
                if protein["organism"][5] == genus:
                    filtered_proteins.append(protein)
            except:
                pass
        new_rxn = rxn
        new_rxn["proteins"] = filtered_proteins
        filtered_rxns.append(new_rxn)

    output["rxn_data"] = filtered_rxns

        # count number of filtered proteins
    # filtered_proteins = len([protein for rxn in filtered_rxns for protein in rxn["proteins"]])
    # output["metadata"]["number_lineage_filtered_enzymes"] = filtered_proteins


    return output
                




def fetch_data(smiles, rigerous, genus):

    InChiKey = get_inchikey(smiles, "smiles")

    if rigerous == 1:
        max_react = 100
        max_genes = 100
        reviewed = False
        lineage_filter = "None"
        max_operons = 100
    elif rigerous == 2:
        max_react = 100
        max_genes = 500
        reviewed = True
        lineage_filter = "Family"
        max_operons = 500
    elif rigerous == 3:
        max_react = 100
        max_genes = 500
        reviewed = True
        lineage_filter = "None"
        max_operons = 1000

    else:
        max_react = 100
        max_genes = 500
        reviewed = True
        lineage_filter = "Family"
        max_operons = 100

    metrics = {}

    # FETCH REACTIONS

    print("fetching reactions for "+str(smiles))

    reactions = fetch_reactions(InChiKey = InChiKey, max_reactions = max_react)
    total_rxns = len(reactions["rxn_data"])
    metrics["RHEA Reactions"] = total_rxns

    print("fetched "+str(total_rxns)+" total reactions")

    if total_rxns > 0:

        # FETCH GENES

        print("fetching genes for "+str(smiles))

        for i in reactions["rxn_data"]:
            associated_proteins = fetch_genes(i["rhea_id"], reviewed, max_genes)
            i["proteins"] = associated_proteins

        total_genes = sum([len(i["proteins"]) for i in reactions["rxn_data"]])
        metrics["Total genes"] = sum([len(i["proteins"]) for i in reactions["rxn_data"]])

        print("fetched "+str(total_genes)+" total genes")

        # Filter homologous genes
        reactions = filter_genes(reactions, lineage_filter_name = lineage_filter)
        total_filtered_genes = sum([len(i["proteins"]) for i in reactions["rxn_data"]])
        metrics["Filtered genes"] = sum([len(i["proteins"]) for i in reactions["rxn_data"]])

        print(str(total_filtered_genes)+" total filtered genes") 

        # Filter genes by genus
        reactions = filter_by_genus(reactions, genus)
        total_filtered_genes = sum([len(i["proteins"]) for i in reactions["rxn_data"]])

        print(str(total_filtered_genes)+" total genes after genus filter") 

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
            print("Total genes: "+str(total_genes))

            for rxn in range(0,len(reactions["rxn_data"])):
                for i in range(0, len(reactions["rxn_data"][rxn]["proteins"])):
                    protein = reactions["rxn_data"][rxn]["proteins"][i]
                    refseq_id = protein["enzyme"]["ncbi_id"]
                    if refseq_id != None:

                        # Limit number of operons evaluated to avoid program taking too long to complete.
                        if operon_counter <= max_operons:
                            protein["context"] = acc2operon(refseq_id)
                            print("appended operon "+str(operon_counter))
                            operon_counter += 1

            metrics["Total operons"] = operon_counter
            print("Total operons: "+str(operon_counter))



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
                print("Total regulators: "+str(len(regulators)))

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

    #file_name = "Second_pass_Ligify_benchmarking_dataset.xlsx"
    # file_name = "Rerun2_of_Incorrectly-predicted_Ligify_benchmarking_dataset.xlsx"
    #file_name = "Remaining_July12_Benchmarking_dataset.xlsx"
    file_name = "Aug6_Benchmarking_dataset.xlsx"

    benchmark(file_name, rigerous=3)

    #data = fetch_data("C=CC(=O)[O-]")
    #print(data)