import requests
import json
from pprint import pprint
import time
# from pubchem import get_inchiKey


# TODO:
# Get Rhea IDs and equations from PubChem. 


def initialize_database():
    with open("all_chemicals.json", "w+") as f:

        url= "https://www.rhea-db.org/rhea?"
        parameter = {
        "query":'uniprot:*',
        "columns":"rhea-id,chebi-id",
        "format":'tsv',
        }
        response = requests.get(url,params=parameter)
        if response.ok:
            data = {}
            e = response.text.split("\n")[1:-1]
            for rhea in e:
                rhea = rhea.split("\t")
                rhea_id = rhea[0]
                data[rhea_id] = []
                chemicals = rhea[1].split(";")
                for chem in chemicals:
                    data[rhea[0]].append(chem.strip("CHEBI:"))
        
        f.write(json.dumps(data))





def format_data():
    with open("all_chemicals.json", "r+") as f:

        data = json.load(f)

        chemDB = {}

        for rhea in data:
            for chem in data[rhea]:
                if chem not in chemDB:
                    chemDB[chem] = [rhea]
                elif chem in chemDB:
                    chemDB[chem].append(rhea)

        with open("formatted.json", "w+") as out:
            out.write(json.dumps(chemDB))







def append_genes():

    url = "https://rest.uniprot.org/uniprotkb/search?format=json&size=25&query=reviewed:true+AND+"


    with open("formatted.json", "r") as f:
            
        db = json.load(f)

        with open("all_proteins.json", "w+") as out:

            db_with_proteins = {}

            counter = 0
            for chemical in db:
                if counter < 10:
                    rhea_ids = [i.strip("RHEA:") for i in db[chemical]]

                    rhea_suffix = "".join("rhea:"+str(id)+"+OR+" for id in rhea_ids)[:-4]

                    # Remove chemicals associated with A BUNCH of reactions.
                        # These are likely co-factors or non-descript chemicals
                    if len(rhea_ids) < 10:

                        response = requests.get(url+rhea_suffix)
                        data = json.loads(response.text)["results"]   

                        proteins = []
                        
                        for entry in data:
                            if entry["organism"]["lineage"][0] == "Bacteria":

                                # Get reference DOIs
                                description = entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
                                dois = []
                                for j in entry['references']:
                                    if "citationCrossReferences" in j["citation"]:
                                        for k in j["citation"]["citationCrossReferences"]:
                                            if k["database"] == "DOI":
                                                dois.append(k["id"])

                                # Get RefSeq ID
                                # The "NCBI_ID" is needed to get the genome context in the next step. Prefer to use RefSeq, but can use EMBL.
                                try:
                                    ncbi_id = [e["id"] for e in entry["uniProtKBCrossReferences"] if e["database"] == "RefSeq"][0]
                                except:
                                    try:
                                        ncbi_id = [e["properties"] for e in entry["uniProtKBCrossReferences"] if e["database"] == "EMBL"][0]
                                        ncbi_id = [e["value"] for e in ncbi_id if e["key"] == "ProteinId"][0]
                                    except:
                                        ncbi_id = None

                                # Get Uniprot ID and Organism
                                uniprotID = entry["primaryAccession"]
                                organism = entry["organism"]["lineage"]

                                # Format protein data into a dictionary
                                protein = {
                                    "organism": organism,
                                    "enzyme": {
                                        "description": description,
                                        "uniprot_id": uniprotID,
                                        "dois": dois,
                                        "ncbi_id": ncbi_id,
                                    }
                                }
                                if ncbi_id != None:
                                    proteins.append(protein)
                

                        counter += 1
                        if len(proteins) != 0:
                            db_with_proteins[chemical] = proteins
                            print("finished "+str(counter)+" of "+(str(len(db))))
                            out.write(json.dumps(db_with_proteins))
                    
        #print(counter)




def filter_genes(output, lineage_filter_name):

    rxns = output["rxn_data"]
        # filter out empties
    rxns = [i for i in rxns if len(i["proteins"]) != 0]

        #DOUBLE LIST COMPREHENSION!!!
    # num_proteins = len([protein for rxn in rxns for protein in rxn["proteins"]])
    # output["metadata"]["number_reviewed_enzymes"] = num_proteins


    #     # have to map name to number because names are more human readable, but numbers are how
    #         # lineages are retrieved programmatically via the Uniprot API.
    map_lineage2number = {"Domain": 0, "Phylum": 1, "Class": 2, "Order": 3, "Family": 4}
    lineage_filter = map_lineage2number[lineage_filter_name]


    # filter out highly similar proteins
    filtered_rxns = []
    for rxn in rxns:
        filtered_proteins = []
        families = []
        for protein in rxn["proteins"]:
                # Sometimes the family name isn't provided in the lineage returned.
            try:
                family = protein["organism"][lineage_filter]
            except:
                family = ""
            if family not in families:
                families.append(family)
                filtered_proteins.append(protein)
        new_rxn = rxn
        new_rxn["proteins"] = filtered_proteins
        filtered_rxns.append(new_rxn)

    output["rxn_data"] = filtered_rxns

        # count number of filtered proteins
    # filtered_proteins = len([protein for rxn in filtered_rxns for protein in rxn["proteins"]])
    # output["metadata"]["number_lineage_filtered_enzymes"] = filtered_proteins


    return output


    # with open("archives/"+chemical_name+".json", "w+") as f:
    #     f.write(json.dumps(output))
    #     print(str(filtered_proteins)+" enzymes for "+chemical_name+" cached in archives")
    


if __name__ == "__main__":

    # pprint(fetch_reactions("C=CC(=O)[O-]"))
    # initialize_database()
    # format_data()
    append_genes()


    # with open("formatted.json", "r") as f:

    #     data = json.load(f)
    #     print(data["4986"])

    # with open("all_chemicals.json", "r") as f:

    #     data = json.load(f)

    #     for rhea in data:
    #         for chem in data[rhea]:
    #             if chem == "4986":
    #                 print(rhea)
        # print(data["4986"])




