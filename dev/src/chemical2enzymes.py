import requests
import json
from pprint import pprint
import os


#InChiKey = "KWIUHFFTVRNATP-UHFFFAOYSA-N"   # Found BetI
#InChiKey = "KBPLFHHGFOOTCA-UHFFFAOYSA-N"


InChiKey = "NIXOWILDQLNWCW-UHFFFAOYSA-M"    # Found AcuR
chemical_name = "acrylate"
domain_filter = "Bacteria"
lineage_filter_name = "Family"


def chem2enzymes(InChiKey: str, chemical_name: str, domain_filter: str, lineage_filter_name: str, reviewed: str):
    
    if os.path.isfile("archives/"+chemical_name+".json") == True:
        print("enzyme data already cached")
        pass

    else:
            # have to map name to number because names are more human readable, but numbers are how
                # lineages are retrieved programmatically via the Uniprot API.
        map_lineage2number = {"Domain": 0, "Phylum": 1, "Class": 2, "Order": 3, "Family": 4}
        lineage_filter = map_lineage2number[lineage_filter_name]


            # get rhea ids from chemical
        url= "https://www.rhea-db.org/rhea?"
        parameter = {
        "query":'InChiKey:'+str(InChiKey),
        "columns":"rhea-id",
        "format":'json',
        }
        response = requests.get(url,params=parameter)

        data = json.loads(response.text)["results"]

        output = {
            "metadata": {
                "chemical name": chemical_name,
                "InChiKey": InChiKey,
                "domain_filter": domain_filter,
                "lineage_similarity_filter": lineage_filter_name
            }
        }

        output["rxn_data"] = [{"rhea_id": i["id"], "equation": i["equation"]} for i in data]

        rxns = output["rxn_data"]
        
            # is there a limit of 25 entries here???
        url = "https://rest.uniprot.org/uniprotkb/search?format=json&query=reviewed:"+reviewed+"+AND+rhea:"


        for i in range(0,len(rxns)):

            response = requests.get(url+str(rxns[i]["rhea_id"]))
            data = json.loads(response.text)["results"]   

            proteins = []
            for entry in data:

                if entry["organism"]["lineage"][0] == "Bacteria":
                    description = entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
                    dois = []
                    for j in entry['references']:
                        if "citationCrossReferences" in j["citation"]:
                            for k in j["citation"]["citationCrossReferences"]:
                                if k["database"] == "DOI":
                                    dois.append(k["id"])
                        # The "NCBI_ID" is needed to get the genome context in the next step. Prefer to use RefSeq, but can use EMBL.
                    try:
                        ncbi_id = [e["id"] for e in entry["uniProtKBCrossReferences"] if e["database"] == "RefSeq"][0]
                    except:
                        try:
                            ncbi_id = [e["properties"] for e in entry["uniProtKBCrossReferences"] if e["database"] == "EMBL"][0]
                            ncbi_id = [e["value"] for e in ncbi_id if e["key"] == "ProteinId"][0]
                        except:
                            ncbi_id = None
                            print("no ncbi id retrived")

                    uniprotID = entry["primaryAccession"]
                    organism = entry["organism"]["lineage"]

                    protein = {
                        "organism": organism,
                        "enzyme": {
                            "description": description,
                            "uniprot_id": uniprotID,
                            "dois": dois,
                            "ncbi_id": ncbi_id,
                        }
                    }
                    proteins.append(protein)
            
            rxns[i]["proteins"] = proteins


            # filter out empties
        rxns = [i for i in rxns if len(i["proteins"]) != 0]

            #DOUBLE LIST COMPREHENSION!!!
        num_proteins = len([protein for rxn in rxns for protein in rxn["proteins"]])
        output["metadata"]["number_reviewed_enzymes"] = num_proteins



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
        filtered_proteins = len([protein for rxn in filtered_rxns for protein in rxn["proteins"]])
        output["metadata"]["number_lineage_filtered_enzymes"] = filtered_proteins


        with open("archives/"+chemical_name+".json", "w+") as f:
            f.write(json.dumps(output))
            print(str(filtered_proteins)+" enzymes for "+chemical_name+" cached in archives")
        
