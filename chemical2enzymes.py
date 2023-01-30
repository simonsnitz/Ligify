import requests
import json
from pprint import pprint

# Input: chemical InChi-Key
# Output: A list of microbial genes known to act on that chemical

# 1. Query Rhea with InChi-Key, extract all enzymes
# 2. Filter out non-microbial enzymes
# 3. Filter by diversity, cluster the genes into groups
# 4. Fetch the operon for each representative gene
# 5. Search each operon for a regulator
# 6. Pull out any relevant literature on that regulator using PaperBLAST
# 7. Predict the operator of the regulator using GroovIO
# 8. Generate a report for each candidate regulator

#InChiKey = "KWIUHFFTVRNATP-UHFFFAOYSA-N"   # Found BetI
#InChiKey = "KBPLFHHGFOOTCA-UHFFFAOYSA-N"
InChiKey = "NIXOWILDQLNWCW-UHFFFAOYSA-M"    # Found AcuR
chemical_name = "acrylate"
domain_filter = "Bacteria"
lineage_filter_name = "Family"

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

url = "https://rest.uniprot.org/uniprotkb/search?format=json&query=reviewed:true+AND+rhea:"

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
            try:
                ncbi_id = [e["id"] for e in entry["uniProtKBCrossReferences"] if e["database"] == "RefSeq"][0]
            except:
                ncbi_id = None

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
    
    # print(proteins)
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
        family = protein["organism"][lineage_filter]
        if family not in families:
            families.append(family)
            filtered_proteins.append(protein)
    new_rxn = rxn
    new_rxn["proteins"] = filtered_proteins
    filtered_rxns.append(new_rxn)

output["rxn_data"] = filtered_rxns

filtered_proteins = len([protein for rxn in filtered_rxns for protein in rxn["proteins"]])
output["metadata"]["number_lineage_filtered_enzymes"] = filtered_proteins

with open(chemical_name+".json", "w+") as f:
    f.write(json.dumps(output))
    print("json file for "+chemical_name+" created")