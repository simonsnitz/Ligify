import requests
import json
import os
from pprint import pprint
import time
from libchebipy._chebi_entity import ChebiEntity
from math import ceil
import sys
sys.path.append("..") # Adds higher directory to python modules path.
from accID2operon import acc2operon
# from pubchem import get_inchiKey


# TODO:
# Get Rhea IDs and equations from PubChem. 


# Create a json file with all rhea reactions and associated chemical IDs for each
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




# Format the chemical database to link from rhea -> chem ID, to chem ID -> rhea
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




# filter chemical database to drug-like molecules that likely can cross the cell membrane
def lipinski_filter():
    with open("formatted.json", "r") as f:

        data = json.load(f)

        num_batches = ceil(len(data)/1000)

        chem_json = {}
        
        for i in range(0, num_batches):
            
            # extract the group of 1000 relevant chemical IDs
            counter = 0
            start = i*1000
            end = (i+1)*1000
            chebis = ""
            for chebi in data:
                if counter >= start and counter < end:
                    chebis += "CHEBI:"+str(chebi)+","
                counter += 1

            chebis = chebis[0:-1]

            # get chemical data from an API request to mychem
            params = {'ids':chebis, 'fields':\
                'pubchem.xlogp,\
                pubchem.formal_charge,\
                pubchem.molecular_weight,\
                pubchem.hydrogen_bond_acceptor_count,\
                pubchem.hydrogen_bond_donor_count,\
                pubchem.topological_polar_surface_area,\
                pubchem.rotatable_bond_count,\
                pubchem.smiles.isomeric,\
                pubchem.smiles.canonical,\
                chebi.name'}

            res = requests.post('http://mychem.info/v1/chem', params)
            con = res.json()
            
            t = 0
            a = 0

            for chem in con:

                chebid = chem['query'].strip("CHEBI:")
                rheas = data[chebid]
                try:
                    name = chem["chebi"]["name"]
                except:
                    try:
                        name = chem["chebi"][0]["name"]
                    except:
                        name = "unknown"

                try:
                    charge = chem['pubchem']['formal_charge']
                    mass = chem['pubchem']['molecular_weight']
                    hbond_acceptor = chem["pubchem"]['hydrogen_bond_acceptor_count']
                    hbond_donor = chem["pubchem"]["hydrogen_bond_donor_count"]
                    rbonds = chem['pubchem']['rotatable_bond_count']
                    polar_sa = chem['pubchem']['topological_polar_surface_area']
                    logp = chem['pubchem']['xlogp']
                    try:
                        smiles = chem["pubchem"]["smiles"]["isomeric"]
                    except:
                        smiles = chem["pubchem"]["smiles"]["canonical"]
                    t += 1

                    if logp < 5 and mass < 500 and (hbond_acceptor+hbond_donor) < 12 and rbonds < 10 and polar_sa < 140 and abs(charge) < 2:
                        # add to new dictionary
                        chem_json[chebid] = { 
                            "rheas": rheas,
                            "name": name,
                            "smiles": smiles }
                        a += 1


                except:
                    try:
                        charge = chem['pubchem']['formal_charge']
                        mass = chem['pubchem']['molecular_weight']
                        hbond_acceptor = chem["pubchem"]['hydrogen_bond_acceptor_count']
                        hbond_donor = chem["pubchem"]["hydrogen_bond_donor_count"]
                        rbonds = chem['pubchem']['rotatable_bond_count']
                        polar_sa = chem['pubchem']['topological_polar_surface_area']
                        try:
                            smiles = chem["pubchem"]["smiles"]["isomeric"]
                        except:
                            smiles = chem["pubchem"]["smiles"]["canonical"]
                        t += 1

                        if mass < 500 and (hbond_acceptor+hbond_donor) < 12 and rbonds < 10 and polar_sa < 140 and abs(charge) < 2:
                            # add to new dictionary
                            chem_json[chebid] = { 
                                "rheas": rheas,
                                "name": name,
                                "smiles": smiles }
                            a += 1                   
                    except:
                        pass                

            print('tested: '+str(t))
            print("added: "+ str(a))

        with open("drug_like.json", "w+") as o:

            out = json.dumps(chem_json)
            o.write(out)
            print('saved drug-like molecules')






def append_genes():

    url = "https://rest.uniprot.org/uniprotkb/search?format=json&size=25&query=reviewed:true+AND+"


    with open("drug_like.json", "r") as f:
            
        db = json.load(f)

        with open("all_proteins.json", "w") as out:

            db_with_proteins = {}

            counter = 0
            for chemical in db:
                if counter < 3000:
                    rhea_ids = [i.strip("RHEA:") for i in db[chemical]["rheas"]]
                    

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
                



                        if len(proteins) != 0:
                            db_with_proteins[chemical] = {
                                "name": db[chemical]["name"],
                                "smiles": db[chemical]["smiles"],
                                "rheas": db[chemical]["rheas"],
                                "proteins":proteins
                            }
                            print("finished "+str(counter)+" of "+(str(len(db))))

                counter += 1
            out.write(json.dumps(db_with_proteins))
            print("saved data to all_proteins.json")





def filter_genes():

    with open("all_proteins.json", "r") as f:

        data = json.load(f)

        print("Total number of chemicals: "+str(len(data)))


        filtered = {}
        for chem in data:
            diverse_proteins = []
            families = []
            for protein in data[chem]["proteins"]:
                    # Sometimes the family name isn't provided in the lineage returned.
                try:
                    family = protein["organism"][3]
                except:
                    family = ""
                if family not in families:
                    families.append(family)
                    diverse_proteins.append(protein)
            new_rxn = data[chem]
            new_rxn["proteins"] = diverse_proteins
            filtered[chem] = new_rxn

        with open("filtered_proteins.json", "w+") as out_file:
            out_file.write(json.dumps(filtered))
            print("saved filtered data in filtered_proteins.json")




def fetch_operons():

    with open("filtered_proteins.json", "r") as f:

        og_data = json.load(f)
        # if os.path.exists("with_operons.json") == False:
        #     print("creating a new with_operons output file")
        #     with open("with_operons.json", "w") as o:
        #         pass


        with open("with_operons.json", "r+") as out:
            data_with_operons = json.load(out)
            
            num_proteins = 0
            for i in og_data:
                for p in og_data[i]["proteins"]:
                    num_proteins += 1


            c = 0
            for chem in og_data:
                protein = og_data[chem]["proteins"]
                for protein_index in range(0, len(og_data[chem]["proteins"])):
                    if c < 3000:
                        if chem in data_with_operons:
                            if "context" in data_with_operons[chem]["proteins"][protein_index]:
                                pass
                            else:
                                refseq = protein[protein_index]["enzyme"]["ncbi_id"]
                                if refseq != None:
                                    operon = acc2operon(refseq)
                                    protein[protein_index]["context"] = operon
                                    data_with_operons[chem] = og_data[chem]
                        else:
                            refseq = protein[protein_index]["enzyme"]["ncbi_id"]
                            if refseq != None:
                                operon = acc2operon(refseq)
                                protein[protein_index]["context"] = operon
                                data_with_operons[chem] = og_data[chem]

                        print("fetched protein "+str(c)+" of "+str(num_proteins))
                        c += 1
                

                out.seek(0)
                out.write(json.dumps(data_with_operons))
                out.truncate()            

                




def pull_regulators(protein, chemical_name, rxn):


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

    # pprint(fetch_reactions("C=CC(=O)[O-]"))

    # initialize_database()
    
    # format_data()

    #lipinski_filter()

    # append_genes()

    # filter_genes()

    # fetch_operons()

    pull_regulators()



    # with open("all_proteins.json", "r") as f:
    #     data = json.load(f)

    #     c = 0
    #     for chem in data:
    #         for protein in data[chem]["proteins"]:
    #             c += 1
    #     print(c)





    #     with open("names.txt", "w+") as out_file:
    #         out = ""
    #         for i in data:
    #             out += str(data[i]["name"])+"\n"

    #         out_file.write(out)


    

