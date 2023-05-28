import streamlit as st


from ligify.predict.chemical2enzymes import fetch_reactions, fetch_genes, filter_genes
from ligify.predict.enzymes2operons import pull_regulators
from ligify.predict.accID2operon import acc2operon



@st.cache_data
def fetch_data(smiles, filters):


    prog_container = st.container()
    prog_spacerL, prog, prog_spacerR = prog_container.columns((1,1,1))
    st.spinner("Processing")



    # FETCH REACTIONS

    prog_bar = prog.progress(0, text="1. Fetching reaction IDs")
    reactions = fetch_reactions(smiles = smiles)



    # FETCH GENES

    total_rxns = len(reactions["rxn_data"])
    prog_bar_increment = 20/int(total_rxns)

    counter = 0
    for i in reactions["rxn_data"]:
        prog_value = int(10+counter*prog_bar_increment)
        prog_bar.progress(prog_value, text=f"2. Fetching genes for reaction {str(counter+1)} of {str(total_rxns)} (rhea:{i['rhea_id']})")
        associated_proteins = fetch_genes(i["rhea_id"], filters["reviewed"])
        i["proteins"] = associated_proteins
        counter += 1

    # Filter homologous genes
    reactions = filter_genes(reactions, lineage_filter_name = filters["lineage"])



    # FETCH OPERONS

    if len(reactions["rxn_data"]) == 0:
        print("No enzymes found for "+str(smiles))
        pass
    else:
        counter = 0

        # Get the total number of valid protein entries for operon fetching
        # This can be shortened using list comprehension
        total_genes = 0
        for rxn in range(0,len(reactions["rxn_data"])):
            for i in range(0, len(reactions["rxn_data"][rxn]["proteins"])):
                protein = reactions["rxn_data"][rxn]["proteins"][i]
                if protein["enzyme"]["ncbi_id"] !=  None:
                    total_genes +=1

        prog_bar_increment = 50/int(total_genes)

        for rxn in range(0,len(reactions["rxn_data"])):
            for i in range(0, len(reactions["rxn_data"][rxn]["proteins"])):
                protein = reactions["rxn_data"][rxn]["proteins"][i]
                refseq_id = protein["enzyme"]["ncbi_id"]
                if refseq_id != None:

                    # Limit number of operons evaluated to avoid program taking too long to complete.
                    if counter <= filters["max_entries"]:
                        prog_value = int(30+counter*prog_bar_increment)

                        prog_bar.progress(prog_value, text=f"3. Fetching operon for gene {str(counter+1)} of {str(total_genes)} ({refseq_id})")
                        protein["context"] = acc2operon(refseq_id)
                        counter += 1



        # FETCH REGULATORS

        if reactions == None:
            prog_bar.progress(100, text="Complete.")
            return None
        
        else:

            total_regs = 0
            for rxn in reactions["rxn_data"]:
                for protein in rxn["proteins"]:
                    total_regs += 1
            prog_bar_increment = 20/int(total_regs)

            counter = 0
            regulators = []
            for rxn in reactions["rxn_data"]:
                for protein in rxn["proteins"]:
                    prog_value = int(80+counter*prog_bar_increment)
                    prog_bar.progress(prog_value, text=f"4. Fetching data for regulator {str(counter+1)} of {str(total_rxns)} ({protein['organism'][-2]}, {protein['organism'][-1]})")

                    regs = pull_regulators(protein, smiles, rxn)
                    for r in regs:
                        regulators.append(r)
                    counter += 1


            prog_bar.progress(100, text="Complete.")

            if regulators == None or len(regulators) == 0:
                return None
            else:
                return regulators

