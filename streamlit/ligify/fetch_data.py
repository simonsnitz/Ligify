import streamlit as st
import json

from ligify.predict.chemical2enzymes import fetch_reactions, fetch_genes, filter_genes
from ligify.predict.enzymes2operons import pull_regulators
from ligify.predict.accID2operon import acc2operon
from ligify.predict.rank import calculate_rank


@st.cache_data
def fetch_data(InChiKey, filters):

    prog_container = st.container()
    prog_spacerL, prog, prog_spacerR = prog_container.columns((1,1,1))
    st.spinner("Processing")

    metrics = {}

    # FETCH REACTIONS

    prog_bar = prog.progress(0, text="1. Fetching reaction IDs")
    reactions = fetch_reactions(InChiKey = InChiKey, max_reactions = filters["max_reactions"])
    total_rxns = len(reactions["rxn_data"])
    metrics["RHEA Reactions"] = total_rxns

    if total_rxns > 0:

        # FETCH GENES

        prog_bar_increment = 20/int(total_rxns)

        counter = 0
        for i in reactions["rxn_data"]:
            prog_value = int(10+counter*prog_bar_increment)
            prog_bar.progress(prog_value, text=f"2. Fetching genes for reaction {str(counter+1)} of {str(total_rxns)} (rhea:{i['rhea_id']})")
            associated_proteins = fetch_genes(i["rhea_id"], filters["reviewed"], filters["proteins_per_reaction"])
            i["proteins"] = associated_proteins
            counter += 1

        metrics["Total genes"] = sum([len(i["proteins"]) for i in reactions["rxn_data"]])

        # Filter homologous genes
        reactions = filter_genes(reactions, lineage_filter_name = filters["lineage"])
        metrics["Filtered genes"] = sum([len(i["proteins"]) for i in reactions["rxn_data"]])


        # FETCH OPERONS

        if len(reactions["rxn_data"]) == 0:
            print("No enzymes found")
            prog_bar.empty()
            return None, None
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

            prog_bar_increment = 50/int(total_genes)

            for rxn in range(0,len(reactions["rxn_data"])):
                for i in range(0, len(reactions["rxn_data"][rxn]["proteins"])):
                    protein = reactions["rxn_data"][rxn]["proteins"][i]
                    refseq_id = protein["enzyme"]["ncbi_id"]
                    if refseq_id != None:

                        # Limit number of operons evaluated to avoid program taking too long to complete.
                        if operon_counter <= filters["max_operons"]:
                            prog_value = int(30+operon_counter*prog_bar_increment)

                            prog_bar.progress(prog_value, text=f"3. Fetching operon for gene {str(operon_counter+1)} of {str(total_genes)} ({refseq_id})")
                            protein["context"] = acc2operon(refseq_id)
                            operon_counter += 1

            metrics["Total operons"] = operon_counter


            # FETCH REGULATORS

            if reactions == None:
                prog_bar.empty()
                return None, None
            
            else:
                total_regs = 0
                for rxn in reactions["rxn_data"]:
                    for protein in rxn["proteins"]:
                        total_regs += 1
                prog_bar_increment = 20/int(total_regs)


                # This is where all of the display data is created
                counter = 0
                regulators = []
                

                for rxn in reactions["rxn_data"]:
                    for protein in rxn["proteins"]:
                        prog_value = int(80+counter*prog_bar_increment)
                        prog_bar.progress(prog_value, text=f"4. Fetching data for regulator {str(counter+1)} of {str(total_regs)} ({protein['organism'][-2]}, {protein['organism'][-1]})")

                        regs = pull_regulators(protein, rxn)
                        for r in regs:
                            regulators.append(r)
                        counter += 1

                metrics["Total regulators"] = len(regulators)
                prog_bar.empty()



                # Filter out duplicate regulators
                refseq_ids = []
                filtered_regulators = []
                for i in regulators:
                    if i["refseq"] not in refseq_ids:
                        filtered_regulators.append(i)
                        refseq_ids.append(i["refseq"])

                # Filter out regulators without a predicted promoter
                filtered_regulators = [i for i in filtered_regulators if i["protein"]["context"]["promoter"] != None]

                # Create a rank for each regulator
                for r in filtered_regulators:
                    rank = calculate_rank(r)
                    r["rank"] = rank


                if filtered_regulators == None or len(filtered_regulators) == 0:
                    return None, None
                else:
                    return filtered_regulators, metrics


    else:
        prog_bar.empty()
        return None, None

