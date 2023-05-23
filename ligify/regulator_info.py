import streamlit as st
import sys 
import pandas as pd


def format_display(data_column):

    if st.session_state.SUBMITTED and not st.session_state.data:

        select_container = data_column.container()
        select_spacerL, please_select, select_spacerR  = select_container.columns([1,2,1])
        please_select.subheader("Please select a regulator")


    elif st.session_state.data:

        # Regulator info
        refseq = st.session_state.data['refseq']
        annotation = st.session_state.data['annotation'] 

        # Enzyme info
        enz_annotation = st.session_state.data['protein']['enzyme']['description']
        enz_uniprot = st.session_state.data['protein']['enzyme']['uniprot_id']
        enz_refseq = st.session_state.data['protein']['enzyme']['ncbi_id']
        equation = st.session_state.data['equation'] 
        rhea_id = st.session_state.data['rhea_id'] 
        references = st.session_state.data['protein']['enzyme']['dois']

        enz_json = {"name": "enzyme attribute",
                    "annotation": enz_annotation,
                    "equation": equation,
                    "uniprot": enz_uniprot,
                    "refseq": enz_refseq,
                    "rhea_id": rhea_id}

        # Organism info
        genome_id = st.session_state.data['protein']['context']['genome'] 
        kingdom = st.session_state.data['protein']['organism'][0]
        phylum = st.session_state.data['protein']['organism'][1]
        org_class = st.session_state.data['protein']['organism'][2]
        order = st.session_state.data['protein']['organism'][3]
        family = st.session_state.data['protein']['organism'][4]
        genus = st.session_state.data['protein']['organism'][5]

        org_json = {"name": "organism attribute",
                    "genome_id": genome_id,
                    "kingdom": kingdom,
                    "phylum": phylum,
                    "class": org_class,
                    "order": order,
                    "family": family,
                    "genus": genus}



        data_column.markdown(f'<h1 style="text-align: center; color: black; margin-top: -50px;">{refseq}</h1>', unsafe_allow_html=True)
        data_column.markdown(f'<h3 style="text-align: center; color: black;">{annotation}</h3>', unsafe_allow_html=True)


        data_column.divider()


        # Spacer
        data_column.text("")
        data_column.text("")


        # Enzyme info
        enzyme_and_org = data_column.container()
        enz, org = enzyme_and_org.columns([1,1])
        enz_df = pd.DataFrame(enz_json, index=[0])
        enz_df.set_index("name", inplace=True)
        enz_df = enz_df.T
        enz.subheader("Associated enzyme")
        enz.table(enz_df)

        enz.text("Literature")
        for i in references:
            enz.markdown(f'<a target="__blank">{"https://doi.org/"+i}</a>', unsafe_allow_html=True)


        # Organism info
        org_df = pd.DataFrame(org_json, index=[0])
        org_df.set_index("name", inplace=True)
        org_df = org_df.T
        org.subheader("Host organsim")
        org.table(org_df)


        # Spacer
        data_column.text("")
        data_column.text("")



        # Operon
        operon = data_column.container()
        operon.subheader("Operon")
        genes = []
        
        for i in st.session_state.data['protein']['context']['operon']:
            try:
                gene = {
                    "alias": i['alias'],
                    "description": i['description'],
                    "refseq": i['accession'],
                    "direction": i['direction'],
                    "start_position": i['start'],
                    "end_position": i['stop']
                }
                genes.append(gene)
            except:
                # Sometimes the alias isn't returned
                gene = {
                    "description": i['description'],
                    "refseq": i['accession'],
                    "direction": i['direction'],
                    "start_position": i['start'],
                    "end_position": i['stop']
                }
                genes.append(gene)

        operon_df = pd.DataFrame(genes)
        operon.table(operon_df)




        # Spacer
        # data_column.text("")
        # data_column.text("")


        # This takes too long to load. It's a 'nice to have' that's not worth the lag it creates.

        # Alternative ligands
            # This takes a while to load if there are a lot of ligands
        # alt_ligands = data_column.container()
        # alt_ligands.subheader("Possible alternative ligands")
        # a_ligands_smiles = []
        # a_ligands_names = []
        
        # for i in st.session_state.data['alt_ligands']:
        #     try:
        #         SMILES = get_smiles(str(i))
        #         url = f'http://hulab.rxnfinder.org/smi2img/{SMILES}/'
        #         if check_url(url):
        #             a_ligands_smiles.append(url)
        #             a_ligands_names.append(str(i))
        #     except:
        #         pass


        # alt_ligands.image(a_ligands_smiles, width=200, caption= a_ligands_names)