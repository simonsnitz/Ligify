import streamlit as st
import sys 
import pandas as pd


def format_display(data_column):

    if st.session_state.data:

        # Header
        refseq = st.session_state.data['refseq']
        data_column.write("")
        data_column.write("")
        data_column.markdown(f'<h1 style="text-align: center; color: black; margin-top: -50px;">{refseq}</h1>', unsafe_allow_html=True)

        data_column.divider()


        # Spacer
        data_column.text("")
        data_column.text("")


        # Regulator info
        reg_ncbi_anno = st.session_state.data["annotation"]
        reg_uniprot_anno = st.session_state.data["uniprot_reg_data"]["annotation"]
        reg_uniprot_id = st.session_state.data["uniprot_reg_data"]["id"]
        reg_length = st.session_state.data["uniprot_reg_data"]["length"]

        reg_json = {
            "name": "Regulator attribute",
            "NCBI annotation": reg_ncbi_anno,
            "Uniprot annotation": reg_uniprot_anno,
            "Uniprot ID": reg_uniprot_id,
            "Protein length": reg_length
        }
        
        regulator_con = data_column.container()
        reg_info, reg_spacer, reg_genbank, reg_spacer2 = regulator_con.columns((8,1,4,1))
        regulator_df = pd.DataFrame(reg_json, index=[0])
        regulator_df.set_index("name", inplace=True)
        regulator_df = regulator_df.T
        reg_info.subheader("Regulator information")
        reg_info.table(regulator_df)

        reg_genbank.form_submit_button(label="Download Plasmid", type="primary")
        reg_genbank.markdown(f'<p style="font-size: 16px">This plasmid is designed to induce GFP expression in the presence of the target molecule using the '+str(refseq)+' regulator</>', unsafe_allow_html=True)

        
        # download_button(
        #     label="Download Plasmid",
        #     data="genbank/pLigifyVprR.gb",
        #     file_name="pLigify"+str(st.session_state.data["refseq"])+".gb",
        #     mime="chemical/seq-na-genbank",
        # )



        # Enzyme info
        enz_annotation = st.session_state.data['protein']['enzyme']['description']
        enz_uniprot = st.session_state.data['protein']['enzyme']['uniprot_id']
        enz_refseq = st.session_state.data['protein']['enzyme']['ncbi_id']
        equation = st.session_state.data['equation'] 
        rhea_id = st.session_state.data['rhea_id'] 
        alt_ligands = st.session_state.data['alt_ligands'][0] 
        references = st.session_state.data['protein']['enzyme']['dois']

        enz_json = {"name": "Enzyme attribute",
                    "annotation": enz_annotation,
                    "equation": equation,
                    "uniprot": enz_uniprot,
                    "refseq": enz_refseq,
                    "rhea_id": rhea_id,
                    "alt_ligands": alt_ligands}

        enzyme_and_org = data_column.container()
        enz, org = enzyme_and_org.columns([1,1])
        enz_df = pd.DataFrame(enz_json, index=[0])
        enz_df.set_index("name", inplace=True)
        enz_df = enz_df.T
        enz.subheader("Associated enzyme")
        enz.table(enz_df)

        enz.text("Literature references")
        for i in references:
            enz.markdown(f'<a target="__blank">{"https://doi.org/"+i}</a>', unsafe_allow_html=True)


        # Organism info
        genome_id = st.session_state.data['protein']['context']['genome'] 
        org_json = {"name": "organism attribute",
            "genome_id": genome_id}
        phylogeny_names = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
        for i in range(0, len(st.session_state.data["protein"]["organism"])):
            org_json[phylogeny_names[i]] = st.session_state.data["protein"]["organism"][i]

        org_df = pd.DataFrame(org_json, index=[0])
        org_df.set_index("name", inplace=True)
        org_df = org_df.T
        org.subheader("Host organism")
        org.table(org_df)


        # Spacer
        data_column.text("")
        data_column.text("")



        # Operon
        operon = data_column.container()
        operon.subheader("Operon")
        genes = []

        # Get the regulator position within the operon
        operon_json = st.session_state.data['protein']['context']['operon']
        reg_index = 0
        for i in operon_json:
            if i["accession"] == refseq:
                break
            else:
                reg_index += 1


        for i in operon_json:
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
        # Highlight the regulator green
        def bg_color_col(col):
            return ['background-color: %s' % "#b3ffb0"
                        if i==reg_index
                        else ''
                    for i,x in col.items()]
        operon_df = operon_df.style.apply(bg_color_col)

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