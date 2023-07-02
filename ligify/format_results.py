import streamlit as st
import sys 
import pandas as pd
import re

def format_results(data_column):

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
        reg_uniprot_id = st.session_state.data["uniprot_reg_data"]["id"]
        reg_length = st.session_state.data["uniprot_reg_data"]["length"]
        organism = st.session_state.data["protein"]["organism"]
        organism_name = str(organism[-2])+", "+str(organism[-1])

        reg_json = {
            "name": "Regulator attribute",
            "NCBI annotation": reg_ncbi_anno,
            "Uniprot ID": reg_uniprot_id,
            "Protein length": reg_length,
            "Organism": organism_name,
        }
        
        regulator_con = data_column.container()
        reg_info, reg_spacer, reg_genbank, reg_spacer2 = regulator_con.columns((8,1,4,1))
        regulator_df = pd.DataFrame(reg_json, index=[0])
        regulator_df.set_index("name", inplace=True)
        regulator_df = regulator_df.T
        reg_info.subheader("Regulator information")
        reg_info.table(regulator_df)

        reg_genbank.header("")
        reg_genbank.header("")
        reg_genbank.form_submit_button(label="Download Plasmid", type="primary")
        reg_genbank.markdown(f'<p style="font-size: 16px">This plasmid is designed to induce GFP expression in the presence of the target molecule via '+str(refseq)+', within E. coli</>', unsafe_allow_html=True)

        # TODO:
            # MAKE A FUNCTIONAL DOWNLOAD BUTTON

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
        references = st.session_state.data['protein']['enzyme']['dois']

        enz_json = {"name": "Enzyme attribute",
                    "Annotation": enz_annotation,
                    "Reaction": equation,
                    "RHEA ID": rhea_id,
                    "RefSeq ID": enz_refseq,
                    }

        enzyme_and_org = data_column.container()
        enz, alt_lig = enzyme_and_org.columns([1,1])
        enz_df = pd.DataFrame(enz_json, index=[0])
        enz_df.set_index("name", inplace=True)
        enz_df = enz_df.T
        enz.subheader("Associated enzyme")
        enz.table(enz_df)

        enzyme_and_org.subheader("Enzyme references")
        for i in references:
            enzyme_and_org.markdown(f'<a target="__blank">{"https://doi.org/"+i}</a>', unsafe_allow_html=True)




        # Alternative ligands
        alt_lig.subheader("Possible alternative ligands")
        alt_lig.table(st.session_state.data['alt_ligands'])
            #TODO: Clean this table up

        # Spacer
        data_column.text("")
        data_column.text("")




        # TODO:
        # Add a description for the rank given.




        # Operon
        operon_seq = data_column.container()
        operon_seq.subheader("Operon")
        operon_data = st.session_state.data['protein']['context']
        operon_json = operon_data['operon']
        
        with operon_seq.expander(label="Click here to see the genetic context for "+str(refseq)):


            # Get the regulator position within the operon
            reg_index = 0
            for i in operon_json:
                if i["accession"] == refseq:
                    break
                else:
                    reg_index += 1
            # Get the enzyme position within the operon
            enz_index = st.session_state.data['protein']['context']['protein_index']

            # Create the operon table
            genes = []
            for i in operon_json:
                gene = {
                    "RefSeq ID": i['accession'],
                    "Description": i['description'],
                    "Direction": i['direction'],
                    "Start position": i['start'],
                    "End position": i['stop']
                }
                genes.append(gene)
            operon_df = pd.DataFrame(genes)

            # Color the operon table
            colors = ["#fcb1b1", "#e6cffc", "#9afcac", "#fcc99d", "#a3fff6", "#fdff9c", "#ccd4fc", "#fcbdf6"]
            def bg_color_col(col):
                color = [colors[i % len(colors)] for i,x in col.items()]
                return ['background-color: %s' % color[i]
                            if col.name=='RefSeq ID' or (col.name=="Description" and i==reg_index) or (col.name=="Description" and i==enz_index)
                            else ''
                        for i,x in col.items()]

            operon_df = operon_df.style.apply(bg_color_col)
            st.table(operon_df)


            # Display the predicted promoter
            st.markdown("<h5>Predicted promoter</h5>", unsafe_allow_html=True)
            st.write(operon_data["promoter"]['regulated_seq'])

            # Create and display the color-annotated genome fragment
            operon_seq = ""

            c = 0
            for seq in operon_data["operon_seq"]:
                sequence = operon_data["operon_seq"][seq]
                if re.compile(r"spacer").search(seq):
                    html = "<span style='color: grey;'>"+str(sequence)+"</span>"
                elif re.compile(r"overlap").search(seq):
                    html = "<span style='color: red;'>"+str(sequence)+"</span>"
                elif re.compile(r"fwd").search(seq):
                    html = f"<b><span style='background: {colors[c % len(colors)]};'>"+str(sequence)+"</span></b>"
                    c += 1
                else:
                    html = f"<span style='background: {colors[c % len(colors)]};'>"+str(sequence)+"</span>"
                    c += 1
                operon_seq += html
            st.markdown("<h5>Full operon sequence</h5>", unsafe_allow_html=True)
            st.markdown(operon_seq, unsafe_allow_html=True)
