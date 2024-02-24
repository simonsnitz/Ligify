import streamlit as st
import pandas as pd
import re
import base64
import streamlit.components.v1 as components
from Bio.Seq import Seq
from ligify.genbank.create_genbank import create_genbank


def download_button(object_to_download, download_filename):
    # Generates a link to download the given object_to_download.
    try:
        # some strings <-> bytes conversions necessary here
        b64 = base64.b64encode(object_to_download.encode()).decode()

    except AttributeError as e:
        b64 = base64.b64encode(object_to_download).decode()

    dl_link = f"""
    <html>
    <head>
    <title>Start Auto Download file</title>
    <script src="http://code.jquery.com/jquery-3.2.1.min.js"></script>
    <script>
    $('<a href="data:text/csv;base64,{b64}" download="{download_filename}">')[0].click()
    </script>
    </head>
    </html>
    """
    return dl_link


def download_df(ligand_name, genbank_con):

    # Figure out which direction the promoter should be facing
    enzyme_direction = st.session_state.data['protein']['context']['enzyme_direction']
    if enzyme_direction == "-":
        promoter_seq = st.session_state.data['protein']['context']["promoter"]["regulated_seq"]
    elif enzyme_direction == "+":
        promoter_seq = st.session_state.data['protein']['context']["promoter"]["regulated_seq"]
        promoter_seq = Seq(promoter_seq).reverse_complement()


    regulator_name = st.session_state.data['refseq']
    regulator_protein_seq = st.session_state.data["reg_protein_seq"]


    

    data = create_genbank(regulator_name, ligand_name, promoter_seq, regulator_protein_seq)
    '''
    Having trouble downloading the GenBank file when not locally hosted, so we'll just display it in the browser.
    '''
    # components.html(
    #     download_button(data, "pLigify_"+str(regulator_name)+".gb"),
    #     height=0,
    # )

    genbank_con.subheader("Biosensor plasmid GenBank")
    genbank_con.code(data)




def format_results(data_column, ligand_name):

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
        reg_protein_seq = st.session_state.data["reg_protein_seq"]
        reg_uniprot_id = st.session_state.data["uniprot_reg_data"]["id"]
        organism = st.session_state.data["protein"]["organism"]
        organism_name = str(organism[-2])+", "+str(organism[-1])

        reg_json = {
            "name": "Regulator attribute",
            "NCBI annotation": reg_ncbi_anno,
            "Uniprot ID": reg_uniprot_id,
            "Protein length": len(reg_protein_seq),
            "Organism": organism_name,
        }
        
            #Regulator container
        regulator_con = data_column.container()

            #GenBank container
        genbank_con = data_column.container()   

        reg_info, reg_spacer, reg_genbank = regulator_con.columns((8,1,5))
        regulator_df = pd.DataFrame(reg_json, index=["name"]).astype(str)
        regulator_df.set_index("name", inplace=True)
        regulator_df = regulator_df.T
        reg_info.subheader("Regulator information")
        reg_info.table(regulator_df)

        reg_genbank.header("")
        reg_genbank.header("")
        reg_genbank.form_submit_button(label="Show Plasmid", type="primary", on_click=download_df, args=(ligand_name,genbank_con))
        reg_genbank.markdown(f'<p style="font-size: 16px">This plasmid is designed to induce GFP expression in the presence of the target molecule via '+str(refseq)+', within E. coli</>', unsafe_allow_html=True)




        # Enzyme info
        enz_annotation = st.session_state.data['protein']['enzyme']['description']
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

        enzyme_and_lig = data_column.container()
        enz, alt_lig = enzyme_and_lig.columns([1,1])
        enz_df = pd.DataFrame(enz_json, index=[0])
        enz_df.set_index("name", inplace=True)
        enz_df = enz_df.T
        enz.subheader("Associated enzyme")
        enz.table(enz_df)


        # Alternative ligands
        alt_lig.subheader("Possible alternative ligands")
        alt_lig_list = st.session_state.data['alt_ligands'] 
        lig_html = "<ul style='height: 180px; overflow-y: scroll'>"
        for lig in alt_lig_list:
            lig_html += "<li>" + lig + "</li>"
        lig_html += "</ul>"
        alt_lig.markdown(lig_html, unsafe_allow_html=True)
            #alt_lig.markdown("- " + lig + "\n")
        #alt_lig.table(st.session_state.data['alt_ligands'])
            #TODO: Clean this table up


        # Enzyme references
        ref_and_rank = data_column.container()
        reference_col, rank_col = ref_and_rank.columns(2)
        reference_col.subheader("Enzyme references")
        for i in references:
            reference_col.markdown(f'<a target="__blank">{"https://doi.org/"+i}</a>', unsafe_allow_html=True)
        


        # Rank metrics
        rank_col.subheader("Rank description")
        rank_metrics = st.session_state.data["rank"]["metrics"]

        rank_df = pd.DataFrame(rank_metrics).T
        #orient='index')
        # rank_df = pd.DataFrame.from_dict({(i,j): rank_metrics[i][j] 
        #                    for i in rank_metrics.keys() 
        #                    for j in rank_metrics[i].keys()},
        #                orient='index')
        rank_col.dataframe(rank_df)





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
            enz_index = st.session_state.data['protein']['context']['enzyme_index']


            # Create the operon table
            genes = []
            for i in operon_json:
                if "description" in i.keys():
                    gene = {
                        "RefSeq ID": i['accession'],
                        "Description": i['description'],
                        "Direction": i['direction'],
                        "Start position": i['start'],
                        "End position": i['stop']
                    }
                else:
                    gene = {
                        "RefSeq ID": i['accession'],
                        "Description": " ",
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
