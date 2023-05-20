import streamlit as st
import sys 
import pandas as pd
import os
import json
from pprint import pprint

from ligify import __version__ as ligify_version
from ligify.predict.pubchem import get_inchiKey
from ligify.predict.chemical2enzymes import chem2enzymes
from ligify.predict.enzymes2operons import append_operons, pull_regulators

def setup_page():
    
    st.set_page_config(page_title="Ligify", layout='wide', initial_sidebar_state='auto')
    sys.tracebacklimit = 0 #removes traceback so code is not shown during errors

    hide_streamlit_style = '''
    <style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    </style>
    '''

    st.markdown(hide_streamlit_style, unsafe_allow_html=True)

    st.markdown(f'<div style="text-align: right; font-size: 0.9em"> {ligify_version} </div>', unsafe_allow_html=True)



    #this removes the full-screen button for various elements
    style_fullscreen_button_css = """
        button[title="View fullscreen"] {
            display: none;
        }
        button[title="View fullscreen"]:hover {
            display: none;
            }
        """
    st.markdown(
        "<style>"
        + style_fullscreen_button_css
        + "</styles>",
        unsafe_allow_html=True,
    )

    # Remove border around the form
    css = r'''
        <style>
            [data-testid="stForm"] {border: 0px}
        </style>
    '''

    st.markdown(css, unsafe_allow_html=True)


    title_alignment="""
        <style>
        #the-title {
        text-align: center
        }
        </style>
        """
    st.markdown(title_alignment, unsafe_allow_html=True)


def display_data(regulators):

    for i in range(0, len(regulators)):
        with st.expander(regulators[i]["refseq"]):

            reg = regulators[i]["protein"]



                #Sensor info
            sensor = [
                        "https://www.ncbi.nlm.nih.gov/protein/"+str(regulators[i]["refseq"]), 
                        regulators[i]["annotation"],
                        ", ".join(reg["organism"])[:-2] 
                    ]
            sensor_columns = ["RefSeq link", "Annotation",  "Organism"]

            s_df = pd.DataFrame(sensor, index= sensor_columns)
            s_df.columns = ["Regulator information"]

            st.dataframe(s_df, width=700)



                #Enzyme info
            enzyme = [
                        regulators[i]["equation"],
                        regulators[i]["rhea_id"], 
                        reg["enzyme"]["description"],
                        "https://www.uniprot.org/uniprotkb/"+str(reg["enzyme"]["uniprot_id"])
                        ]
            enzyme_columns = ["Reaction", "Rhea ID", "Annotation", "Uniprot link"]

            for ref in range(0,len(reg["enzyme"]["dois"])):
                link = "https://doi.org/"+str(reg["enzyme"]["dois"][ref])
                text = "Reference "+str(ref+1)
                enzyme.append(link)
                enzyme_columns.append(text)

            e_df = pd.DataFrame(enzyme, index= enzyme_columns)
            e_df.columns = ["Enzyme information"]
            st.dataframe(e_df, width=700)


            alt_ligands = [lig for lig in regulators[i]["alt_ligands"]]
            l_df = pd.DataFrame(alt_ligands, columns=["Alternative ligand name"])
            st.dataframe(l_df, width=700)











def run_streamlit():

    setup_page()





    # #Create and name sidebar
    # st.sidebar.header('Advanced parameters')

    # def advanced_options():
    #     lineage_filter_name = st.sidebar.select_slider("Domain filter stringency", options=["Domain", "Phylum", "Class", "Order", "Family", "None"], value="Family")
    #     reviewed = st.sidebar.checkbox("Reviewed?", value=True)

    # domain_filter = "Bacteria"
    # # lineage_filter_name = st.select_slider("Domain filter stringency", options=["Domain", "Phylum", "Class", "Order", "Family", "None"], value="Family")
    # reviewed = st.checkbox("Reviewed?", value=True)


    # options = st.container()
    # options_input = options.expander("advanced options")


    # output = st.container()
    




    # col1.form("my_form")


    # HEADER
    head = st.container()
    head1, head2, head3 = head.columns(3)

    head2.markdown("<h1 style='text-align: center; color: black;'>Ligify</h1>", unsafe_allow_html=True)

    head2.subheader('Predict sensors responsive to an input ligand')
    chemical_name = head2.text_input("Chemical name", "Acrylate")


    # OPTIONS
    options = st.container()
    col1, col2, col3 = options.columns((1,3,1))

    with col2.expander("Advanced options"):

        reviewed = st.checkbox("Reviewed?", value=True)
        lineage_filter_name = st.select_slider("Domain filter stringency", options=["Domain", "Phylum", "Class", "Order", "Family", "None"], value="Family")
        max_entries = st.number_input("Max number of entries surveyed", value=25)
        filters = {"reviewed": reviewed, "lineage": lineage_filter_name, "max_entries": max_entries}


    # RESULTS
    results = st.container()
    regulator_column, data = results.columns([1,3])


    # FORM
    with st.form('ligify'):
        submit = st.container()
        submit_spacer_1, submit_button, submit_spacer_2 = submit.columns([5,1,5])

        submitted = submit_button.form_submit_button("Submit", use_container_width=True)
        if submitted:
    
            # submit_spacer_1.write("stuff")
            progress = st.container()
            prog_spacerL, prog, prog_spacerR = progress.columns([1,3,1])

            run_ligify(regulator_column, prog, chemical_name, filters)










def run_ligify(regulator_column, progress, chemical_name, filters):

    InChiKey = get_inchiKey(str(chemical_name), "name")

    st.spinner("Processing")
    my_bar = progress.progress(0, text="Fetching enzymes ...")

    # col1.write("InchiKey: "+str(InChiKey))
    
    if os.path.exists("./ligify/temp/"+str(chemical_name)+".json"):
        with open("./ligify/temp/"+str(chemical_name)+".json", "r") as f:
            regulators = json.load(f)
            print("loaded cached reg data")

            display_data(regulators)

    # else:
    # if 1 == 1:

    #     data = chem2enzymes(InChiKey = InChiKey,
    #         domain_filter = "Bacteria",
    #         lineage_filter_name = filters["lineage"], 
    #         reviewed_bool = filters["reviewed"])

    #     my_bar.progress(40, text="Fetching operons ...")

    #     data = append_operons(data, chemical_name)

    #     my_bar.progress(95, text="Fetching regulators ...")


    #     if data == None:
    #         col1.write("No regulators found")
        
    #     else:

    #         regulators = pull_regulators(data, chemical_name)
            
    #         with open("./ligify/temp/"+str(chemical_name)+".json", "w+") as f:
    #             f.write(json.dumps(regulators))
    #             print("cached regulator data")


    #         if regulators == None or len(regulators) == 0:
    #             col1.write("No regulators found")


    #         else:
    #             display_data(regulators)


            
    my_bar.progress(100, text="Complete.")


            
