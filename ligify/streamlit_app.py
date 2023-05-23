import streamlit as st
from streamlit_ketcher import st_ketcher
import sys 
import pandas as pd

from ligify import __version__ as ligify_version

from ligify.regulator_info import format_display
from ligify.predict.chemical2enzymes import fetch_reactions, fetch_genes
from ligify.predict.enzymes2operons import pull_regulators
from ligify.predict.accID2operon import acc2operon

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

    st.markdown(f'<div style="text-align: right; font-size: 0.9em"> Version: {ligify_version} </div>', unsafe_allow_html=True)



    # Removes the full-screen button for various elements
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

    # Removes border around forms
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








def run_streamlit():

    setup_page()


    # The below code really just sets up the input section at the top


    # Initialize state variables
    if "data" not in st.session_state:
        st.session_state.data = False

    if 'SUBMITTED' not in st.session_state:
        st.session_state.SUBMITTED =  False


    def _connect_form_cb(connect_status):
        st.session_state.SUBMITTED = connect_status
        st.session_state.data = False


    # HEADER
    head = st.container()
    head1, head2, head3 = head.columns(3)

    head2.markdown("<h1 style='text-align: center; color: black;'>Ligify</h1>", unsafe_allow_html=True)
    head2.subheader('Predict sensors responsive to an input ligand')

    chemical_smiles = head2.text_input("Molecule", "C=CC(=O)[O-]")


    # OPTIONS
    options = st.container()
    col1, col2, col3 = options.columns((1,3,1))

    with col2:
        smiles_code = st_ketcher(chemical_smiles, height=400)

    with col2.expander("Advanced options"):

        adv_options = st.container()
        option1, option2, option3 = adv_options.columns((1,3,2))
        reviewed = option1.checkbox("Reviewed only", value=True)
        lineage_filter_name = option2.selectbox("Domain filter stringency", options=["Domain", "Phylum", "Class", "Order", "Family", "None"], index=4)
        max_entries = option3.number_input("Max number of entries surveyed", value=20)
        filters = {"reviewed": reviewed, "lineage": lineage_filter_name, "max_entries": max_entries}





    # FORM
    with st.form(key='ligify'):

        # SUBMIT BUTTON
        submit = st.container()
        submit_spacer_1, submit_button, submit_spacer_2 = submit.columns([5,1,5])
        submitted = submit_button.form_submit_button("Submit", use_container_width=True, on_click=_connect_form_cb, args=(True,))

        # PROGRESS BAR
        progress = st.container()
        prog_spacerL, prog, prog_spacerR = progress.columns([1,3,1])

        # CHEMICAL_STRUCTURE
        structure = st.container()
        chem_spacerL, chem, chem_spacerR = structure.columns([2,1,2])

        # RESULTS
        results = st.container()
        regulator_column, data_column = results.columns([1,3])



        return chem, regulator_column, data_column, prog, chemical_smiles, filters






@st.cache_data
def fetch_data(chemical_smiles, filters):


    prog_container = st.container()
    prog_spacerL, prog, prog_spacerR = prog_container.columns((1,1,1))
    st.spinner("Processing")


    prog_bar = prog.progress(0, text="1. Fetching reaction IDs ...")
    reactions = fetch_reactions(smiles = chemical_smiles)


    prog_bar.progress(10, text="2. Fetching associated genes ...")
    data = fetch_genes(reactions, lineage_filter_name = filters["lineage"], reviewed_bool = filters["reviewed"])


    if len(data["rxn_data"]) == 0:
        print("No enzymes found for "+str(chemical_smiles))
        pass
    else:
        counter = 0

        # Get the total number of valid protein entries for operon fetching
        # This can be shortened using list comprehension
        total_entries = 0
        for rxn in range(0,len(data["rxn_data"])):
            for i in range(0, len(data["rxn_data"][rxn]["proteins"])):
                protein = data["rxn_data"][rxn]["proteins"][i]
                if protein["enzyme"]["ncbi_id"] !=  None:
                    total_entries +=1
        prog_bar_increment = 60/int(total_entries)

        for rxn in range(0,len(data["rxn_data"])):
            for i in range(0, len(data["rxn_data"][rxn]["proteins"])):
                protein = data["rxn_data"][rxn]["proteins"][i]
                refseq_id = protein["enzyme"]["ncbi_id"]
                if refseq_id != None:

                    # Limit number of operons evaluated to avoid program taking too long to complete.
                    if counter <= filters["max_entries"]:
                        prog_value = int(30+counter*prog_bar_increment)
                        prog_bar.progress(prog_value, text=f"3. Fetching data for operon {str(counter+1)} of {str(total_entries)}")
                        protein["context"] = acc2operon(refseq_id)
                        counter += 1


        prog_bar.progress(90, text="Fetching regulators ...")

        if data == None:
            prog_bar.progress(100, text="Complete.")
            return None
        
        else:
            regulators = pull_regulators(data, chemical_smiles)
            prog_bar.progress(100, text="Complete.")

            if regulators == None or len(regulators) == 0:
                return None
            else:
                return regulators








def run_ligify(chem, regulator_column, data_column, progress, chemical_smiles, filters):

    if st.session_state.SUBMITTED:


        # SMILES = get_smiles(str(chemical_smiles))
        # chem.image(f'http://hulab.rxnfinder.org/smi2img/{SMILES}/', width=300)



        regulators = fetch_data(chemical_smiles, filters)
        # if os.path.exists("./ligify/temp/"+str(chemical_name)+".json"):
        #     with open("./ligify/temp/"+str(chemical_name)+".json", "r") as f:
        #         regulators = json.load(f)
        #         print("loaded cached reg data")

        format_display(data_column)

        regulator_column.subheader(f'{chemical_smiles} sensor candidates')
        regulator_column.divider()

        if regulators == None:
            regulator_column.write("No regulators found")
            
        else:
            for i in range(0, len(regulators)):
                name = "var"+str(i)
                name = regulator_column.form_submit_button(regulators[i]['refseq'])
                if name:
                    st.session_state.data = regulators[i]
                    st.experimental_rerun()

