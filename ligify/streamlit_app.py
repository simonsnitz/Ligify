import streamlit as st
from streamlit_ketcher import st_ketcher
import sys 

from ligify import __version__ as ligify_version

from ligify.regulator_info import format_display
from ligify.fetch_data import fetch_data



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


# Add a sidebar for documentation?
    # with st.sidebar:
    #     st.header("How to use")
    #     st.divider()
    #     st.header("About")
    #     st.divider()
    #     st.header("FAQ")


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

    input_mode = head2.radio(label="Select an input mode", options=["SMILES", "Name", "Draw"], horizontal=True)

    # OPTIONS
    options = st.container()
    col1, col2, col3 = options.columns((1,3,1))


    if input_mode == 'SMILES':
        smiles = head2.text_input("", "C=CC(=O)[O-]", label_visibility="collapsed")
    elif input_mode == 'Name':
        smiles = head2.text_input("", "Acrylate", label_visibility="collapsed")
    elif input_mode == 'Draw':
        with col2:
            smiles = st_ketcher("C=CC(=O)[O-]", height=400)




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





        return chem, regulator_column, data_column, prog, smiles, filters










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

