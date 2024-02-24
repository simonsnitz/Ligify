import streamlit as st
from streamlit_ketcher import st_ketcher
import sys 

from ligify import __version__ as ligify_version

from ligify.format_results import format_results
from ligify.fetch_data import fetch_data
from ligify.predict.database.chemical_db import blast_chemical
from ligify.predict.pubchem import get_smiles, get_inchikey, get_name



def setup_page():
    
    st.set_page_config(page_title="Ligify", layout='wide', initial_sidebar_state='auto', page_icon="images/Ligify_Favicon.png")
    #sys.tracebacklimit = 0 #removes traceback so code is not shown during errors

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
    head1, head2, head3 = head.columns((1,2,1))
    head2.image("images/Ligify_Logo.png", use_column_width=True)
    head2.subheader('Predict sensors responsive to an input ligand')

    input_mode = head2.radio(label="Select an input mode", options=["SMILES", "Name", "Draw"], horizontal=True)

    # OPTIONS
    options = st.container()
    col1, col2, col3 = options.columns((1,3,1))


    if input_mode == 'SMILES':
        smiles = head2.text_input(label="Chemical SMILES", value="C=CC(=O)[O-]", label_visibility="collapsed")
        chemical_name = get_name(smiles, "smiles")
        InChiKey = get_inchikey(smiles, "smiles")
        chemical = {"name": chemical_name, "smiles": smiles, "InChiKey": InChiKey}

    elif input_mode == 'Name':
        chemical_name = head2.text_input(label="Chemical name", value="Isoeugenol", label_visibility="collapsed")
        smiles = get_smiles(chemical_name)
        InChiKey = get_inchikey(chemical_name, "name")
        chemical = {"name": chemical_name, "smiles": smiles, "InChiKey": InChiKey}

    # elif input_mode == 'InChiKey':
    #     InChiKey = head2.text_input(label="InChiKey", value="BJIOGJUNALELMI-ONEGZZNKSA-N", label_visibility="collapsed")
    #     chemical_name = get_name(InChiKey, "inchikey")
    #     smiles = get_smiles(chemical_name)
    #     chemical = {"name": chemical_name, "smiles": smiles, "InChiKey": InChiKey}

    elif input_mode == 'Draw':
        with col2:
            smiles = st_ketcher("C=CC(=O)[O-]", height=400)
            chemical_name = get_name(smiles, "smiles")
            InChiKey = get_inchikey(chemical_name, "name")
            chemical = {"name": chemical_name, "smiles": smiles, "InChiKey": InChiKey}




    with st.sidebar:

        st.write("Bacteria respond to chemical cues using protein transcription regulators.")

        st.write("Ligify is a search tool used to mine regulators responsive to user-defined metabolites.")

        # GitHub and Email links
        st.markdown("<p style='font-size: 12px'>If you have any questions or would like to report any bugs, please contact us via <a href='mailto: simonsnitz@gmail.com'>Email</a>. \
            Our code is publically available on <a href='https://github.com/simonsnitz/Ligify'>GitHub</a>.</p>", unsafe_allow_html=True)

        #st.markdown("<div style='font-size: 12px;'>d'Oelsnitz S., Stofel S.K., and Ellington A.D. (2023) Snowprint: a predictive tool for genetic biosensor discovery. \
        #            <i>bioRxiv</i> <b>DOI:</b><a href='https://www.biorxiv.org/content/10.1101/2023.04.29.538814v1'>10.1101/2023.04.29.538814v1</a></div> <br>", unsafe_allow_html=True)

        st.markdown("<p style='font-size: 12px'>Ligify development was supported by the National Institute of Standards and Technology (70NANB21H100)", unsafe_allow_html=True)

        st.divider()

        adv_options = st.container()
        adv_options.header("Advanced options")

        adv_options.write("Editing these changes processing time and data returned")
        adv_options.divider()

        adv_options.write("Fetch reactions")
        max_reactions = adv_options.number_input("Max number of reactions evaluated", value=20)
        adv_options.divider()

        adv_options.write("Fetch genes")
        proteins_per_reaction = adv_options.number_input("Max number of proteins fetched per reaction", value=20)
        reviewed = adv_options.checkbox("Reviewed only", value=True)
        adv_options.divider()

        adv_options.write("Fetch operons")
        max_operons = adv_options.number_input("Max number of operons evaluated", value=20)
        lineage_filter_name = adv_options.selectbox("Domain filter stringency", options=["Domain", "Phylum", "Class", "Order", "Family", "Genus", "None"], index=4)
        adv_options.divider()

        adv_options.write("Fetch regulators")
        max_alt_chems = adv_options.number_input("Max suggested chemicals if failed", value=10)

        filters = { 
                    "max_reactions": max_reactions,
                    "proteins_per_reaction": proteins_per_reaction, 
                    "reviewed": reviewed, 
                    "lineage": lineage_filter_name, 
                    "max_operons": max_operons,
                    "max_alt_chems": max_alt_chems,
                    }





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


        return chem, results, prog, chemical, filters








# This is essentially the frontend component of the Ligify Web applications.

def run_ligify(chem, results, progress, chemical, filters):

    m_spacer1, metrics_col, m_spacer2 = results.container().columns((1,3,1))
    regulator_column, data_column = results.columns([1,3])

    if st.session_state.SUBMITTED:

        if chemical["smiles"] == None:
            data_column.subheader("Chemical input was not recognized. Please try a different input method.")

        else:

            # SMILES = str(chemical_smiles)
            # chem.image(f'http://hulab.rxnfinder.org/smi2img/{SMILES}/', width=200)

            regulators, metrics = fetch_data(chemical["InChiKey"], filters)

            select_spacerL, please_select, select_spacerR  = data_column.container().columns([1,2,1])     

            format_results(data_column, chemical["name"])

            regulator_column.header('')
            regulator_column.subheader('Sensor candidates')
            regulator_column.divider()

            # If no regulators are returned, suggest alternative queries
            if regulators == None:
                similar_chemicals = blast_chemical(chemical["smiles"], filters["max_alt_chems"])
                regulator_column.write("No regulators found")
                please_select.subheader("No associated reactions   :pensive:") 
                please_select.write("Consider an alternative query")   
                data_column.dataframe(similar_chemicals)
                
            # If regulators are returned, format display
            else:

                # Metrics data
                metrics_col.subheader("Search metrics")
                m_rhea, m_genes, m_filtered, m_operons, m_regs = metrics_col.columns(5)
                m_rhea.metric("Rhea reactions", metrics["RHEA Reactions"])
                m_genes.metric("Bacterial genes", metrics["Total genes"])
                m_filtered.metric("Filtered genes", metrics["Filtered genes"])
                m_operons.metric("Operons", metrics["Total operons"])
                m_regs.metric("Regulators", metrics["Total regulators"])
                metrics_col.divider()

                if not st.session_state.data:
                    please_select.subheader("Please select a regulator") 

                reg_acc_col, rank_col = regulator_column.columns((2,1))
                reg_acc_col.markdown("<h5>Regulator</h5>", unsafe_allow_html=True)
                rank_col.markdown("<h5>Rank</h5>", unsafe_allow_html=True)

                for i in range(0, len(regulators)):
                    name = "var"+str(i)
                    rank = regulators[i]["rank"]["rank"]
                    color = regulators[i]["rank"]["color"]
                    rank_col.markdown(f"<p style='font-size:20px; font-weight: 600; color: {color};'>{rank}</p>", unsafe_allow_html=True)
                    name = reg_acc_col.form_submit_button(regulators[i]['refseq'])
                    if name:
                        st.session_state.data = regulators[i]
                        st.experimental_rerun()

