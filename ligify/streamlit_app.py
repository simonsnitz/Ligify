import streamlit as st
import sys 
import pandas as pd


from ligify import __version__ as ligify_version
from ligify.predict.pubchem import get_inchiKey, get_smiles, check_url
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




        # data_column.write(st.session_state.data)










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
    chemical_name = head2.text_input("Chemical name", "Acrylate")




    # OPTIONS
    options = st.container()
    col1, col2, col3 = options.columns((1,3,1))

    with col2.expander("Advanced options"):

        adv_options = st.container()
        option1, option2, option3 = adv_options.columns((1,3,2))
        reviewed = option1.checkbox("Reviewed only", value=True)
        lineage_filter_name = option2.selectbox("Domain filter stringency", options=["Domain", "Phylum", "Class", "Order", "Family", "None"], index=4)
        max_entries = option3.number_input("Max number of entries surveyed", value=25)
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



        return chem, regulator_column, data_column, prog, chemical_name, filters






@st.cache_data
def fetch_data(chemical_name, InChiKey, filters):
    # if os.path.exists("./ligify/temp/"+str(chemical_name)+".json"):
    #     with open("./ligify/temp/"+str(chemical_name)+".json", "r") as f:
    #         regulators = json.load(f)
    #         print("loaded cached reg data")
            
    #     return regulators
    
    if 1 == 1:

        prog_container = st.container()
        prog_spacerL, prog, prog_spacerR = prog_container.columns((1,1,1))
        st.spinner("Processing")
        prog_bar = prog.progress(0, text="Fetching enzymes ...")

        data = chem2enzymes(InChiKey = InChiKey,
            domain_filter = "Bacteria",
            lineage_filter_name = filters["lineage"], 
            reviewed_bool = filters["reviewed"])


        prog_bar.progress(40, text="Fetching operons ...")

        data = append_operons(data, chemical_name)

        prog_bar.progress(95, text="Fetching regulators ...")


        if data == None:
            prog_bar.progress(100, text="Complete.")
            return None
            # _regulator_column.write("No regulators found")
        
        else:

            regulators = pull_regulators(data, chemical_name)
            prog_bar.progress(100, text="Complete.")

            
            # with open("./ligify/temp/"+str(chemical_name)+".json", "w+") as f:
            #     f.write(json.dumps(regulators))
            #     print("cached regulator data")

            if regulators == None or len(regulators) == 0:
                return None
                # _regulator_column.write("No regulators found")


            else:
                return regulators








def run_ligify(chem, regulator_column, data_column, progress, chemical_name, filters):

    if st.session_state.SUBMITTED:


        InChiKey = get_inchiKey(str(chemical_name), "name")

        SMILES = get_smiles(str(chemical_name))
        chem.image(f'http://hulab.rxnfinder.org/smi2img/{SMILES}/', width=300)



        regulators = fetch_data(chemical_name, InChiKey, filters)
        # if os.path.exists("./ligify/temp/"+str(chemical_name)+".json"):
        #     with open("./ligify/temp/"+str(chemical_name)+".json", "r") as f:
        #         regulators = json.load(f)
        #         print("loaded cached reg data")

        format_display(data_column)

        regulator_column.subheader(f'{chemical_name} sensor candidates')
        regulator_column.divider()

        for i in range(0, len(regulators)):
            name = "var"+str(i)
            name = regulator_column.form_submit_button(regulators[i]['refseq'])
            if name:
                st.session_state.data = regulators[i]
                st.experimental_rerun()

