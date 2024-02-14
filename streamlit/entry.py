
import streamlit.web.cli
import sys

from ligify.streamlit_app import run_streamlit, run_ligify


def streamlit_run():
    chem, results, prog, smiles, filters = run_streamlit()
    run_ligify(chem, results, prog, smiles, filters)

streamlit_run()