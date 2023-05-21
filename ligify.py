
import click
import streamlit.web.cli
import sys

from ligify.streamlit_app import run_streamlit, run_ligify


def streamlit_run():
    chem, regulator_column, data_column, prog, chemical_name, filters = run_streamlit()
    run_ligify(chem, regulator_column, data_column, prog, chemical_name, filters)

if __name__ == '__main__':
    streamlit_run()