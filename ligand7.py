
import click
import streamlit.web.cli
import sys

from ligand7.streamlit_app import run_streamlit


def streamlit_run():
    run_streamlit()

if __name__ == '__main__':
    streamlit_run()