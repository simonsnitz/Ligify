import streamlit as st
import pandas as pd
import numpy as np
from accID2operon import acc2operon


st.markdown("# Operon")


st.write("""
### Display the local genetic context of a protein
***
""")

st.header("Input protein accession ID")


with st.form("my_form"):
    seq = st.text_area("Accession ID", "WP_009944749.1", height= 50)
    checkbox_val = st.checkbox("Form checkbox")

    # Every form must have a submit button.
    submitted = st.form_submit_button("Submit")
    if submitted:
        data = acc2operon(seq)
        st.write(data)
