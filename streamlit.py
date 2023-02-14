import streamlit as st
import pandas as pd
import numpy as np
from accID2operon import acc2operon


st.title('GroovIO')

st.write("""
### Predict a transcription factor's DNA-binding sequence
***
""")

st.header("Input regulator accession ID")


with st.form("my_form"):
    seq = st.text_area("Accession ID", "WP_000", height= 50)
    checkbox_val = st.checkbox("Form checkbox")

    # Every form must have a submit button.
    submitted = st.form_submit_button("Submit")
    if submitted:
        data = acc2operon(seq)
        #protein = accID2sequence(seq)
        st.write(data)

st.write("Outside the form")