import streamlit as st
from ligify.streamlit_app import run_streamlit, run_ligify
from celery_worker import run_ligify_task, app
import time

def streamlit_run():
    chem, results, prog, chemical, filters = run_streamlit()

    if st.session_state.SUBMITTED:
        # Start Celery task
        task = run_ligify_task.delay(chemical, filters)
        st.session_state.task_id = task.id
        st.session_state.SUBMITTED = False
        st.experimental_rerun()  # Trigger a rerun to avoid immediate execution

    if 'task_id' in st.session_state and st.session_state.task_id is not None:
        task_id = st.session_state.task_id
        task = app.AsyncResult(task_id)

        if task.state == 'SUCCESS':
            st.session_state.result = task.result
            st.session_state.task_id = None
            st.experimental_rerun()  # Force rerun to update UI
        elif task.state == 'PENDING':
            st.write('Task is in the queue, please wait...')
            time.sleep(1)  # Sleep for a second to avoid too many reruns
            st.experimental_rerun()  # Check again
        elif task.state == 'STARTED':
            st.write('Task is currently being processed...')
            time.sleep(1)  # Sleep for a second to avoid too many reruns
            st.experimental_rerun()  # Check again
        elif task.state == 'FAILURE':
            st.write('Task failed, please try again.')
            st.session_state.task_id = None

    if 'result' in st.session_state:
        regulators = st.session_state.result['regulators']
        metrics = st.session_state.result['metrics']

        # Update the UI with the results
        run_ligify(chem, results, prog, chemical, filters, regulators, metrics)

streamlit_run()