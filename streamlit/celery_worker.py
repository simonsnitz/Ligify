# streamlit/celery_worker.py
import time
from celery import Celery
from ligify.fetch_data import fetch_data

app = Celery('celery_worker', broker='redis://redis:6379/0', backend='redis://redis:6379/0')

app.conf.update(
    result_backend='redis://redis:6379/0',
    task_serializer='json',
    accept_content=['json'],
    result_serializer='json',
    timezone='UTC',
    enable_utc=True,
    worker_concurrency=1,  # Ensure only one task is processed at a time
)

@app.task
def run_ligify_task(chemical, filters):
    try:
        # self.update_state(state="STARTED")
        # Fetch data
        regulators, metrics = fetch_data(chemical["InChiKey"], filters)

        # self.update_state(state="SUCCESS")
        # time.sleep(1000) # Display the message
        
        return {"regulators": regulators, "metrics": metrics}
    except Exception as e:
        print(e)
        # Optionally, you can store the error details in the result backend
        # self.update_state(
        #     state='FAILURE',
        #     meta={'exc_type': type(e).__name__, 'exc_message': str(e)}
        # )
