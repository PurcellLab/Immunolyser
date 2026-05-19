from app import app
from app.email_registry import init_email_registry
from app.job_registry import init_job_registry

# Initialize both tables
init_email_registry()
init_job_registry()
