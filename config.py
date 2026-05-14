import os

# Use the environment variable as the DB path
DB_PATH = os.environ.get('IMMUNOLYSER_DATA')
if not DB_PATH:
    raise RuntimeError("IMMUNOLYSER_DATA environment variable is not set!")

# If DB_PATH is a folder, append the database filename
if os.path.isdir(DB_PATH):
    DB_PATH = os.path.join(DB_PATH, 'results.sqlite')
class Config(object):
    SECRET_KEY = os.environ.get('SECRET_KEY')
    if not SECRET_KEY:
        raise RuntimeError("SECRET_KEY environment variable must be set!")

    # Location to store all the data
    IMMUNOLYSER_DATA = os.environ.get('IMMUNOLYSER_DATA')

    # Task id of demo job
    DEMO_TASK_ID = os.environ.get('DEMO_TASK_ID')

    # Dynamically set the CELERY_BROKER_URL based on the environment
    if os.environ.get('IS_DOCKER') == 'true':
        CELERY_BROKER_URL = 'redis://redis:6379/0'  # Use redis container name in Docker
    else:
        CELERY_BROKER_URL = 'redis://localhost:6379/0'  # Use localhost on the server
    CELERY_RESULT_BACKEND = f'db+sqlite:///{DB_PATH}'
    CELERY_DEFAULT_QUEUE='celery'  # Ensure all tasks are routed to 'celery' queue
    # Redis re-delivers tasks whose visibility timeout is exceeded. Long jobs (many alleles)
    # can run for several hours, so set this well above the worst-case job duration.
    CELERY_BROKER_TRANSPORT_OPTIONS = {'visibility_timeout': 43200}  # 12 hours
    DEBUG = os.environ.get('DEBUG', 'False').lower() == 'true'


    BASE_URL = os.environ.get('BASE_URL', 'https://immunolyser.erc.monash.edu')

    # Days before job data is deleted from the server
    DATA_RETENTION_DAYS = int(os.environ.get('DATA_RETENTION_DAYS', 30))

    # Job input limites saved by variable. Used by both server and the client.
    SAMPLE_NAME_MAX_LENGTH = 30
    MAX_SAMPLES = 10
    MAX_TOTAL_PEPTIDES = 300000000000000000000000
    MAX_ALLELES = 6