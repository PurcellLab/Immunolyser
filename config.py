import os

# Use the environment variable as the DB path
DB_PATH = os.environ.get('IMMUNOLYSER_DATA')
if not DB_PATH:
    raise RuntimeError("IMMUNOLYSER_DATA environment variable is not set!")

# If DB_PATH is a folder, append the database filename
if os.path.isdir(DB_PATH):
    DB_PATH = os.path.join(DB_PATH, 'results.sqlite')
class Config(object):
    SECRET_KEY = os.environ.get('SECRET_KEY') or b'6\xe9\xda\xead\x81\xf7\x8d\xbbH\x87\xe8m\xdd3%'

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
    DEBUG = True
    PIN = '123'

    # Job input limites saved by variable. Used by both server and the client.
    SAMPLE_NAME_MAX_LENGTH = 30
    MAX_SAMPLES = 10
    MAX_TOTAL_PEPTIDES = 300000000000000000000000
    MAX_ALLELES = 6