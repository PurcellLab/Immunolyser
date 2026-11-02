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
    # Redis re-delivers tasks whose visibility timeout is exceeded. Must stay comfortably
    # above the longest possible task runtime (see LONG_JOB_TIME_LIMIT below) — otherwise
    # Redis can re-deliver a still-running long job to a second worker mid-run.
    CELERY_BROKER_TRANSPORT_OPTIONS = {'visibility_timeout': 100000}  # ~27.8 hours
    DEBUG = os.environ.get('DEBUG', 'False').lower() == 'true'


    BASE_URL = os.environ.get('BASE_URL', 'https://immunolyser.erc.monash.edu')

    # Days before job data is deleted from the server
    DATA_RETENTION_DAYS = int(os.environ.get('DATA_RETENTION_DAYS', 30))

    # Job input limites saved by variable. Used by both server and the client.
    SAMPLE_NAME_MAX_LENGTH = 30
    MAX_SAMPLES = 10
    MAX_TOTAL_PEPTIDES = 300000000000000000000000
    MAX_ALLELES = 6

    # Non-blocking heads-up shown per-sample above this peptide count (GibbsCluster
    # runtime scales worse than linearly with peptide count — see the ~9,840-peptide/
    # ~3h vs ~27,900-peptide/~9h23m comparison that motivated this). Informational
    # only, does not block submission; separate from the MAX_TOTAL_PEPTIDES hard cap
    # above, which sums across all samples rather than checking any one sample.
    PEPTIDE_WARNING_THRESHOLD = int(os.environ.get('PEPTIDE_WARNING_THRESHOLD', 10000))

    # reCAPTCHA v3 (leave blank to disable — useful for local dev)
    RECAPTCHA_SITE_KEY = os.environ.get('RECAPTCHA_SITE_KEY', '')
    RECAPTCHA_SECRET_KEY = os.environ.get('RECAPTCHA_SECRET_KEY', '')

    # Emails allowed to run jobs with an extended time limit (see LONG_JOB_* below),
    # for datasets too large to finish within the default limit. Comma-separated,
    # matched case-insensitively. Kept out of git — set only in the server's .env.
    # Deliberately not exposed anywhere in the UI or API responses.
    LONG_JOB_ALLOWED_EMAILS = {
        e.strip().lower()
        for e in os.environ.get('LONG_JOB_ALLOWED_EMAILS', '').split(',')
        if e.strip()
    }
    LONG_JOB_SOFT_TIME_LIMIT = 72000  # 20 hours
    LONG_JOB_TIME_LIMIT = 72300       # 20 hours 5 minutes

    # How long dispatch_pending_jobs waits after submission before dispatching a
    # job with the default (non-extended) time limit, giving the submitter a
    # window to add an allowlisted email via the confirmation page first. See
    # /initialiser's PENDING_JOBS_DIR comment in app/routes.py for why dispatch
    # is deferred at all rather than immediate.
    PENDING_DISPATCH_GRACE_SECONDS = int(os.environ.get('PENDING_DISPATCH_GRACE_SECONDS', 90))

    # --- MoDec (MHC class II motif deconvolution; alternative to GibbsCluster) ---
    # Binary-only, academic non-commercial licence (LICR). It is fetched at
    # `docker build` time, never committed — see app/tools/README.md. The option
    # is offered on the submit form only when this file exists, so a local install
    # built without it keeps working.
    MODEC_BIN = os.environ.get(
        'MODEC_BIN',
        os.path.join(os.path.dirname(os.path.abspath(__file__)), 'app', 'tools', 'MoDec', 'MoDec_unix'),
    )
    # Version of the licence text the submitter accepted; stored with the job so a
    # later wording change is distinguishable from an earlier acceptance.
    MODEC_LICENCE_VERSION = '1.2'
    # Motif counts searched. K above 6 is not offered by Gibbs here either, and the
    # report's K-switcher is sized for 6.
    MODEC_KMIN = int(os.environ.get('MODEC_KMIN', 1))
    MODEC_KMAX = int(os.environ.get('MODEC_KMAX', 6))
    # Random restarts per K. PROVISIONAL: MoDec's default is 20, and --MHC2 (which we
    # always pass) multiplies that by 5 via --specInits. A 5,987-peptide replicate took
    # 2m22s at nruns=2 on a laptop (K<=5), so 20 would be tens of minutes per replicate
    # and blow the 2h default job limit. To be re-set from a Linux benchmark.
    MODEC_NRUNS = int(os.environ.get('MODEC_NRUNS', 3))
    # Per-replicate wall-clock cap. On expiry the whole MoDec process group is killed
    # (a Celery time limit does NOT reap grandchildren — see the orphaned GibbsCluster
    # incident, ops/ORPHANED_GIBBSCLUSTER_2026-11-02 in the private notes).
    MODEC_TIME_LIMIT_SECONDS = int(os.environ.get('MODEC_TIME_LIMIT_SECONDS', 3000))