# app/job_registry.py

import sqlite3
import os
from datetime import datetime

# Use the environment variable as the DB path
DB_PATH = os.environ.get('IMMUNOLYSER_DATA')
if not DB_PATH:
    raise RuntimeError("IMMUNOLYSER_DATA environment variable is not set!")

# If DB_PATH is a folder, append the database filename
if os.path.isdir(DB_PATH):
    DB_PATH = os.path.join(DB_PATH, 'results.sqlite')

def init_job_registry():
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS job_registry (
                job_id TEXT PRIMARY KEY,
                ip_address TEXT,
                user_agent TEXT,
                mhc_class TEXT,
                species TEXT,
                alleles TEXT,
                status TEXT,
                submission_time TEXT,
                completed_time TEXT,
                error_message TEXT,
                country TEXT,
                referrer TEXT
            )
        ''')
        # Columns added after the table first shipped. Existing databases migrate in
        # place; NULL clustering_method means a job from before MoDec existed, i.e.
        # GibbsCluster. The licence columns record that the submitter accepted
        # MoDec's academic-non-commercial terms (the licence requires every user of
        # a server running it to have agreed) and which wording they saw.
        for column in ('clustering_method', 'modec_licence_accepted_at', 'modec_licence_version'):
            try:
                cursor.execute(f'ALTER TABLE job_registry ADD COLUMN {column} TEXT')
            except sqlite3.OperationalError:
                pass  # already present
        conn.commit()

def insert_job(job_id, country, mhc_class, species, alleles,
               user_agent=None, referrer=None, status="PENDING",
               clustering_method=None, modec_licence_version=None):
    now = datetime.utcnow().isoformat()
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            INSERT OR REPLACE INTO job_registry
            (job_id, country, mhc_class, species, alleles, user_agent, referrer, status, submission_time,
             clustering_method, modec_licence_accepted_at, modec_licence_version)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            job_id, country, mhc_class, species, alleles,
            user_agent, referrer, status, now,
            clustering_method,
            now if modec_licence_version else None,
            modec_licence_version,
        ))
        conn.commit()

def get_job_clustering_method(job_id):
    """Return the stored clustering method for a job, or None (pre-MoDec job)."""
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute('SELECT clustering_method FROM job_registry WHERE job_id = ?', (job_id,))
        row = cursor.fetchone()
        return row[0] if row else None

def get_completed_time(job_id):
    """Return the completed_time ISO string for a job, or None if not found."""
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute('SELECT completed_time FROM job_registry WHERE job_id = ?', (job_id,))
        row = cursor.fetchone()
        return row[0] if row else None

def get_jobs_older_than(cutoff_time):
    """Return job_ids with completed_time older than cutoff_time (ISO string) and not yet EXPIRED."""
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            SELECT job_id FROM job_registry
            WHERE status IN ('SUCCESS', 'FAILURE')
            AND completed_time <= ?
        ''', (cutoff_time,))
        return [row[0] for row in cursor.fetchall()]

def get_job_error(job_id):
    """Return the stored error_message for a job, or None if not found."""
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute('SELECT error_message FROM job_registry WHERE job_id = ?', (job_id,))
        row = cursor.fetchone()
        return row[0] if row else None

def update_job_status(job_id, status, error_message=None, logger=None):
    try:
        if logger:
            logger.info(f"Updating job status for job_id={job_id} to '{status}'.")
        with sqlite3.connect(DB_PATH) as conn:
            cursor = conn.cursor()
            cursor.execute('''
                UPDATE job_registry
                SET status = ?, completed_time = ?, error_message = ?
                WHERE job_id = ?
            ''', (
                status, datetime.utcnow().isoformat(), error_message, job_id
            ))
            conn.commit()
        if logger:
            logger.info(f"Job status updated successfully for job_id={job_id}.")
    except Exception as e:
        if logger:
            logger.error(f"Failed to update job status for job_id={job_id}: {e}", exc_info=True)
        raise
