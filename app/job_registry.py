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
        conn.commit()

def insert_job(job_id, country, mhc_class, species, alleles,
               user_agent=None, referrer=None, status="PENDING"):
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            INSERT OR REPLACE INTO job_registry
            (job_id, country, mhc_class, species, alleles, user_agent, referrer, status, submission_time)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            job_id, country, mhc_class, species, alleles,
            user_agent, referrer, status, datetime.utcnow().isoformat()
        ))
        conn.commit()

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
