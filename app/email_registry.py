# app/email_registry.py

import sqlite3
import os

# Use the environment variable as the DB path
DB_PATH = os.environ.get('IMMUNOLYSER_DATA')
if not DB_PATH:
    raise RuntimeError("IMMUNOLYSER_DATA environment variable is not set!")

# If DB_PATH is a folder, append the database filename
if os.path.isdir(DB_PATH):
    DB_PATH = os.path.join(DB_PATH, 'results.sqlite')

def init_email_registry():
    """Create the email_registry table if it doesn't exist."""
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS email_registry (
                job_id TEXT PRIMARY KEY,
                email TEXT,
                job_name TEXT,
                email_sent INTEGER DEFAULT 0
            )
        ''')
        # Migrate existing databases that pre-date the email_sent column
        try:
            cursor.execute('ALTER TABLE email_registry ADD COLUMN email_sent INTEGER DEFAULT 0')
        except sqlite3.OperationalError:
            pass
        try:
            cursor.execute('ALTER TABLE email_registry ADD COLUMN warning_sent INTEGER DEFAULT 0')
        except sqlite3.OperationalError:
            pass
        conn.commit()

def save_email(job_id, email, job_name=None):
    """Save or update an email and optional job name associated with a job ID."""
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            REPLACE INTO email_registry (job_id, email, job_name)
            VALUES (?, ?, ?)
        ''', (job_id, email, job_name))
        conn.commit()

def get_email(job_id):
    """Retrieve the email associated with a job ID."""
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute('SELECT email FROM email_registry WHERE job_id = ?', (job_id,))
        row = cursor.fetchone()
        return row[0] if row else None

def get_job_name(job_id):
    """Retrieve the job name associated with a job ID."""
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute('SELECT job_name FROM email_registry WHERE job_id = ?', (job_id,))
        row = cursor.fetchone()
        return row[0] if row else None

def get_jobs_due_for_warning(cutoff_time):
    """Return (job_id, email, job_name) rows where warning not yet sent and email exists."""
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            SELECT job_id, email, job_name FROM email_registry
            WHERE warning_sent = 0 AND email IS NOT NULL AND job_id IN (
                SELECT job_id FROM job_registry
                WHERE status IN ('SUCCESS', 'FAILURE')
                AND completed_time <= ?
            )
        ''', (cutoff_time,))
        return cursor.fetchall()

def claim_warning_send(job_id):
    """Atomically mark the warning email as sent. Returns True if this caller won the race."""
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute(
            'UPDATE email_registry SET warning_sent = 1 WHERE job_id = ? AND warning_sent = 0',
            (job_id,)
        )
        conn.commit()
        return cursor.rowcount == 1

def claim_email_send(job_id):
    """Atomically mark the email as sent. Returns True if this caller won the race, False if already sent."""
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute(
            'UPDATE email_registry SET email_sent = 1 WHERE job_id = ? AND email_sent = 0',
            (job_id,)
        )
        conn.commit()
        return cursor.rowcount == 1