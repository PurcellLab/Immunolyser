# app/email_registry.py

import sqlite3
import os

project_root = os.path.dirname(os.path.realpath(os.path.join(__file__, "..")))

DB_PATH = os.path.join(project_root, 'results.sqlite')

def init_email_registry():
    """Create the email_registry table if it doesn't exist."""
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS email_registry (
                job_id TEXT PRIMARY KEY,
                email TEXT
            )
        ''')
        conn.commit()

def save_email(job_id, email):
    """Save or update an email associated with a job ID."""
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            REPLACE INTO email_registry (job_id, email)
            VALUES (?, ?)
        ''', (job_id, email))
        conn.commit()

def get_email(job_id):
    """Retrieve the email associated with a job ID."""
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute('SELECT email FROM email_registry WHERE job_id = ?', (job_id,))
        row = cursor.fetchone()
        return row[0] if row else None
