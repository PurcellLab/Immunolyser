# app/job_registry.py

import sqlite3
import os
from datetime import datetime

project_root = os.path.dirname(os.path.realpath(os.path.join(__file__, "..")))
DB_PATH = os.path.join(project_root, 'results.sqlite')

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
                error_message TEXT
            )
        ''')
        conn.commit()

def insert_job(job_id, ip_address, mhc_class, species, alleles,
               user_agent=None, referrer=None, status="PENDING"):
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            INSERT OR REPLACE INTO job_registry
            (job_id, ip_address, mhc_class, species, alleles, user_agent, referrer, status, submission_time)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            job_id, ip_address, mhc_class, species, alleles,
            user_agent, referrer, status, datetime.utcnow().isoformat()
        ))
        conn.commit()

def update_job_status(job_id, status, error_message=None):
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
