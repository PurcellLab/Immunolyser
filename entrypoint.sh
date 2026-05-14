#!/bin/bash
set -e

# Activate virtual environment
source lenv/bin/activate

# DB path
DB_FILE="${IMMUNOLYSER_DATA}/results.sqlite"

# Only Flask initializes the database
if [ "$1" = "flask" ]; then
    if [ ! -f "$DB_FILE" ]; then
        echo "Creating new SQLite database at $DB_FILE"
        sqlite3 "$DB_FILE" <<SQL
CREATE TABLE IF NOT EXISTS email_registry (
    job_id TEXT PRIMARY KEY,
    email TEXT,
    job_name TEXT
);

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
);
SQL
    else
        echo "Database already exists at $DB_FILE, skipping creation"
    fi
fi

# Run the service
if [ "$1" = "flask" ]; then
    gunicorn --workers 2 --bind 0.0.0.0:5000 --timeout 120 firstdemo:app
elif [ "$1" = "celery" ]; then
    celery -A app.celery worker --loglevel=info --concurrency=1
else
    exec "$@"
fi
