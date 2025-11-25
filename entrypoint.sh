#!/bin/bash
set -e

# Activate the virtual environment
source lenv/bin/activate

# Set Flask environment variables
export FLASK_APP=firstdemo.py
export FLASK_ENV=development
export IMMUNOLYSER_DATA=${IMMUNOLYSER_DATA}

# Ensure SQLite database exists and tables are created
DB_FILE="${IMMUNOLYSER_DATA}/results.sqlite"

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

# Run the appropriate process
if [ "$1" = "flask" ]; then
    flask run --host=0.0.0.0 --port=5000
elif [ "$1" = "celery" ]; then
    celery -A app.celery worker --loglevel=info --concurrency=1
else
    exec "$@"
fi
