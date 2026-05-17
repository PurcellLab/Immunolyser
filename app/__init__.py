from flask import Flask, render_template
from config import Config
from celery import Celery
from celery.schedules import crontab
import os

def make_celery(app):
    celery = Celery(
        app.import_name,
        backend=app.config['CELERY_RESULT_BACKEND'],
        broker=app.config['CELERY_BROKER_URL']
    )
    celery.conf.update(app.config)
    
    class ContextTask(celery.Task):
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return self.run(*args, **kwargs)
    celery.Task = ContextTask
    return celery

app = Flask(__name__)

@app.route("/")
def index():
    return render_template("index.html", index=True)

app.config.from_object(Config)

celery = make_celery(app)

from app import routes

celery.conf.CELERYBEAT_SCHEDULE = {
    'warn-expiring-jobs': {
        'task': 'app.routes.warn_expiring_jobs',
        'schedule': crontab(hour=8, minute=0),
    },
    'cleanup-expired-jobs': {
        'task': 'app.routes.cleanup_expired_jobs',
        'schedule': crontab(hour=9, minute=0),
    },
}

@app.after_request
def security_headers(response):
    response.headers['X-Frame-Options'] = 'DENY'
    response.headers['X-Content-Type-Options'] = 'nosniff'
    response.headers['X-XSS-Protection'] = '1; mode=block'
    return response