import os
import tempfile
import pytest

# Set required env vars before any app imports.
# On the server the developer's shell already has these; this fallback
# lets the test suite run in a plain "pytest" invocation from the repo root.
_test_data_dir = tempfile.mkdtemp(prefix="immunolyser_test_")
os.environ.setdefault("IMMUNOLYSER_DATA", _test_data_dir)
os.environ.setdefault("SECRET_KEY", "test-secret-key-not-for-production")
os.environ.setdefault("DEMO_TASK_ID", "00000000-0000-0000-0000-000000000000")

from app.job_registry import init_job_registry
init_job_registry()

from app import app as _flask_app
_flask_app.config["TESTING"] = True


@pytest.fixture(scope="session")
def app():
    yield _flask_app


@pytest.fixture
def client(app):
    return app.test_client()
