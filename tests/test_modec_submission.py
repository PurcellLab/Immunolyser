"""POST /initialiser: clustering-method validation and MoDec licence acceptance.

The licence (LICR, academic non-commercial) requires every user of a server running
MoDec to have agreed to its terms, so the server-side check is the one that counts —
these tests post to the route directly, bypassing the browser tick box.
"""
import io
import json
import os
import sqlite3

import pytest

from constants import ClusteringMethod
from app import routes
from app.job_registry import DB_PATH, get_job_clustering_method


@pytest.fixture(autouse=True)
def _open_submission(app, monkeypatch):
    """Skip reCAPTCHA (its DEBUG path) and CSRF, which the test client has no token for."""
    monkeypatch.setitem(app.config, 'DEBUG', True)
    monkeypatch.setitem(app.config, 'RECAPTCHA_SECRET_KEY', '')
    monkeypatch.setitem(app.config, 'WTF_CSRF_ENABLED', False)


def _post(client, mhc_class, **extra):
    """Submit one sample. As in the real form, the sample name and its file input
    share a field name ('file-1'): one part carries the name, one the upload."""
    data = {
        'file-1': ['Sample1', (io.BytesIO(b'Peptide\nAAAAAAAAAAAAA\n'), 'pep.csv')],
        'mhc_class': mhc_class,
        'species': 'Human',
        'motif_length': '9',
        'alleles': '',
        'useFullDB': 'no',
    }
    data.update(extra)
    return client.post('/initialiser', data=data, content_type='multipart/form-data')


def _job_id(resp):
    assert resp.status_code == 302, resp.data
    return resp.headers['Location'].rstrip('/').split('/')[-1]


def _pending(job_id):
    path = os.path.join(routes.data_mount, routes.PENDING_JOBS_DIR, f'{job_id}.json')
    with open(path) as f:
        return json.load(f)


def _licence_columns(job_id):
    with sqlite3.connect(DB_PATH) as conn:
        return conn.execute(
            'SELECT modec_licence_version, modec_licence_accepted_at FROM job_registry WHERE job_id = ?',
            (job_id,)).fetchone()


class TestMethodValidation:
    def test_unknown_method_rejected(self, client):
        assert _post(client, 'II', clustering_method='bogus').status_code == 400

    def test_modec_for_class_one_rejected(self, client, monkeypatch):
        monkeypatch.setattr(routes, 'modec_available', lambda _bin: True)
        resp = _post(client, 'I', clustering_method=ClusteringMethod.MODEC, modec_licence_accepted='yes')
        assert resp.status_code == 400

    def test_modec_requires_installed_binary(self, client, monkeypatch):
        monkeypatch.setattr(routes, 'modec_available', lambda _bin: False)
        resp = _post(client, 'II', clustering_method=ClusteringMethod.MODEC, modec_licence_accepted='yes')
        assert resp.status_code == 400

    @pytest.mark.parametrize('tick', [None, '', 'no', 'on', 'true'])
    def test_modec_without_exact_acceptance_rejected(self, client, monkeypatch, tick):
        monkeypatch.setattr(routes, 'modec_available', lambda _bin: True)
        extra = {'clustering_method': ClusteringMethod.MODEC}
        if tick is not None:
            extra['modec_licence_accepted'] = tick
        assert _post(client, 'II', **extra).status_code == 400


class TestAcceptedSubmissions:
    def test_missing_method_defaults_to_gibbs_for_class_two(self, client):
        job_id = _job_id(_post(client, 'II'))
        assert _pending(job_id)['clustering_method'] == ClusteringMethod.GIBBS
        assert get_job_clustering_method(job_id) == ClusteringMethod.GIBBS

    def test_class_one_job_is_gibbs(self, client):
        job_id = _job_id(_post(client, 'I'))
        assert _pending(job_id)['clustering_method'] == ClusteringMethod.GIBBS

    def test_accepted_modec_job_is_queued_and_acceptance_recorded(self, client, monkeypatch):
        monkeypatch.setattr(routes, 'modec_available', lambda _bin: True)
        job_id = _job_id(_post(client, 'II', clustering_method=ClusteringMethod.MODEC,
                               modec_licence_accepted='yes'))
        assert _pending(job_id)['clustering_method'] == ClusteringMethod.MODEC
        assert get_job_clustering_method(job_id) == ClusteringMethod.MODEC
        version, accepted_at = _licence_columns(job_id)
        assert version == routes.app.config['MODEC_LICENCE_VERSION']
        assert accepted_at

    def test_gibbs_job_records_no_licence_acceptance(self, client):
        job_id = _job_id(_post(client, 'II'))
        assert _licence_columns(job_id) == (None, None)


class TestForm:
    def test_option_hidden_when_binary_missing(self, client, monkeypatch):
        monkeypatch.setattr(routes, 'modec_available', lambda _bin: False)
        assert b'name="clustering_method"' not in client.get('/initialiser').data

    def test_option_shown_when_binary_present(self, client, monkeypatch):
        monkeypatch.setattr(routes, 'modec_available', lambda _bin: True)
        html = client.get('/initialiser').data
        assert b'name="clustering_method"' in html
        assert b'modec_licence_accepted' in html
