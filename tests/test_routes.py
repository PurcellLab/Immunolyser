"""Smoke tests for Flask routes using the test client."""
import json
import pytest


def test_initialiser_loads(client):
    resp = client.get("/initialiser")
    assert resp.status_code == 200
    assert b"Sample name" in resp.data


def test_privacy_page_loads(client):
    resp = client.get("/privacy")
    assert resp.status_code == 200


def test_check_status_invalid_uuid(client):
    resp = client.get("/check_status/not-a-valid-uuid")
    assert resp.status_code == 404


def test_check_status_unknown_valid_uuid(client):
    resp = client.get("/check_status/aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee")
    assert resp.status_code == 200
    data = json.loads(resp.data)
    assert "status" in data


def test_download_core_missing_job_returns_404(client):
    resp = client.get(
        "/download_gibbscluster_core"
        "/aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee"
        "/SampleA/rep1/1of3"
    )
    assert resp.status_code == 404
