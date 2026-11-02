"""MoDec in the report: the shared section partial (live page and offline export) and
the per-motif peptide download route."""
import os
import uuid

import pytest

from app import app as flask_app, routes, utils
from app import modec
from test_modec import make_output_tree

TASK = str(uuid.uuid4())


@pytest.fixture
def job(tmp_path, monkeypatch):
    """A finished MoDec replicate under <tmp>/app/static/images/<task>/S1/modec/repA
    plus the user's original upload, with the app pointed at the tmp tree."""
    root = tmp_path / 'proj'
    images = root / 'app' / 'static' / 'images'
    rep = images / TASK / 'S1' / modec.MODEC_DIRNAME / 'repA'
    make_output_tree(str(rep), kmax=3)
    (rep / modec.LOGO_DIRNAME).mkdir()
    for k in (1, 2, 3):
        for mid in [str(i) for i in range(1, k + 1)] + ['flat']:
            (rep / modec.LOGO_DIRNAME / f'{modec.logo_basename(mid, k)}-001.jpg').write_bytes(b'x')

    data = tmp_path / 'data'
    (data / TASK / 'S1').mkdir(parents=True)
    # peptides PEPTIDE000XYZ.. match the fixture's bestPepResp rows
    rows = ['Peptide,Intensity']
    rows += [f'PEPTIDE{i:03d}XYZ,{1000 + i}' for i in range(40)]
    (data / TASK / 'S1' / 'repA.csv').write_text('\n'.join(rows) + '\n')

    monkeypatch.setattr(routes, 'project_root', str(root))
    monkeypatch.setattr(utils, 'data_mount', str(data))
    return {'root': str(root), 'images': str(images), 'rep': str(rep)}


def _render(job, **extra):
    data = modec.report_data(job['images'], TASK, {'S1': ['repA']})
    ctx = dict(modecData=data, taskId=TASK)
    ctx.update(extra)
    with flask_app.test_request_context():
        from flask import render_template
        return render_template('includes/modec_section.html', **ctx)


class TestSectionPartial:
    def test_live_page_has_logos_counts_and_download_links(self, job):
        html = _render(job)
        assert 'Peptide Motif Deconvolution (MoDec)' in html
        assert '/static/images/' in html
        assert 'modec_logos_2of2-001.jpg' in html
        assert f'/download_modec_core/{TASK}/S1/repA/2/1' in html
        assert 'download peptides' in html

    def test_flat_motif_is_labelled_and_last(self, job):
        html = _render(job)
        assert 'Flat motif' in html
        assert html.index('Motif 1') < html.index('Flat motif')

    def test_only_the_lowest_aic_k_is_visible_initially(self, job):
        html = _render(job)
        # fixture: lowest AIC at K=2 -> K=1 and K=3 blocks start hidden, K=2 does not
        import re
        blocks = re.findall(r'<div class="modec-k[^"]*" data-k="(\d)"\s*([^>]*)>', html)
        hidden = {k for k, attrs in blocks if 'display: none' in attrs}
        assert hidden == {'1', '3'}
        assert 'selected' in html and '(lowest AIC)' in html

    def test_graph_shown_beside_logos_when_available(self, job):
        plot = os.path.join(job['rep'], modec.AIC_PLOT_FILENAME)
        with open(plot, 'wb') as f:
            f.write(b'x')
        html = _render(job)
        assert 'alt="AIC by number of motifs"' in html
        assert 'aic_plot.png' in html
        assert 'K=2: AIC' in html          # values are in the tooltip
        assert '<th>AIC</th>' not in html  # table is only the fallback

    def test_table_fallback_for_jobs_without_a_graph(self, job):
        html = _render(job)
        assert 'alt="AIC by number of motifs"' not in html
        assert '<th>AIC</th>' in html

    def test_attribution_and_licence_credit_present(self, job):
        html = _render(job)
        assert 'Racle' in html
        assert 'Ludwig Institute for Cancer Research' in html

    def test_offline_mode_inlines_images_and_omits_downloads(self, job):
        image_map = {'/static/' + m['image']: 'data:image/jpeg;base64,AAAA'
                     for m in modec.report_data(job['images'], TASK, {'S1': ['repA']})['S1']['repA']['k'][2]
                     if m['image']}
        html = _render(job, image_map=image_map)
        assert 'data:image/jpeg;base64,AAAA' in html
        assert '/download_modec_core/' not in html
        assert 'download peptides' not in html

    def test_empty_data_says_so(self, job):
        with flask_app.test_request_context():
            from flask import render_template
            html = render_template('includes/modec_section.html', modecData={}, taskId=TASK)
        assert 'No MoDec results were produced' in html


class TestDownloadRoute:
    def _get(self, client, k, motif, sample='S1', replicate='repA', task=TASK):
        return client.get(f'/download_modec_core/{task}/{sample}/{replicate}/{k}/{motif}')

    def test_returns_enriched_csv_for_a_motif(self, client, job):
        resp = self._get(client, 2, 1)
        assert resp.status_code == 200
        assert resp.headers['Content-Type'].startswith('text/csv')
        assert 'repA_modec_K2_motif1.csv' in resp.headers['Content-Disposition']
        lines = resp.data.decode().strip().split('\n')
        header = lines[0].split(',')
        assert {'StrippedPeptide', 'Core', 'Core start (P1 offset)', 'Responsibility', 'Motif', 'K'} <= set(header)
        # enriched with the user's original upload column
        assert 'Intensity' in header
        assert len(lines) - 1 == 10  # K=2 motif 1 has 10 peptides in the fixture

    def test_flat_motif_downloadable(self, client, job):
        resp = self._get(client, 2, 'flat')
        assert resp.status_code == 200
        assert len(resp.data.decode().strip().split('\n')) - 1 == 2

    def test_unknown_k_is_404(self, client, job):
        assert self._get(client, 9, 1).status_code == 404

    @pytest.mark.parametrize('k,motif', [('0', '1'), ('x', '1'), ('2', 'abc'), ('2', '../1'), ('100', '1'), ('2', '0')])
    def test_malformed_k_or_motif_rejected(self, client, job, k, motif):
        assert self._get(client, k, motif).status_code in (400, 404)

    def test_invalid_task_id(self, client, job):
        assert self._get(client, 2, 1, task='not-a-uuid').status_code == 400

    @pytest.mark.parametrize('sample,replicate', [('..', 'repA'), ('S1', '..'), ('S1/x', 'repA')])
    def test_path_traversal_rejected(self, client, job, sample, replicate):
        # '..' passes the character check (dots are allowed) and is stopped by the
        # real-path containment check with 403; slashes never match the route (404).
        assert self._get(client, 2, 1, sample=sample, replicate=replicate).status_code in (400, 403, 404)

    def test_missing_original_upload_falls_back_to_unenriched(self, client, job):
        os.remove(os.path.join(utils.data_mount, TASK, 'S1', 'repA.csv'))
        resp = self._get(client, 2, 1)
        assert resp.status_code == 200
        assert 'Intensity' not in resp.data.decode().split('\n')[0]
