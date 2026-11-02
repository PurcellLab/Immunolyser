"""Tests for MoDec support: method constants, command building, output parsing,
and the process-group-safe runner (app/modec.py).

Fixtures are synthetic files in MoDec's exact output format (checked against a
real run on the tool's own HLA-II test set) so no third-party output is committed.
"""
import math
import os
import stat
import sys
import time
import types

import pytest

from constants import ClusteringMethod, MHC_Class
from app import modec

POSIX_ONLY = pytest.mark.skipif(os.name != 'posix', reason='process groups need POSIX')


# --------------------------------------------------------------------------
# Fixture builders
# --------------------------------------------------------------------------

def write_pwm(path, n_peptides, positions=9, hot=None):
    """A MoDec PWM file; every column sums to 1. `hot` = {(aa, pos): value} overrides."""
    aas = modec.MODEC_AA_ORDER
    cols = [[1.0 / len(aas)] * len(aas) for _ in range(positions)]
    for (aa, pos), val in (hot or {}).items():
        col = cols[pos]
        col[aas.index(aa)] = val
        rest = (1.0 - val) / (len(aas) - 1)
        for i in range(len(aas)):
            if i != aas.index(aa):
                col[i] = rest
    with open(path, 'w') as f:
        f.write(f'#\t{os.path.basename(path)[:-4]}_n_pep_{n_peptides}\t{n_peptides}\n')
        f.write('\t' + '\t'.join(str(p + 1) for p in range(positions)) + '\n')
        for i, aa in enumerate(aas):
            f.write(aa + '\t' + '\t'.join(f'{cols[p][i]:.6g}' for p in range(positions)) + '\n')


def write_resp(path, motifs):
    """bestPepResp file with one peptide per entry of `motifs` (the Motif_1st)."""
    with open(path, 'w') as f:
        f.write('Peptide\tLength\tScore\tPep_weight\tresp_1st\tMotif_1st\tP1_1st\tCore_1st\t'
                'resp_2nd\tMotif_2nd\tP1_2nd\tCore_2nd\n')
        for i, m in enumerate(motifs):
            f.write(f'PEPTIDE{i:03d}XYZ\t13\t9.5\t1\t0.98\t{m}\t2\tPEPTIDEXY\t0.01\t0\t2\tPEPTIDEXY\n')


def make_output_tree(root, kmax=3):
    """A MoDec result folder for K=1..kmax with consistent counts."""
    for sub in ('PWMs', 'Responsibilities', 'Summary'):
        os.makedirs(os.path.join(root, sub))
    aic_lines = ['n\tlogL\tlogL/n\tAIC']
    for k in range(1, kmax + 1):
        sizes = [10 * (i + 1) for i in range(k)]
        for i, size in enumerate(sizes, start=1):
            # real motifs always carry information (a uniform PWM is skipped by render_logos)
            write_pwm(os.path.join(root, 'PWMs', f'PWM_K{k}_{i}.txt'), size, hot={('F', 0): 0.5})
        # a real flat motif carries a little information (a compositional bias), unlike an
        # empty one; see test_empty_uniform_flat_motif_is_skipped_not_failed
        write_pwm(os.path.join(root, 'PWMs', f'PWM_K{k}_flat.txt'), 2, hot={('A', 0): 0.2})
        motifs = [i for i, size in enumerate(sizes, start=1) for _ in range(size)] + [0, 0]
        write_resp(os.path.join(root, 'Responsibilities', f'bestPepResp_K{k}.txt'), motifs)
        # AIC is lowest (best) at K=2
        aic_lines.append(f'{k}\t{100.0 * k}\t{50.0 * k}\t{-1000 + 40 * abs(k - 2)}')
    with open(os.path.join(root, 'Summary', f'nMotScores_K1-{kmax}.txt'), 'w') as f:
        f.write('\n'.join(aic_lines) + '\n')


def make_fake_binary(path, body):
    """An executable python script standing in for MoDec."""
    with open(path, 'w') as f:
        f.write(f'#!{sys.executable}\n{body}\n')
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)


# --------------------------------------------------------------------------
# ClusteringMethod
# --------------------------------------------------------------------------

class TestClusteringMethod:
    def test_default_is_gibbs(self):
        assert ClusteringMethod.DEFAULT == ClusteringMethod.GIBBS

    @pytest.mark.parametrize('raw', [None, '', '   '])
    def test_blank_normalises_to_default(self, raw):
        assert ClusteringMethod.normalise(raw) == ClusteringMethod.GIBBS

    def test_normalise_lowercases_and_does_not_validate(self):
        assert ClusteringMethod.normalise(' MoDec ') == ClusteringMethod.MODEC
        assert ClusteringMethod.normalise('bogus') == 'bogus'

    def test_modec_allowed_for_class_two_only(self):
        assert ClusteringMethod.is_allowed(ClusteringMethod.MODEC, MHC_Class.Two)
        assert not ClusteringMethod.is_allowed(ClusteringMethod.MODEC, MHC_Class.One)

    def test_gibbs_allowed_for_both_classes(self):
        assert ClusteringMethod.is_allowed(ClusteringMethod.GIBBS, MHC_Class.One)
        assert ClusteringMethod.is_allowed(ClusteringMethod.GIBBS, MHC_Class.Two)

    def test_unknown_method_or_class_rejected(self):
        assert not ClusteringMethod.is_allowed('bogus', MHC_Class.Two)
        assert not ClusteringMethod.is_allowed(ClusteringMethod.GIBBS, 'III')


# --------------------------------------------------------------------------
# Command building
# --------------------------------------------------------------------------

class TestBuildCommand:
    def _cmd(self, **kw):
        args = dict(kmin=1, kmax=6, nruns=3, motif_length=9, species='Human')
        args.update(kw)
        return modec.build_command('/x/MoDec_unix', '/in/p.txt', '/out/p', **args)

    def test_core_flags(self):
        cmd = self._cmd()
        assert cmd[:5] == ['/x/MoDec_unix', '-i', '/in/p.txt', '-o', '/out/p']
        assert '--MHC2' in cmd and '--log' in cmd
        assert cmd[cmd.index('--Kmax') + 1] == '6'
        assert cmd[cmd.index('--nruns') + 1] == '3'
        assert cmd[cmd.index('-L') + 1] == '9'

    def test_human_keeps_default_background(self):
        assert '--theta_norm' not in self._cmd(species='Human')

    def test_mouse_uses_mouse_background(self):
        cmd = self._cmd(species='Mouse')
        assert cmd[cmd.index('--theta_norm') + 1] == '4'

    def test_unknown_species_leaves_default(self):
        assert '--theta_norm' not in self._cmd(species=None)

    def test_no_shell_string_all_args_are_str(self):
        assert all(isinstance(a, str) for a in self._cmd())

    @pytest.mark.parametrize('kmin,kmax', [(0, 5), (4, 3), (-1, 2)])
    def test_bad_k_range_rejected(self, kmin, kmax):
        with pytest.raises(ValueError):
            self._cmd(kmin=kmin, kmax=kmax)

    def test_zero_runs_rejected(self):
        with pytest.raises(ValueError):
            self._cmd(nruns=0)


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------

class TestParsePwm:
    def test_reads_counts_shape_and_columns_sum_to_one(self, tmp_path):
        p = tmp_path / 'PWM_K3_1.txt'
        write_pwm(p, 2045)
        pwm = modec.parse_pwm(p)
        assert pwm['n_peptides'] == 2045
        assert pwm['positions'] == 9
        assert set(pwm['freq']) == set(modec.MODEC_AA_ORDER)
        for pos in range(9):
            assert sum(pwm['freq'][aa][pos] for aa in pwm['freq']) == pytest.approx(1.0, abs=1e-4)

    def test_rejects_non_pwm_file(self, tmp_path):
        p = tmp_path / 'x.txt'
        p.write_text('hello\nworld\n')
        with pytest.raises(ValueError):
            modec.parse_pwm(p)

    def test_rejects_wrong_amino_acid_rows(self, tmp_path):
        p = tmp_path / 'PWM_K1_1.txt'
        p.write_text('#\tn\t5\n\t1\t2\nA\t0.5\t0.5\nC\t0.5\t0.5\n')
        with pytest.raises(ValueError):
            modec.parse_pwm(p)


class TestGibbsOrientation:
    def test_transposes_to_positions_by_gibbs_amino_acid_order(self, tmp_path):
        p = tmp_path / 'PWM_K1_1.txt'
        write_pwm(p, 10, hot={('F', 0): 0.6, ('L', 8): 0.5})
        rows = modec.pwm_as_gibbs_rows(modec.parse_pwm(p))
        assert len(rows) == 9 and all(len(r) == 20 for r in rows)
        assert rows[0][modec.GIBBS_AA_ORDER.index('F')] == pytest.approx(0.6, abs=1e-4)
        assert rows[8][modec.GIBBS_AA_ORDER.index('L')] == pytest.approx(0.5, abs=1e-4)
        # a row (position) is still a distribution
        assert sum(rows[3]) == pytest.approx(1.0, abs=1e-3)

    def test_orders_are_permutations_of_each_other(self):
        assert sorted(modec.MODEC_AA_ORDER) == sorted(modec.GIBBS_AA_ORDER)


class TestParseResponsibilitiesAndAic:
    def test_best_pep_resp_types(self, tmp_path):
        p = tmp_path / 'bestPepResp_K2.txt'
        write_resp(p, [1, 2, 0])
        rows = modec.parse_best_pep_resp(p)
        assert [r['Motif_1st'] for r in rows] == [1, 2, 0]
        assert isinstance(rows[0]['resp_1st'], float)
        assert isinstance(rows[0]['Length'], int)
        assert rows[0]['Core_1st'] == 'PEPTIDEXY'

    def test_aic_sorted_and_suggested_k_is_lowest(self, tmp_path):
        p = tmp_path / 'nMotScores_K1-3.txt'
        p.write_text('n\tlogL\tlogL/n\tAIC\n3\t1\t1\t-10443\n1\t1\t1\t-6954\n2\t1\t1\t-10505\n')
        rows = modec.parse_aic(p)
        assert [r['n'] for r in rows] == [1, 2, 3]
        assert modec.suggested_k(rows) == 2

    def test_suggested_k_none_when_empty(self):
        assert modec.suggested_k([]) is None


class TestCollectResults:
    def test_assembles_all_ks(self, tmp_path):
        make_output_tree(str(tmp_path / 'out'), kmax=3)
        res = modec.collect_results(str(tmp_path / 'out'))
        assert list(res['k']) == [1, 2, 3]
        assert res['suggested_k'] == 2
        assert set(res['k'][3]['motifs']) == {'1', '2', '3', 'flat'}

    def test_assignment_counts_match_pwm_headers(self, tmp_path):
        """MoDec's per-motif peptide count in the PWM header equals the number of
        peptides whose best motif it is (verified on real output)."""
        make_output_tree(str(tmp_path / 'out'), kmax=3)
        res = modec.collect_results(str(tmp_path / 'out'))
        for k, entry in res['k'].items():
            counts = {}
            for row in entry['assignments']:
                counts[str(row['Motif_1st'])] = counts.get(str(row['Motif_1st']), 0) + 1
            for motif_id, pwm in entry['motifs'].items():
                key = '0' if motif_id == 'flat' else motif_id
                assert counts[key] == pwm['n_peptides'], (k, motif_id)

    def test_missing_summary_raises(self, tmp_path):
        os.makedirs(tmp_path / 'out' / 'PWMs')
        os.makedirs(tmp_path / 'out' / 'Responsibilities')
        os.makedirs(tmp_path / 'out' / 'Summary')
        with pytest.raises(FileNotFoundError):
            modec.collect_results(str(tmp_path / 'out'))


# --------------------------------------------------------------------------
# CPU detection
# --------------------------------------------------------------------------

class TestAvailableCpus:
    def test_cgroup_v2_quota_caps_cpu_count(self, tmp_path):
        (tmp_path / 'cpu.max').write_text('800000 100000\n')
        assert modec._cgroup_cpu_quota(str(tmp_path)) == 8.0
        assert modec.available_cpus(str(tmp_path)) <= 8

    def test_cgroup_v2_unlimited(self, tmp_path):
        (tmp_path / 'cpu.max').write_text('max 100000\n')
        assert modec._cgroup_cpu_quota(str(tmp_path)) is None

    def test_cgroup_v1_quota(self, tmp_path):
        d = tmp_path / 'cpu'
        d.mkdir()
        (d / 'cpu.cfs_quota_us').write_text('400000\n')
        (d / 'cpu.cfs_period_us').write_text('100000\n')
        assert modec._cgroup_cpu_quota(str(tmp_path)) == 4.0

    def test_no_cgroup_files_means_no_quota(self, tmp_path):
        assert modec._cgroup_cpu_quota(str(tmp_path)) is None

    def test_fractional_quota_rounds_up_and_floor_is_one(self, tmp_path):
        (tmp_path / 'cpu.max').write_text('50000 100000\n')
        assert modec.available_cpus(str(tmp_path)) == 1

    def test_never_exceeds_host_count(self):
        assert 1 <= modec.available_cpus('/nonexistent') <= (os.cpu_count() or 1)


# --------------------------------------------------------------------------
# Runner
# --------------------------------------------------------------------------

@POSIX_ONLY
class TestRunModec:
    def test_success(self, tmp_path):
        out = tmp_path / 'img' / 'modec' / 'rep'
        modec.run_modec([sys.executable, '-c', 'pass'], str(out), 30)

    def test_nonzero_exit_raises_with_stderr(self, tmp_path):
        with pytest.raises(RuntimeError, match='boom'):
            modec.run_modec(
                [sys.executable, '-c', 'import sys; sys.stderr.write("boom"); sys.exit(3)'],
                str(tmp_path / 'o' / 'rep'), 30)

    def test_existing_output_folder_refused(self, tmp_path):
        out = tmp_path / 'rep'
        out.mkdir()
        with pytest.raises(FileExistsError):
            modec.run_modec([sys.executable, '-c', 'pass'], str(out), 30)

    def test_timeout_kills_the_whole_process_group(self, tmp_path):
        """The regression this guards: a timeout must also kill grandchildren, or
        they outlive the job (the orphaned-GibbsCluster incident)."""
        pidfile = tmp_path / 'child.pid'
        parent = (
            'import subprocess, sys, time\n'
            f'c = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])\n'
            f'open({str(pidfile)!r}, "w").write(str(c.pid))\n'
            'time.sleep(60)\n'
        )
        with pytest.raises(TimeoutError):
            modec.run_modec([sys.executable, '-c', parent], str(tmp_path / 'o' / 'rep'), 2)
        child_pid = int(pidfile.read_text())
        for _ in range(50):  # allow the kernel a moment to reap
            try:
                os.kill(child_pid, 0)
            except ProcessLookupError:
                break
            time.sleep(0.1)
        else:
            pytest.fail('grandchild survived the timeout')


class TestReadJobMethod:
    def test_missing_marker_means_gibbs(self, tmp_path):
        assert modec.read_job_method(str(tmp_path), 'job1') == ClusteringMethod.GIBBS

    def test_reads_marker(self, tmp_path):
        (tmp_path / 'job1').mkdir()
        (tmp_path / 'job1' / 'clustering_method.txt').write_text('modec\n')
        assert modec.read_job_method(str(tmp_path), 'job1') == ClusteringMethod.MODEC

    def test_blank_marker_means_gibbs(self, tmp_path):
        (tmp_path / 'job1').mkdir()
        (tmp_path / 'job1' / 'clustering_method.txt').write_text('  \n')
        assert modec.read_job_method(str(tmp_path), 'job1') == ClusteringMethod.GIBBS


class TestReportData:
    def _tree(self, tmp_path, with_logos=True):
        images = tmp_path / 'static' / 'images'
        rep = images / 'job1' / 'S1' / modec.MODEC_DIRNAME / 'repA'
        make_output_tree(str(rep), kmax=3)
        if with_logos:
            (rep / modec.LOGO_DIRNAME).mkdir()
            for k in (1, 2, 3):
                for mid in [str(i) for i in range(1, k + 1)] + ['flat']:
                    (rep / modec.LOGO_DIRNAME / f'{modec.logo_basename(mid, k)}-001.jpg').write_bytes(b'x')
        return str(images)

    def test_structure_and_ordering(self, tmp_path):
        images = self._tree(tmp_path)
        data = modec.report_data(images, 'job1', {'S1': ['repA']})
        rep = data['S1']['repA']
        assert rep['suggested_k'] == 2
        assert list(rep['k']) == [1, 2, 3]
        assert [m['id'] for m in rep['k'][3]] == ['1', '2', '3', 'flat']  # flat last
        assert rep['k'][3][-1]['is_flat'] and not rep['k'][3][0]['is_flat']

    def test_counts_and_fractions(self, tmp_path):
        images = self._tree(tmp_path)
        motifs = modec.report_data(images, 'job1', {'S1': ['repA']})['S1']['repA']['k'][2]
        # fixture: K=2 -> motif 1 has 10, motif 2 has 20, flat 2 => 32 peptides
        assert [m['n_peptides'] for m in motifs] == [10, 20, 2]
        assert motifs[1]['fraction'] == pytest.approx(20 / 32, abs=1e-3)
        assert sum(m['fraction'] for m in motifs) == pytest.approx(1.0, abs=2e-3)

    def test_image_paths_are_relative_to_static(self, tmp_path):
        images = self._tree(tmp_path)
        motif = modec.report_data(images, 'job1', {'S1': ['repA']})['S1']['repA']['k'][2][0]
        assert motif['image'] == 'images/job1/S1/modec/repA/logos/modec_logos_1of2-001.jpg'

    def test_aic_plot_path_when_present_and_none_when_not(self, tmp_path):
        images = self._tree(tmp_path)
        rep = modec.report_data(images, 'job1', {'S1': ['repA']})['S1']['repA']
        assert rep['aic_plot'] is None
        plot = os.path.join(images, 'job1', 'S1', modec.MODEC_DIRNAME, 'repA', modec.AIC_PLOT_FILENAME)
        with open(plot, 'wb') as f:
            f.write(b'x')
        rep = modec.report_data(images, 'job1', {'S1': ['repA']})['S1']['repA']
        assert rep['aic_plot'] == 'images/job1/S1/modec/repA/aic_plot.png'

    def test_missing_logo_is_none_not_an_error(self, tmp_path):
        images = self._tree(tmp_path, with_logos=False)
        motifs = modec.report_data(images, 'job1', {'S1': ['repA']})['S1']['repA']['k'][1]
        assert all(m['image'] is None for m in motifs)

    def test_replicate_without_output_is_skipped(self, tmp_path):
        images = self._tree(tmp_path)
        data = modec.report_data(images, 'job1', {'S1': ['repA', 'repB'], 'S2': ['x']})
        assert list(data['S1']) == ['repA']
        assert 'S2' not in data


class TestAicPlot:
    def test_writes_a_png_next_to_the_results(self, tmp_path):
        pytest.importorskip('matplotlib')
        out = str(tmp_path / 'out')
        make_output_tree(out, kmax=3)
        path = modec.render_aic_plot(out)
        assert path == os.path.join(out, modec.AIC_PLOT_FILENAME)
        with open(path, 'rb') as f:
            assert f.read(8)[1:4] == b'PNG'  # PNG signature

    def test_missing_summary_raises(self, tmp_path):
        pytest.importorskip('matplotlib')
        os.makedirs(tmp_path / 'out' / 'Summary')
        with pytest.raises(FileNotFoundError):
            modec.render_aic_plot(str(tmp_path / 'out'))

    def test_read_aic_matches_collect_results(self, tmp_path):
        out = str(tmp_path / 'out')
        make_output_tree(out, kmax=3)
        assert modec.read_aic(out) == modec.collect_results(out)['aic']


class TestInformationContent:
    def test_uniform_pwm_has_no_information(self, tmp_path):
        p = tmp_path / 'PWM_K1_flat.txt'
        write_pwm(p, 0)
        assert modec.pwm_information_bits(modec.parse_pwm(p)) < modec.MIN_LOGO_INFORMATION_BITS

    def test_peaked_pwm_has_information(self, tmp_path):
        p = tmp_path / 'PWM_K1_1.txt'
        write_pwm(p, 10, hot={('F', 0): 0.6})
        assert modec.pwm_information_bits(modec.parse_pwm(p)) > 0.5


class TestMhctpAdapter:
    def _pwm(self, tmp_path, **kw):
        p = tmp_path / 'PWM_K1_1.txt'
        write_pwm(p, 10, **kw)
        return modec.parse_pwm(p)

    def test_background_sums_to_one_and_covers_all_residues(self):
        bg = modec.background_frequencies()
        assert set(bg) == set(modec.GIBBS_AA_ORDER)
        assert sum(bg.values()) == pytest.approx(1.0)

    def test_halfbits_are_log_odds_against_the_background(self, tmp_path):
        pwm = self._pwm(tmp_path, hot={('F', 0): 0.6})
        rows = modec.pwm_to_halfbit_rows(pwm)
        assert len(rows) == 9 and all(len(r) == 20 for r in rows)
        f = modec.GIBBS_AA_ORDER.index('F')
        bg_f = modec.background_frequencies()['F']
        assert rows[0][f] == pytest.approx(2 * math.log2(0.6 / bg_f), abs=1e-3)
        assert rows[0][f] > 0  # enriched over background

    def test_scores_are_floored(self, tmp_path):
        pwm = self._pwm(tmp_path, hot={('F', 0): 0.99999})
        assert min(min(r) for r in modec.pwm_to_halfbit_rows(pwm)) >= modec.HALFBIT_FLOOR

    def test_matrix_text_matches_what_mhctp_parses(self, tmp_path):
        """MHC-TP finds the header line by its amino-acid string and reads each data row
        as '<pos> <letter> <20 numbers>'."""
        pwm = self._pwm(tmp_path, hot={('F', 0): 0.6})
        lines = modec.halfbit_matrix_text(pwm).strip().split('\n')
        assert 'halfbits' in lines[0]
        assert 'A R N D C Q E G H I L K M F P S T W Y V' in lines[1]
        data = [ln.split() for ln in lines[2:]]
        assert len(data) == 9
        for pos, parts in enumerate(data, start=1):
            assert len(parts) == 22 and parts[0] == str(pos)
            [float(x) for x in parts[2:]]
        assert data[0][1] == 'F'  # consensus letter of position 1

    def test_kld_is_zero_for_background_and_positive_for_a_motif(self, tmp_path):
        bg = modec.background_frequencies()
        flat = {'positions': 2, 'freq': {aa: [bg[aa]] * 2 for aa in modec.GIBBS_AA_ORDER}}
        assert modec.pwm_kld(flat) == pytest.approx(0.0, abs=1e-9)
        assert modec.pwm_kld(self._pwm(tmp_path, hot={('F', 0): 0.6})) > 0

    def test_kld_table_layout(self, tmp_path):
        out = str(tmp_path / 'out')
        make_output_tree(out, kmax=3)
        lines = modec.kld_table_text(modec.collect_results(out)).strip().split('\n')
        assert lines[0].split('\t') == ['1', '2', '3', '4', '5', '6']
        rows = [ln.split('\t') for ln in lines[1:]]
        assert [r[0] for r in rows] == ['1', '2', '3']
        assert all(len(r) == 7 for r in rows)              # K + six groups
        assert [float(x) for x in rows[0][2:]] == [0.0] * 5  # groups 2-6 are padding for K=1

    def test_writes_matrices_for_real_motifs_only(self, tmp_path):
        out = str(tmp_path / 'out')
        make_output_tree(out, kmax=3)
        n = modec.write_mhctp_inputs(out)
        names = sorted(os.listdir(os.path.join(out, modec.MHCTP_MATRIX_DIRNAME)))
        assert n == 1 + 2 + 3 == len(names)
        assert 'gibbs.2of3.mat' in names and 'gibbs.1of1.mat' in names
        assert not any('flat' in name for name in names)
        assert os.path.isfile(os.path.join(out, modec.MHCTP_IMAGES_DIRNAME, modec.MHCTP_KLD_FILENAME))

    def test_cluster_id_mhctp_derives_is_i_of_k(self, tmp_path):
        """MHC-TP takes the cluster id from the second dot-field of the file name."""
        out = str(tmp_path / 'out')
        make_output_tree(out, kmax=2)
        modec.write_mhctp_inputs(out)
        for name in os.listdir(os.path.join(out, modec.MHCTP_MATRIX_DIRNAME)):
            assert name.split('.')[1] in {'1of1', '1of2', '2of2'}


class TestExportResults:
    def _full(self, tmp_path):
        out = tmp_path / 'out'
        make_output_tree(str(out), kmax=2)
        (out / 'log_K1-2.txt').write_text('log')
        (out / modec.AIC_PLOT_FILENAME).write_bytes(b'png')
        logos = out / modec.LOGO_DIRNAME
        logos.mkdir()
        (logos / 'modec_logos_1of1-001.jpg').write_bytes(b'jpg')
        (logos / 'modec_logos_1of1.eps').write_bytes(b'eps')
        (logos / 'modec_logos_1of1_freq.mat').write_text('mat')
        (out / modec.LOGO_TABLE_DIRNAME).mkdir()
        (out / modec.LOGO_TABLE_DIRNAME / 'x.txt').write_text('t')
        return str(out)

    def test_copies_results_log_graph_and_logo_images(self, tmp_path):
        out = self._full(tmp_path)
        dst = str(tmp_path / 'export' / 'modec' / 'S1' / 'repA')
        modec.export_results(out, dst)
        assert os.path.isfile(os.path.join(dst, 'PWMs', 'PWM_K2_1.txt'))
        assert os.path.isfile(os.path.join(dst, 'Responsibilities', 'bestPepResp_K1.txt'))
        assert os.path.isfile(os.path.join(dst, 'Summary', 'nMotScores_K1-2.txt'))
        assert os.path.isfile(os.path.join(dst, 'log_K1-2.txt'))
        assert os.path.isfile(os.path.join(dst, modec.AIC_PLOT_FILENAME))
        assert os.path.isfile(os.path.join(dst, 'logos', 'modec_logos_1of1-001.jpg'))

    def test_leaves_out_logo_intermediates(self, tmp_path):
        out = self._full(tmp_path)
        dst = str(tmp_path / 'export')
        modec.export_results(out, dst)
        assert not os.path.exists(os.path.join(dst, 'logos', 'modec_logos_1of1.eps'))
        assert not os.path.exists(os.path.join(dst, 'logos', 'modec_logos_1of1_freq.mat'))
        assert not os.path.exists(os.path.join(dst, modec.LOGO_TABLE_DIRNAME))

    def test_is_safe_to_run_twice(self, tmp_path):
        out = self._full(tmp_path)
        dst = str(tmp_path / 'export')
        modec.export_results(out, dst)
        modec.export_results(out, dst)  # no error on an existing folder


class TestLogoTable:
    def test_table_is_transposed_with_symbol_header(self, tmp_path):
        p = tmp_path / 'PWM_K1_1.txt'
        write_pwm(p, 10, hot={('F', 0): 0.6, ('L', 8): 0.5})
        lines = modec.pwm_to_seq2logo_table(modec.parse_pwm(p)).strip().split('\n')
        assert lines[0].split() == list(modec.MODEC_AA_ORDER)
        assert len(lines) == 1 + 9  # header + one row per position
        assert all(len(row.split()) == 20 for row in lines[1:])
        first = [float(x) for x in lines[1].split()]
        assert first[modec.MODEC_AA_ORDER.index('F')] == pytest.approx(0.6, abs=1e-4)
        last = [float(x) for x in lines[9].split()]
        assert last[modec.MODEC_AA_ORDER.index('L')] == pytest.approx(0.5, abs=1e-4)

    def test_command_uses_shannon_type_and_one_based_positions(self):
        cmd = modec.build_logo_command('/t/in.txt', '/o/logo', 'title')
        assert cmd[0] == 'python2'
        assert cmd[cmd.index('-I') + 1] == '1'
        assert cmd[cmd.index('-S') + 1] == '2'
        assert cmd[cmd.index('--format') + 1] == '[JPEG]'
        assert cmd[cmd.index('-o') + 1] == '/o/logo'

    def test_logo_basename_mirrors_gibbs_naming(self):
        assert modec.logo_basename('2', 4) == 'modec_logos_2of4'
        assert modec.logo_basename('flat', 4) == 'modec_logos_flatof4'


@POSIX_ONLY
class TestRenderLogos:
    FAKE = ('import sys\n'
            'a = sys.argv\n'
            'out = a[a.index("-o") + 1]\n'
            'if "FAIL" in a[a.index("-t") + 1]: sys.exit(5)\n'
            'open(out + "-001.jpg", "wb").write(b"jpg")\n')

    def test_one_logo_per_motif_including_flat(self, tmp_path):
        out = str(tmp_path / 'out')
        make_output_tree(out, kmax=3)
        script = tmp_path / 'fake_seq2logo.py'
        script.write_text(self.FAKE)
        n = modec.render_logos(out, seq2logo=str(script), python=sys.executable)
        assert n == (1 + 1) + (2 + 1) + (3 + 1)  # motifs + flat, for K=1..3
        logo_dir = os.path.join(out, modec.LOGO_DIRNAME)
        assert os.path.isfile(os.path.join(logo_dir, 'modec_logos_2of3-001.jpg'))
        assert os.path.isfile(os.path.join(logo_dir, 'modec_logos_flatof3-001.jpg'))

    def test_empty_uniform_flat_motif_is_skipped_not_failed(self, tmp_path):
        """MoDec can output a flat motif with zero peptides and an exactly uniform PWM
        (seen on a mouse job, K=5). Seq2Logo cannot draw it; the job must not fail."""
        out = str(tmp_path / 'out')
        make_output_tree(out, kmax=2)
        write_pwm(os.path.join(out, 'PWMs', 'PWM_K2_flat.txt'), 0)  # uniform, n_pep 0
        script = tmp_path / 'fake_seq2logo.py'
        script.write_text(self.FAKE)
        n = modec.render_logos(out, seq2logo=str(script), python=sys.executable)
        assert n == (1 + 1) + (2 + 0)  # the empty flat motif at K=2 is the only one skipped
        logo_dir = os.path.join(out, modec.LOGO_DIRNAME)
        assert not os.path.exists(os.path.join(logo_dir, 'modec_logos_flatof2-001.jpg'))
        assert os.path.isfile(os.path.join(logo_dir, 'modec_logos_2of2-001.jpg'))

    def test_failure_raises_naming_the_motif(self, tmp_path):
        out = str(tmp_path / 'out')
        make_output_tree(out, kmax=1)
        script = tmp_path / 'fake_seq2logo.py'
        script.write_text(self.FAKE.replace('"FAIL" in', '"motif 1" in'))
        with pytest.raises(RuntimeError, match='modec_logos_1of1'):
            modec.render_logos(out, seq2logo=str(script), python=sys.executable)


@POSIX_ONLY
class TestRunJob:
    @pytest.fixture(autouse=True)
    def _no_logos(self, monkeypatch):
        # The fake binary writes no PWMs; logo rendering has its own tests above.
        monkeypatch.setattr(modec, 'render_logos', lambda out_dir, **kw: 0)
        monkeypatch.setattr(modec, 'render_aic_plot', lambda out_dir: None)
        monkeypatch.setattr(modec, 'export_results', lambda out_dir, export_dir: [])
        monkeypatch.setattr(modec, 'write_mhctp_inputs', lambda out_dir: 0)

    def _cfg(self, binary):
        return types.SimpleNamespace(
            MODEC_BIN=str(binary), MODEC_KMIN=1, MODEC_KMAX=2, MODEC_NRUNS=1,
            MODEC_TIME_LIMIT_SECONDS=30)

    def _job(self, tmp_path, reps=('a', 'b')):
        data = tmp_path / 'data'
        sample = data / 'task1' / 'S1'
        sample.mkdir(parents=True)
        for r in reps:
            (sample / f'{r}{modec.INPUT_SUFFIX}').write_text('PEPTIDEPEPTIDE\n')
            (sample / f'{r}.csv').write_text('Peptide\nPEPTIDEPEPTIDE\n')  # must be ignored
        return str(data), str(tmp_path / 'images')

    def test_unavailable_binary_fails_loudly(self, tmp_path):
        data, images = self._job(tmp_path)
        with pytest.raises(RuntimeError, match='not installed'):
            modec.run_job('task1', data, 9, 'Human', cfg=self._cfg(tmp_path / 'nope'), images_root=images)

    def test_runs_every_replicate_under_the_job_dir(self, tmp_path):
        binary = tmp_path / 'fake_modec'
        make_fake_binary(binary, 'import os, sys\n'
                                 'out = sys.argv[sys.argv.index("-o") + 1]\n'
                                 'os.makedirs(out)\n'
                                 'open(os.path.join(out, "done"), "w").write("ok")\n')
        data, images = self._job(tmp_path)
        modec.run_job('task1', data, 9, 'Human', cfg=self._cfg(binary), images_root=images)
        for rep in ('a', 'b'):
            out = os.path.join(images, 'task1', 'S1', modec.MODEC_DIRNAME, rep)
            assert os.path.isfile(os.path.join(out, 'done'))
            # all output sits under <images>/<taskId>/ so retention cleanup covers it
            assert os.path.realpath(out).startswith(os.path.realpath(os.path.join(images, 'task1')) + os.sep)

    def test_failure_names_the_replicate(self, tmp_path):
        binary = tmp_path / 'fake_modec'
        make_fake_binary(binary, 'import sys\nsys.stderr.write("bad input")\nsys.exit(2)')
        data, images = self._job(tmp_path, reps=('a',))
        with pytest.raises(RuntimeError, match=r'S1/a.*bad input'):
            modec.run_job('task1', data, 9, 'Human', cfg=self._cfg(binary), images_root=images)

    def test_no_input_files_is_an_error(self, tmp_path):
        binary = tmp_path / 'fake_modec'
        make_fake_binary(binary, 'pass')
        data = tmp_path / 'data'
        (data / 'task1' / 'S1').mkdir(parents=True)
        with pytest.raises(RuntimeError, match='No'):
            modec.run_job('task1', str(data), 9, 'Human', cfg=self._cfg(binary), images_root=str(tmp_path / 'i'))


class TestModecAvailable:
    def test_missing_path(self, tmp_path):
        assert not modec.modec_available(str(tmp_path / 'nope'))
        assert not modec.modec_available('')
        assert not modec.modec_available(None)

    @POSIX_ONLY
    def test_non_executable_file(self, tmp_path):
        f = tmp_path / 'MoDec_unix'
        f.write_text('x')
        f.chmod(0o644)
        assert not modec.modec_available(str(f))
        f.chmod(0o755)
        assert modec.modec_available(str(f))
