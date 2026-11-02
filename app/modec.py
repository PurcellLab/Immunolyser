"""MoDec (MHC class II motif deconvolution) — run it and read its output.

MoDec (Gfeller Lab; Racle et al., Nat Biotechnol 2019) is a binary-only,
academic-non-commercial tool, so nothing here ships or vendors it: the Dockerfile
fetches it at build time and `modec_available()` gates the whole feature.

It is the class II alternative to GibbsCluster. Like gibbscluster.py this file is
both a library (the parsers/command builder are pure and unit-tested without the
binary) and a python3 shim that `submit_job` calls as a subprocess:

    python3 app/modec.py <taskId> <data_mount> <motif_length> <species>
"""
import csv
import math
import os
import signal
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

# Output sub-folder under app/static/images/<taskId>/<sample>/, sibling of
# 'gibbscluster'. Everything MoDec writes lives under the job id, so the existing
# retention sweep (cleanup_expired_jobs) removes it with the job.
MODEC_DIRNAME = 'modec'

# n-mer file the class II pipeline already writes per replicate (12-20-mers,
# duplicates removed); same input GibbsCluster gets.
INPUT_SUFFIX = '_12to20mer.txt'

# MoDec's own amino-acid order in PWM files (alphabetical one-letter codes).
MODEC_AA_ORDER = 'ACDEFGHIKLMNPQRSTVWY'
# Order MHC-TP / GibbsCluster matrices use.
GIBBS_AA_ORDER = 'ARNDCQEGHILKMFPSTWYV'

# --theta_norm background-frequency choice. 1 = HLA-II peptidomics (MoDec's
# default, kept for human), 4 = mouse MHC-II peptidomics.
THETA_NORM_BY_SPECIES = {'human': 1, 'mouse': 4}


# --------------------------------------------------------------------------
# Environment helpers
# --------------------------------------------------------------------------

def modec_available(binary):
    """True when the MoDec binary is installed and executable."""
    return bool(binary) and os.path.isfile(binary) and os.access(binary, os.X_OK)


def _cgroup_cpu_quota(root='/sys/fs/cgroup'):
    """CPUs granted by the container's cgroup quota, or None if unlimited/unknown."""
    try:  # cgroup v2: "<quota> <period>" or "max <period>"
        with open(os.path.join(root, 'cpu.max')) as f:
            quota, period = f.read().split()
        if quota != 'max':
            return int(quota) / int(period)
        return None
    except (OSError, ValueError):
        pass
    try:  # cgroup v1
        with open(os.path.join(root, 'cpu', 'cpu.cfs_quota_us')) as f:
            quota = int(f.read())
        with open(os.path.join(root, 'cpu', 'cpu.cfs_period_us')) as f:
            period = int(f.read())
        if quota > 0:
            return quota / period
    except (OSError, ValueError):
        pass
    return None


def available_cpus(cgroup_root='/sys/fs/cgroup'):
    """CPUs this process can actually use.

    os.cpu_count() reports the *host's* cores inside a container. In production it
    said 16 against an 8-CPU quota, so GibbsCluster was started with -k 16 and
    oversubscribed the worker by 2x. Honour affinity and the cgroup quota instead.
    """
    try:
        n = len(os.sched_getaffinity(0))
    except AttributeError:  # not available on Windows/macOS
        n = os.cpu_count() or 1
    quota = _cgroup_cpu_quota(cgroup_root)
    if quota:
        n = min(n, max(1, math.ceil(quota)))
    return max(1, n)


# --------------------------------------------------------------------------
# Command construction and execution
# --------------------------------------------------------------------------

def build_command(binary, peptide_file, out_dir, *, kmin, kmax, nruns,
                  motif_length, species=None):
    """The MoDec argv for one replicate. Pure, so it is testable without MoDec.

    --MHC2 is always passed: it sets --pepLmin 12 (drops class I contaminants) and
    --specInits (hydrophobic-P1-biased starts), both recommended for HLA-II. It
    also multiplies the real run count by 5, which is why nruns is kept small.
    --log sends MoDec's progress to a file in the result folder instead of a pipe
    we would have to drain.
    """
    if not 1 <= int(kmin) <= int(kmax):
        raise ValueError(f'invalid K range {kmin}..{kmax}')
    if int(nruns) < 1:
        raise ValueError('nruns must be >= 1')
    cmd = [
        binary,
        '-i', peptide_file,
        '-o', out_dir,
        '--MHC2',
        '--Kmin', str(int(kmin)),
        '--Kmax', str(int(kmax)),
        '--nruns', str(int(nruns)),
        '-L', str(int(motif_length)),
        '--log',
    ]
    theta_norm = THETA_NORM_BY_SPECIES.get((species or '').strip().lower())
    if theta_norm is not None and theta_norm != 1:
        cmd += ['--theta_norm', str(theta_norm)]
    return cmd


def run_modec(cmd, out_dir, time_limit):
    """Run one MoDec command; raise on failure; never leave children behind.

    MoDec refuses to run if its output folder exists, so a leftover folder is an
    error rather than something to overwrite. The process runs in its own session
    so that on timeout (or if the calling Celery task is cancelled) the whole
    process group can be killed — Celery's time limit kills only the task process,
    and orphaned tool children kept eating CPU for weeks once before.
    """
    if os.path.exists(out_dir):
        raise FileExistsError(f'MoDec output folder already exists: {out_dir}')
    os.makedirs(os.path.dirname(out_dir), exist_ok=True)

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    try:
        _, stderr = proc.communicate(timeout=time_limit)
    except subprocess.TimeoutExpired:
        _kill_group(proc)
        raise TimeoutError(f'MoDec exceeded {time_limit}s and was killed: {out_dir}')
    except BaseException:  # includes Celery's SoftTimeLimitExceeded
        _kill_group(proc)
        raise
    if proc.returncode != 0:
        raise RuntimeError(f'MoDec failed (rc={proc.returncode}): {(stderr or "")[-500:] or "no error output"}')


def _kill_group(proc):
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        pass
    proc.communicate()


# --------------------------------------------------------------------------
# Output parsing
# --------------------------------------------------------------------------

def parse_pwm(path):
    """Read one PWM_K<k>_<i|flat>.txt.

    Layout: a '#<name>\\t<n_peptides>' line, a tab-separated header of position
    numbers, then one row per amino acid with one frequency per position. Returns
    {'n_peptides': int, 'positions': int, 'freq': {aa: [per-position values]}}.
    """
    with open(path) as f:
        lines = [ln.rstrip('\n') for ln in f if ln.strip()]
    if len(lines) < 3 or not lines[0].startswith('#'):
        raise ValueError(f'not a MoDec PWM file: {path}')
    n_peptides = int(lines[0].split('\t')[-1])
    freq = {}
    for ln in lines[2:]:
        cells = ln.split('\t')
        freq[cells[0]] = [float(x) for x in cells[1:]]
    positions = len(next(iter(freq.values())))
    if set(freq) != set(MODEC_AA_ORDER):
        raise ValueError(f'unexpected amino-acid rows in {path}: {sorted(freq)}')
    return {'n_peptides': n_peptides, 'positions': positions, 'freq': freq}


def pwm_as_gibbs_rows(pwm):
    """PWM as a positions x 20 table in GibbsCluster/MHC-TP amino-acid order."""
    return [[pwm['freq'][aa][pos] for aa in GIBBS_AA_ORDER] for pos in range(pwm['positions'])]


def parse_best_pep_resp(path):
    """Read Responsibilities/bestPepResp_K<k>.txt into a list of typed dicts.

    One row per peptide: its best (and second best) motif, offset and 9-mer core.
    Motif 0 is MoDec's flat (contaminant) motif.
    """
    int_cols = {'Length', 'Motif_1st', 'P1_1st', 'Motif_2nd', 'P1_2nd'}
    float_cols = {'Score', 'Pep_weight', 'resp_1st', 'resp_2nd'}
    rows = []
    with open(path, newline='') as f:
        for rec in csv.DictReader(f, delimiter='\t'):
            row = {}
            for key, val in rec.items():
                if key in int_cols:
                    row[key] = int(val)
                elif key in float_cols:
                    row[key] = float(val)
                else:
                    row[key] = val
            rows.append(row)
    return rows


def parse_aic(path):
    """Read Summary/nMotScores_K<min>-<max>.txt -> [{'n', 'logL', 'AIC'}], by K.

    Lower AIC is better. MoDec's authors advise AIC as a guide only and manual
    inspection of the motifs for the final K.
    """
    out = []
    with open(path, newline='') as f:
        for rec in csv.DictReader(f, delimiter='\t'):
            out.append({'n': int(rec['n']), 'logL': float(rec['logL']), 'AIC': float(rec['AIC'])})
    return sorted(out, key=lambda r: r['n'])


def suggested_k(aic_rows):
    """K with the lowest AIC (None if there are no scores)."""
    return min(aic_rows, key=lambda r: r['AIC'])['n'] if aic_rows else None


def collect_results(out_dir):
    """Everything the report needs from one replicate's MoDec output folder.

    {'aic': [...], 'suggested_k': int|None,
     'k': {K: {'motifs': {id: {'n_peptides', 'positions', 'freq'}}, 'assignments': [...]}}}
    where a motif id is its number as a string, or 'flat'.
    """
    pwm_dir = os.path.join(out_dir, 'PWMs')
    resp_dir = os.path.join(out_dir, 'Responsibilities')
    summary_dir = os.path.join(out_dir, 'Summary')

    aic_files = sorted(f for f in os.listdir(summary_dir) if f.startswith('nMotScores'))
    if not aic_files:
        raise FileNotFoundError(f'no nMotScores file in {summary_dir}')
    aic = parse_aic(os.path.join(summary_dir, aic_files[0]))

    by_k = {}
    for fname in sorted(os.listdir(pwm_dir)):
        # PWM_K<k>_<id>.txt
        stem = fname[:-4] if fname.endswith('.txt') else None
        parts = stem.split('_') if stem else []
        if len(parts) != 3 or parts[0] != 'PWM' or not parts[1].startswith('K'):
            continue
        k = int(parts[1][1:])
        by_k.setdefault(k, {'motifs': {}, 'assignments': []})
        by_k[k]['motifs'][parts[2]] = parse_pwm(os.path.join(pwm_dir, fname))
    for k, entry in by_k.items():
        entry['assignments'] = parse_best_pep_resp(os.path.join(resp_dir, f'bestPepResp_K{k}.txt'))

    return {'aic': aic, 'suggested_k': suggested_k(aic), 'k': dict(sorted(by_k.items()))}


# --------------------------------------------------------------------------
# Whole-job runner (the python3 shim body)
# --------------------------------------------------------------------------

def _settings():
    """Config values, read lazily so importing this module never needs the env."""
    root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    if root not in sys.path:
        sys.path.insert(0, root)
    from config import Config
    return Config


def replicate_jobs(task_id, data_mount, images_root):
    """(input_file, out_dir, label) for every class II n-mer file of the job."""
    task_path = os.path.join(data_mount, task_id)
    jobs = []
    for sample in sorted(os.listdir(task_path)):
        sample_path = os.path.join(task_path, sample)
        if not os.path.isdir(sample_path):
            continue
        for fname in sorted(os.listdir(sample_path)):
            if fname.endswith(INPUT_SUFFIX):
                stem = fname[:-len(INPUT_SUFFIX)]
                jobs.append((
                    os.path.join(sample_path, fname),
                    os.path.join(images_root, task_id, sample, MODEC_DIRNAME, stem),
                    f'{sample}/{stem}',
                ))
    return jobs


def run_job(task_id, data_mount, motif_length, species, cfg=None, images_root=None):
    """Run MoDec on every replicate of a job; raise if any replicate fails.

    Replicates run concurrently, bounded by the CPUs this container really has.
    Failures are collected and raised together so the job is marked FAILURE with a
    reason (GibbsCluster failures, by contrast, were silent empty output).
    """
    cfg = cfg or _settings()
    if not modec_available(cfg.MODEC_BIN):
        raise RuntimeError('MoDec was requested but its binary is not installed on this server.')
    images_root = images_root or os.path.join('app', 'static', 'images')

    jobs = replicate_jobs(task_id, data_mount, images_root)
    if not jobs:
        raise RuntimeError(f'No {INPUT_SUFFIX} files found for job {task_id}.')

    def one(job):
        peptide_file, out_dir, label = job
        cmd = build_command(
            cfg.MODEC_BIN, peptide_file, out_dir,
            kmin=cfg.MODEC_KMIN, kmax=cfg.MODEC_KMAX, nruns=cfg.MODEC_NRUNS,
            motif_length=motif_length, species=species,
        )
        print('MoDec command:', ' '.join(cmd), flush=True)
        run_modec(cmd, out_dir, cfg.MODEC_TIME_LIMIT_SECONDS)
        return label

    errors = []
    with ThreadPoolExecutor(max_workers=min(len(jobs), available_cpus())) as pool:
        futures = [(job[2], pool.submit(one, job)) for job in jobs]
        for label, fut in futures:
            try:
                fut.result()
            except Exception as exc:  # collect, report every replicate
                errors.append(f'{label}: {exc}')
    if errors:
        raise RuntimeError('MoDec failed for ' + '; '.join(errors))


if __name__ == '__main__':
    run_job(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4] if len(sys.argv) > 4 else None)
