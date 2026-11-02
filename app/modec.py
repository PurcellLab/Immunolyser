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
import shutil
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

    returncode, stderr = _run_group(cmd, time_limit, what=f'MoDec ({out_dir})')
    if returncode != 0:
        raise RuntimeError(f'MoDec failed (rc={returncode}): {(stderr or "")[-500:] or "no error output"}')


def _run_group(cmd, time_limit, what):
    """Run cmd in its own session with a wall-clock cap; return (returncode, stderr).

    On timeout, or if the caller is interrupted, the whole process group is killed
    so no child outlives the call. Shared by MoDec and the logo renderer.
    """
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
        raise TimeoutError(f'{what} exceeded {time_limit}s and was killed')
    except BaseException:  # includes Celery's SoftTimeLimitExceeded
        _kill_group(proc)
        raise
    return proc.returncode, stderr


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


def read_aic(out_dir):
    """The AIC table of one replicate's MoDec output folder (see parse_aic)."""
    summary_dir = os.path.join(out_dir, 'Summary')
    aic_files = sorted(f for f in os.listdir(summary_dir) if f.startswith('nMotScores'))
    if not aic_files:
        raise FileNotFoundError(f'no nMotScores file in {summary_dir}')
    return parse_aic(os.path.join(summary_dir, aic_files[0]))


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

    aic = read_aic(out_dir)

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
# Logos (rendered here with the Seq2Logo already in the image; we do not use
# MoDec's own R report, which would need ggseqlogoMOD/htmlTable installed)
# --------------------------------------------------------------------------

# Logo files sit beside the PWMs, named like GibbsCluster's so the report code can
# treat them alike: modec_logos_<motif>of<K>-001.jpg, motif = 1..K or 'flat'.
LOGO_DIRNAME = 'logos'
LOGO_TABLE_DIRNAME = 'logo_tables'
SEQ2LOGO_SCRIPT = os.path.join('app', 'tools', 'seq2logo-2.1', 'Seq2Logo.py')
# Seq2Logo is python2 (see CLAUDE.md: python2.7 is a hard runtime dependency).
SEQ2LOGO_PYTHON = 'python2'


def logo_basename(motif_id, k):
    return f'modec_logos_{motif_id}of{k}'


def pwm_to_seq2logo_table(pwm):
    """A MoDec PWM as a Seq2Logo frequency table: a header of amino-acid symbols,
    then one row per position (the transpose of MoDec's own layout)."""
    lines = [' '.join(MODEC_AA_ORDER)]
    for pos in range(pwm['positions']):
        lines.append(' '.join(f"{pwm['freq'][aa][pos]:.6g}" for aa in MODEC_AA_ORDER))
    return '\n'.join(lines) + '\n'


def build_logo_command(table_path, out_prefix, title, seq2logo=SEQ2LOGO_SCRIPT, python=SEQ2LOGO_PYTHON):
    """Seq2Logo argv for one motif.

    -I 1 (Shannon): letter height = frequency x information content, no
       background model, so letters never dip below zero. This is the closest match
       to MoDec's own logos (checked against its reference output); the
       Kullback-Leibler type used for Gibbs logos draws weak letters below the axis.
    -S 2: number positions from 1 (Seq2Logo's default labels them from 0).
    """
    return [python, seq2logo, '-f', table_path, '-o', out_prefix,
            '-I', '1', '-S', '2', '--format', '[JPEG]', '-t', title]


AIC_PLOT_FILENAME = 'aic_plot.png'


def render_aic_plot(out_dir):
    """Draw AIC against the number of motifs for one replicate (the MoDec counterpart
    of GibbsCluster's KLD bar plot: it is what you look at to choose how many motifs
    to show). Lower AIC is better; the lowest is marked. Drawn here, as an image,
    so the offline export can inline it like every other figure.

    Returns the path written. matplotlib is imported lazily with the headless
    backend so importing this module stays cheap and display-free.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    rows = read_aic(out_dir)
    best = suggested_k(rows)
    ks = [r['n'] for r in rows]
    aics = [r['AIC'] for r in rows]

    fig, ax = plt.subplots(figsize=(5, 4), dpi=100)
    ax.plot(ks, aics, marker='o', color='#4a6fa5', linewidth=2, zorder=2)
    for r in rows:
        is_best = r['n'] == best
        if is_best:
            ax.scatter([r['n']], [r['AIC']], s=170, color='#2e8b57', zorder=3)
        ax.annotate(f"{r['AIC']:.0f}", (r['n'], r['AIC']), textcoords='offset points',
                    xytext=(0, 10), ha='center', fontsize=8,
                    fontweight='bold' if is_best else 'normal')
    ax.set_xticks(ks)
    ax.set_xlabel('Number of motifs (K)')
    ax.set_ylabel('AIC (lower is better)')
    ax.set_title(f'MoDec model selection: lowest AIC at K={best}', fontsize=10)
    ax.margins(y=0.15)
    ax.grid(axis='y', alpha=0.3)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    fig.tight_layout()

    path = os.path.join(out_dir, AIC_PLOT_FILENAME)
    fig.savefig(path)
    plt.close(fig)
    return path


def pwm_information_bits(pwm):
    """Total information content (bits, summed over positions) of a PWM against a
    uniform background: 0 for a motif that is exactly flat."""
    total = 0.0
    for pos in range(pwm['positions']):
        entropy = 0.0
        for aa in MODEC_AA_ORDER:
            p = pwm['freq'][aa][pos]
            if p > 0:
                entropy -= p * math.log2(p)
        total += math.log2(len(MODEC_AA_ORDER)) - entropy
    return total


# A motif below this has nothing to draw. MoDec's flat motif can be exactly uniform with
# zero peptides assigned (seen on a mouse job at K=5); Seq2Logo then aborts with "No
# PSSM data was found", which must not fail the whole job.
MIN_LOGO_INFORMATION_BITS = 1e-6


def render_logos(out_dir, time_limit=120, seq2logo=SEQ2LOGO_SCRIPT, python=SEQ2LOGO_PYTHON):
    """Draw a logo for every motif of every K in one replicate's MoDec output.

    Returns the number of logos written. Raises if Seq2Logo fails on a motif that has
    something to draw, so a broken renderer fails the job instead of leaving a report
    with silently missing images. A motif with no information at all (an empty, exactly
    uniform flat motif) is skipped: the report shows "No peptides assigned" for it.
    """
    results = collect_results(out_dir)
    table_dir = os.path.join(out_dir, LOGO_TABLE_DIRNAME)
    logo_dir = os.path.join(out_dir, LOGO_DIRNAME)
    os.makedirs(table_dir, exist_ok=True)
    os.makedirs(logo_dir, exist_ok=True)

    written = 0
    for k, entry in results['k'].items():
        for motif_id, pwm in entry['motifs'].items():
            if pwm_information_bits(pwm) < MIN_LOGO_INFORMATION_BITS:
                continue
            base = logo_basename(motif_id, k)
            table_path = os.path.join(table_dir, f'{base}.txt')
            with open(table_path, 'w') as f:
                f.write(pwm_to_seq2logo_table(pwm))
            title = f'K={k} motif {motif_id} ({pwm["n_peptides"]} peptides)'
            cmd = build_logo_command(table_path, os.path.join(logo_dir, base), title, seq2logo, python)
            rc, stderr = _run_group(cmd, time_limit, what=f'Seq2Logo ({base})')
            jpg = os.path.join(logo_dir, f'{base}-001.jpg')
            if rc != 0 or not os.path.isfile(jpg):
                raise RuntimeError(f'Seq2Logo failed for {base} (rc={rc}): {(stderr or "")[-300:]}')
            written += 1
    return written


# --------------------------------------------------------------------------
# MHC-TP input (allotype prediction for MoDec motifs)
# --------------------------------------------------------------------------
#
# MHC-TP (clust-search) reads a GibbsCluster-shaped folder: matrices/gibbs.<i>of<K>.mat,
# images/gibbs.KLDvsClusters.tab and logos/*<i>of<K>*.jpg. The logos we already have
# (their file names contain the <i>of<K> id it matches on); the other two are written
# here so MHC-TP runs on MoDec output unchanged. The cluster id it reports is read
# from the second dot-field of the matrix file name, hence the "gibbs." prefix.
#
# A GibbsCluster .mat is a log-odds PSSM in HALFBITS (2*log2(p/background)), and MHC-TP
# correlates it against reference motifs of that same kind. MoDec's PWMs are
# frequencies, so they are converted; correlating raw frequencies against halfbit
# scores would compare a linear with a log scale.

MHCTP_MATRIX_DIRNAME = 'matrices'
MHCTP_IMAGES_DIRNAME = 'images'
MHCTP_KLD_FILENAME = 'gibbs.KLDvsClusters.tab'

# UniProtKB/Swiss-Prot amino-acid composition (%), used as the background. GibbsCluster
# bakes its own pre-calculated UniProt background into its binary, so it cannot be
# reused as a file; these published values are the closest stand-in.
BACKGROUND_PERCENT = {
    'A': 8.25, 'R': 5.53, 'N': 4.06, 'D': 5.45, 'C': 1.37, 'Q': 3.93, 'E': 6.75,
    'G': 7.07, 'H': 2.27, 'I': 5.96, 'L': 9.66, 'K': 5.84, 'M': 2.42, 'F': 3.86,
    'P': 4.70, 'S': 6.56, 'T': 5.34, 'W': 1.08, 'Y': 2.92, 'V': 6.87,
}
# GibbsCluster PSSMs bottom out near -7.5 halfbits (its pseudocounts keep rare residues
# off zero). An unfloored log of a near-zero MoDec frequency would reach about -17 and
# dominate the per-position Pearson correlation MHC-TP averages.
HALFBIT_FLOOR = -8.0


def background_frequencies():
    """Background amino-acid frequencies, normalised to sum to 1."""
    total = sum(BACKGROUND_PERCENT.values())
    return {aa: v / total for aa, v in BACKGROUND_PERCENT.items()}


def pwm_to_halfbit_rows(pwm):
    """Per position, 20 log-odds scores in halfbits in GIBBS_AA_ORDER, floored."""
    bg = background_frequencies()
    rows = []
    for pos in range(pwm['positions']):
        rows.append([
            max(HALFBIT_FLOOR, 2 * math.log2(max(pwm['freq'][aa][pos], 1e-12) / bg[aa]))
            for aa in GIBBS_AA_ORDER
        ])
    return rows


def pwm_kld(pwm):
    """Kullback-Leibler divergence (bits) of a motif from the background, summed over
    positions: the quantity GibbsCluster reports per group and MHC-TP displays."""
    bg = background_frequencies()
    total = 0.0
    for pos in range(pwm['positions']):
        for aa in GIBBS_AA_ORDER:
            p = pwm['freq'][aa][pos]
            if p > 0:
                total += p * math.log2(p / bg[aa])
    return total


def halfbit_matrix_text(pwm):
    """A GibbsCluster-style .mat: comment line, amino-acid header, then one row per
    position as '<pos> <consensus> <20 scores>'."""
    lines = ['Last position-specific scoring matrix computed, values are in halfbits',
             ' '.join(GIBBS_AA_ORDER)]
    for pos, scores in enumerate(pwm_to_halfbit_rows(pwm), start=1):
        consensus = GIBBS_AA_ORDER[scores.index(max(scores))]
        lines.append(f'{pos} {consensus} ' + ' '.join(f'{v:.3f}' for v in scores))
    return '\n'.join(lines) + '\n'


def kld_table_text(results, groups=6):
    """gibbs.KLDvsClusters.tab: a header of group numbers, then one row per motif count
    K as 'K<TAB>kld_of_group_1<TAB>...' padded with zeros to `groups` columns."""
    lines = ['\t'.join(str(g) for g in range(1, groups + 1))]
    for k, entry in results['k'].items():
        klds = [pwm_kld(entry['motifs'][str(i)]) for i in range(1, k + 1) if str(i) in entry['motifs']]
        klds += [0] * (groups - len(klds))
        lines.append('\t'.join([str(k)] + [f'{v:.6g}' for v in klds[:groups]]))
    return '\n'.join(lines) + '\n'


def write_mhctp_inputs(out_dir):
    """Write the matrices and KLD table MHC-TP needs into one replicate's MoDec folder.

    The flat (contaminant) motif is not written: it is not a binding motif and would
    only produce meaningless allotype matches. Returns the number of matrices written.
    """
    results = collect_results(out_dir)
    matrix_dir = os.path.join(out_dir, MHCTP_MATRIX_DIRNAME)
    images_dir = os.path.join(out_dir, MHCTP_IMAGES_DIRNAME)
    os.makedirs(matrix_dir, exist_ok=True)
    os.makedirs(images_dir, exist_ok=True)

    written = 0
    for k, entry in results['k'].items():
        for motif_id, pwm in entry['motifs'].items():
            if motif_id == 'flat':
                continue
            with open(os.path.join(matrix_dir, f'gibbs.{motif_id}of{k}.mat'), 'w') as f:
                f.write(halfbit_matrix_text(pwm))
            written += 1
    with open(os.path.join(images_dir, MHCTP_KLD_FILENAME), 'w') as f:
        f.write(kld_table_text(results))
    return written


# --------------------------------------------------------------------------
# Raw-results export
# --------------------------------------------------------------------------

# Where a job's downloadable files live. zip_job_exports() zips this folder for the
# "Download Data (ZIP)" link, and it sits under the job id, so the retention sweep
# deletes it with the job.
EXPORT_DIRNAME = 'export'
EXPORT_SUBDIR = 'modec'


def export_results(out_dir, export_dir):
    """Copy what a user would want from one replicate's MoDec run into export_dir.

    The PWMs, per-peptide responsibilities and AIC summary are MoDec's own results;
    the log, the AIC graph and the logo images make the folder self-explanatory.
    The logos' intermediate .eps/.txt/.mat files and the logo input tables are left
    out. Returns the list of top-level entries written.
    """
    os.makedirs(export_dir, exist_ok=True)
    written = []
    for name in ('PWMs', 'Responsibilities', 'Summary'):
        src = os.path.join(out_dir, name)
        if os.path.isdir(src):
            shutil.copytree(src, os.path.join(export_dir, name), dirs_exist_ok=True)
            written.append(name)
    for fname in sorted(os.listdir(out_dir)):
        if fname.startswith('log_') or fname == AIC_PLOT_FILENAME:
            shutil.copy2(os.path.join(out_dir, fname), os.path.join(export_dir, fname))
            written.append(fname)
    logo_src = os.path.join(out_dir, LOGO_DIRNAME)
    if os.path.isdir(logo_src):
        logo_dst = os.path.join(export_dir, LOGO_DIRNAME)
        os.makedirs(logo_dst, exist_ok=True)
        for fname in sorted(os.listdir(logo_src)):
            if fname.endswith('.jpg'):
                shutil.copy2(os.path.join(logo_src, fname), os.path.join(logo_dst, fname))
        written.append(LOGO_DIRNAME)
    return written


# --------------------------------------------------------------------------
# Report data (shared by the report page and the offline export)
# --------------------------------------------------------------------------

def read_job_method(images_root, task_id):
    """The clustering method recorded for a job; GibbsCluster if none was recorded
    (jobs from before MoDec existed have no marker file)."""
    from constants import ClusteringMethod
    try:
        with open(os.path.join(images_root, task_id, 'clustering_method.txt')) as f:
            return ClusteringMethod.normalise(f.read())
    except OSError:
        return ClusteringMethod.DEFAULT


def motif_peptides(rep_dir, k, motif_id):
    """The bestPepResp rows whose best motif (at this K) is `motif_id`.

    motif_id is a motif number or 'flat'; MoDec numbers its flat (contaminant)
    motif 0 in the responsibilities files. Raises FileNotFoundError if that K
    was not run.
    """
    target = 0 if motif_id == 'flat' else int(motif_id)
    path = os.path.join(rep_dir, 'Responsibilities', f'bestPepResp_K{int(k)}.txt')
    return [row for row in parse_best_pep_resp(path) if row['Motif_1st'] == target]


def _motif_sort_key(motif_id):
    # numbered motifs first in order, the flat (contaminant) motif last
    return (motif_id == 'flat', int(motif_id) if motif_id != 'flat' else 0)


def report_data(images_root, task_id, replicates_by_sample):
    """What the MoDec report section renders, per sample and replicate.

    replicates_by_sample: {sample: [replicate stem, ...]}.
    Returns {sample: {replicate: {'aic': [...], 'suggested_k': int|None,
                                  'k': {K: [motif, ...]}}}} where each motif is
    {'id', 'is_flat', 'n_peptides', 'fraction', 'image'} and 'image' is relative
    to app/static (None if the logo is missing). Replicates MoDec produced nothing
    for are skipped, like the Gibbs report skips replicates without a bar plot.
    """
    out = {}
    for sample in sorted(replicates_by_sample):
        for stem in sorted(replicates_by_sample[sample]):
            rep_dir = os.path.join(images_root, task_id, sample, MODEC_DIRNAME, stem)
            if not os.path.isdir(os.path.join(rep_dir, 'PWMs')):
                continue
            res = collect_results(rep_dir)
            by_k = {}
            for k, entry in res['k'].items():
                total = len(entry['assignments'])
                motifs = []
                for motif_id in sorted(entry['motifs'], key=_motif_sort_key):
                    n = entry['motifs'][motif_id]['n_peptides']
                    logo = os.path.join(rep_dir, LOGO_DIRNAME, f'{logo_basename(motif_id, k)}-001.jpg')
                    rel = os.path.relpath(logo, os.path.join(images_root, os.pardir)) if os.path.isfile(logo) else None
                    motifs.append({
                        'id': motif_id,
                        'is_flat': motif_id == 'flat',
                        'n_peptides': n,
                        'fraction': round(n / total, 3) if total else 0.0,
                        'image': rel.replace(os.sep, '/') if rel else None,
                    })
                by_k[k] = motifs
            plot = os.path.join(rep_dir, AIC_PLOT_FILENAME)
            plot_rel = os.path.relpath(plot, os.path.join(images_root, os.pardir)) if os.path.isfile(plot) else None
            out.setdefault(sample, {})[stem] = {
                'aic': res['aic'], 'suggested_k': res['suggested_k'], 'k': by_k,
                'aic_plot': plot_rel.replace(os.sep, '/') if plot_rel else None}
    return out


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
        render_logos(out_dir)
        render_aic_plot(out_dir)
        write_mhctp_inputs(out_dir)
        sample, stem = label.split('/', 1)
        export_results(out_dir, os.path.join(images_root, task_id, EXPORT_DIRNAME, EXPORT_SUBDIR, sample, stem))
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
