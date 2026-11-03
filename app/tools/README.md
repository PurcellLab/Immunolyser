# Bioinformatics Tools

This directory holds the bioinformatics prediction tools used by Immunolyser.

The Dockerfile / `docker-compose.yml` in the repo root are a single setup path
that works for both a local/self-hosted install and our own server deployment
— only the `.env` values (`DATA_VOLUME`, `HLA_PEPCLUST_REF_DATA`) differ. A
local user runs `docker compose up --build` with the defaults; a server
deployment points those two variables at persistent host paths instead.

## Tools included in the repository (open-source)

The following tool archives are committed to this repo and extracted automatically during `docker build`:

| Tool | Archive |
|------|---------|
| seq2logo 2.1 | `seq2logo-2.1.all.tar.gz` |
| GibbsCluster 2.0f | `gibbscluster-2.0f.Linux.tar.gz` |

MixMHCpred, MixMHC2pred, and HLA-PepClust (MHC-TP) are downloaded from their public GitHub releases at build time.

## MoDec (optional, academic non-commercial licence)

MoDec 1.2 (Gfeller Lab; Racle et al., *Nat Biotechnol* 2019;37:1283-1286) is the
MHC class II motif-deconvolution alternative to GibbsCluster, selectable on the
submit form for class II jobs. It is **binary-only** and licensed by the Ludwig
Institute for Cancer Research for **academic, non-commercial research only**
(<https://github.com/GfellerLab/MoDec/blob/v1.2.1/LICENSE>). For-profit use needs a
separate licence from LICR (contact in the MoDec README).

- It is **downloaded at `docker build` time** (the `MODEC_TAG` / `MODEC_SHA256` build args in
  the Dockerfile) into `app/tools/MoDec/MoDec_unix`, never committed. We may not
  redistribute it; each deployer fetches their own copy.
- A failed download only warns. The app then hides the MoDec option (it checks
  `MODEC_BIN` exists and is executable), so a local install without it still works.
  A **checksum mismatch fails the build**.
- The licence requires every user of a server running MoDec to have agreed to its
  terms. The submit form therefore requires a tick box when MoDec is chosen, enforced
  server-side and recorded per job. Keep that if you host this publicly.
- No R packages are needed: we render the motif logos ourselves from MoDec's PWM files
  and do not use MoDec's own `--makeReport` HTML report.
- **MHC-TP runs on MoDec output unchanged.** MHC-TP reads a GibbsCluster-shaped folder, so
  `app/modec.py` writes `matrices/gibbs.<i>of<K>.mat` and `images/gibbs.KLDvsClusters.tab`
  into each replicate's MoDec folder (the logos it also needs are already there). A
  GibbsCluster matrix is a log-odds PSSM in halfbits, so MoDec's frequency PWMs are
  converted (`2*log2(p/background)` against UniProt/Swiss-Prot frequencies, floored at
  -8 as GibbsCluster's own matrices are). On the same 1,500 HLA-II peptides MoDec and
  GibbsCluster give the same allotype families (DRB1*04, DRB1*01, DRB4) at comparable
  correlations. The flat motif is not submitted to MHC-TP. The reference data for class II
  (`Gibbs_motifs_*_classII`, `*_classii.db`, a warmed `numba_cache`) must be in
  `HLA_PEPCLUST_REF_DATA`, as for GibbsCluster jobs.
- To use another build, set `MODEC_BIN` in `.env`.

## Licensed tools — you must download these before building

**netMHCpan 4.2** and **netMHCIIpan 4.3** require a free academic licence from DTU Health Tech.

1. Go to https://services.healthtech.dtu.dk/software.php
2. Register and accept the licence for each tool
3. Download the Linux tarballs:
   - `netMHCpan-4.2c.Linux.tar.gz`
   - `netMHCIIpan-4.3j.Linux.tar.gz`
4. Place them in this directory (`app/tools/`) before running `docker build`

The Dockerfile detects whether the tarballs are present and prints a warning (not an error) if they are missing. Without them, NetMHCpan-based predictions will not be available.

## Ref data for HLA-PepClust

The large reference databases (~1 GB each) are NOT in the repository and are NOT
baked into the image. `docker-compose.yml` bind-mounts a host directory to
`app/tools/HLA-PepClust/data/ref_data` inside the container — by default
`./hla-pepclust-ref-data` next to `docker-compose.yml`, or set
`HLA_PEPCLUST_REF_DATA=/path/on/host` in `.env` to point elsewhere (e.g. a data
volume on a server). Because it's a bind mount, whatever you download lands on
the host and survives `docker compose down`, image rebuilds, and container
recreation — including the warmed Numba cache (`ref_data/numba_cache/`), which
must NOT be rebuilt on every fresh container (see warm-up step below).

Download the reference data once, after the first `docker compose up`:

```bash
# Inside the running flask_app container:
docker compose exec flask_app bash
cd app/tools/HLA-PepClust/data/ref_data

# Class I human ref data
python3 -m gdown 'https://drive.google.com/uc?id=1iAAvir1woMOnURkP46zr_ETqpW2oUgGD'
unzip Gibbs_motifs_human.zip && rm Gibbs_motifs_human.zip

# Class II human ref data
python3 -m gdown 'https://drive.google.com/uc?id=1Dd9_63vCHezgculDycqfbCDzXsARsvbi'
unzip Gibbs_motifs_human_classII.zip && rm Gibbs_motifs_human_classII.zip

# Build the Class II database
source ../hlapepclust-env/bin/activate
python -m cli.main --database scripts/classII_config.json
deactivate
```

After the first Class II job completes, warm the Numba cache once (or MHC-TP output will be empty):
```bash
# Use the gibbscluster output dir from a completed Class II job:
clust-search /pvol/<taskId>/<sample>/gibbscluster/<replicate>/<subdir>/ \
  data/ref_data -im --output /tmp/cache_warmup \
  --NumbaDB data/ref_data -s human_classii -t 0.1 --processes 4
rm -rf /tmp/cache_warmup
```
