FROM python:3.12-slim-bookworm

WORKDIR /app

# ── Layer 1: OS dependencies ────────────────────────────────────────────────
# apt-get upgrade -y resolves the bulk of SCA/SCD CVE findings from Checkmarx
# curl deliberately excluded — a security scan flagged it as unnecessary
# attack surface on a prior image; wget already covers everything we need,
# including the HEALTHCHECK below.
RUN apt-get update && apt-get upgrade -y && apt-get install -y --no-install-recommends \
    bash git wget unzip tar \
    build-essential libssl-dev libbz2-dev libreadline-dev libsqlite3-dev zlib1g-dev \
    sqlite3 r-base ncompress tcsh perl \
    fonts-urw-base35 gsfonts \
  && rm -rf /var/lib/apt/lists/*

# ── Layer 2: Python 2.7 (required by seq2logo and some legacy tool helpers) ─
# Debian bookworm dropped python2 from main repos — compile from source
RUN wget -q https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz \
  && tar xf Python-2.7.18.tgz \
  && cd Python-2.7.18 \
  && ./configure --prefix=/usr/local/python2 > /dev/null \
  && make -j"$(nproc)" > /dev/null \
  && make install > /dev/null \
  && cd .. && rm -rf Python-2.7.18 Python-2.7.18.tgz \
  && ln -sf /usr/local/python2/bin/python2 /usr/local/bin/python2 \
  && ln -sf /usr/local/python2/bin/python2.7 /usr/local/bin/python2.7

RUN wget -q https://bootstrap.pypa.io/pip/2.7/get-pip.py \
  && python2 get-pip.py > /dev/null && rm get-pip.py \
  && python2 -m pip install --quiet numpy matplotlib

# ── Layer 4: App code ───────────────────────────────────────────────────────
COPY . .

# ── Layer 3b: Ghostscript 9.53.3 pre-built binary ──────────────────────────
# System gs (9.06) has a font rendering bug — must use 9.53.3. Must come
# after COPY: extracts into app/tools/, which doesn't exist until the repo
# is copied in.
RUN wget -q https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs9533/ghostscript-9.53.3-linux-x86_64.tgz \
  && tar -xzf ghostscript-9.53.3-linux-x86_64.tgz -C /app/app/tools/ \
  && rm ghostscript-9.53.3-linux-x86_64.tgz \
  && chmod +x /app/app/tools/ghostscript-9.53.3-linux-x86_64/gs-9533-linux-x86_64

# ── Layer 5: seq2logo and GibbsCluster (tarballs in repo) ──────────────────
RUN tar -xzf app/tools/seq2logo-2.1.all.tar.gz -C app/tools/ \
  && rm app/tools/seq2logo-2.1.all.tar.gz

RUN tar -xvf app/tools/gibbscluster-2.0f.Linux.tar.gz -C app/tools/ \
  && rm app/tools/gibbscluster-2.0f.Linux.tar.gz

# Patch GibbsCluster GIBBS path
RUN sed -i 's|setenv GIBBS .*|setenv GIBBS /app/app/tools/gibbscluster-2.0|' \
    app/tools/gibbscluster-2.0/gibbscluster

# Patch seq2logo to use bundled Ghostscript 9.53.3
RUN sed -i "s|gsPath='gs'|gsPath='/app/app/tools/ghostscript-9.53.3-linux-x86_64/gs-9533-linux-x86_64'|" \
    app/tools/seq2logo-2.1/Seq2Logo.py

# Apply GibbsCluster SA script patches (output path + seqlogo command)
RUN sed -i \
    -e 's|^\(\s*\$resdir .= "/\$prefix";\)|# \1  # removed prefix from output path|' \
    -e 's|^\(my \$barplot = "\$resdir/images/\$prefix.gibbs.KLDvsCluster.barplot.png";\)|my \$barplot = "\$resdir/images/gibbsKLDvsCluster.barplot.JPG";|' \
    -e '530s@.*@$cmd .= "$seq2logo -f $corefile -o $logofile -I 2 --format [JPEG] -b $wlc -C 2 -S 2 -t $title \&>/dev/null; ";@' \
    app/tools/gibbscluster-2.0/GibbsCluster-2.0e_SA.pl

# ── Layer 6: Licensed tools (netMHCpan / netMHCIIpan) ──────────────────────
# IMPORTANT: These tools require a DTU academic license.
# Download from https://services.healthtech.dtu.dk/software.php and place the
# tarballs in app/tools/ BEFORE running docker build.
# If the tarballs are absent the build continues and those prediction methods
# will be unavailable at runtime.
# The vendor scripts use a literal TAB between "setenv" and "NMHOME" (not
# spaces) on the NMHOME line specifically, while other setenv lines use
# spaces — a space-only sed pattern silently fails to match NMHOME, leaving
# it pointed at the vendor's own build path (/tools/src/...) and breaking
# every prediction at runtime with "no binaries found". Match either.
RUN if [ -f app/tools/netMHCpan-4.2c.Linux.tar.gz ]; then \
      gunzip -c app/tools/netMHCpan-4.2c.Linux.tar.gz | tar xf - -C app/tools/ \
      && mkdir -p app/tools/netMHCpan-4.2/tmp \
      && sed -i 's|setenv[ \t]*NMHOME[ \t]*[^ \t]*|setenv NMHOME /app/app/tools/netMHCpan-4.2|g' \
             app/tools/netMHCpan-4.2/netMHCpan \
      && sed -i 's|setenv[ \t]*TMPDIR[ \t]*[^ \t]*|setenv TMPDIR /app/app/tools/netMHCpan-4.2/tmp|g' \
             app/tools/netMHCpan-4.2/netMHCpan \
      && rm app/tools/netMHCpan-4.2c.Linux.tar.gz; \
    else echo "WARNING: netMHCpan-4.2c.Linux.tar.gz not found — netMHCpan will be unavailable"; fi

RUN if [ -f app/tools/netMHCIIpan-4.3j.Linux.tar.gz ]; then \
      tar -xvf app/tools/netMHCIIpan-4.3j.Linux.tar.gz -C app/tools/ \
      && mkdir -p app/tools/netMHCIIpan-4.3/tmp \
      && sed -i 's|setenv[ \t]*NMHOME[ \t]*[^ \t]*|setenv NMHOME /app/app/tools/netMHCIIpan-4.3|g' \
             app/tools/netMHCIIpan-4.3/netMHCIIpan \
      && rm app/tools/netMHCIIpan-4.3j.Linux.tar.gz; \
    else echo "WARNING: netMHCIIpan-4.3j.Linux.tar.gz not found — netMHCIIpan will be unavailable"; fi

# ── Layer 7: MixMHCpred and MixMHC2pred ────────────────────────────────────
RUN wget -q https://github.com/GfellerLab/MixMHCpred/archive/refs/tags/v3.0.tar.gz -O /tmp/mixmhcpred.tar.gz \
  && tar -xzf /tmp/mixmhcpred.tar.gz -C app/tools/ \
  && mv app/tools/MixMHCpred-3.0 app/tools/MixMHCpred \
  && chmod +x app/tools/MixMHCpred/MixMHCpred \
  && rm /tmp/mixmhcpred.tar.gz

# The release zip extracts flat (MixMHC2pred_unix, PWMdef/, README.md, ... at
# archive root, no top-level MixMHC2pred-2.0/ folder) — extract straight into
# the target dir we create, not app/tools/ itself (would collide with our
# own app/tools/README.md).
#
# Using the v2.1 beta1 release, not the stable v2.0.2.2 one: v2.0.2.2 only
# bundles human (HLA) allele PWM definitions. Non-human alleles (e.g. mouse
# H-2) normally come from a separate download at mixmhc2pred.gfellerlab.org/PWMdef,
# but that site is no longer live (domain now resolves to an unrelated host,
# cert doesn't even match). v2.1.beta1 bundles PWMdef for all supported species
# in the same package as its binary — the binary and PWMdef format are
# versioned together (pan_v29), so don't mix v2.0.2.2's binary with v2.1's
# PWMdef data or vice versa. The directory is still called MixMHC2pred-2.0
# because that string is also the tool-name key used to look up allele names
# in app/static/Immunolyser2.0_Allele_Dictionary.csv — keep it as-is.
#
# 2026-11: upstream deleted the original v2.1.beta1 asset (MixMHC2pred-2.1.beta1.zip
# now 404s and, with wget -q, failed every fresh build silently at this layer) and
# re-published it as release v2.1.beta1.2 / MixMHC2pred-2.1-beta1.zip (note -beta1,
# not .beta1). Checked before switching, on all 6,519 MixMHC2pred allele names in the
# allele dictionary: the new binary runs every one, and on the 867 the old build could
# also run, the output (columns and every score) is identical, only the version banner
# differs (vpan_v29 -> v2.1-beta1). The old build lacked definitions for the other
# 5,652 names. The sha256 below pins the asset so upstream cannot swap it unnoticed.
ARG MIXMHC2PRED_URL=https://github.com/GfellerLab/MixMHC2pred/releases/download/v2.1.beta1.2/MixMHC2pred-2.1-beta1.zip
ARG MIXMHC2PRED_SHA256=fe86e0390c96ca7b4e7b8b68d563c717ce43f4091fe460fc283a33ccd716be74
RUN wget -q "${MIXMHC2PRED_URL}" -O /tmp/mixmhc2pred.zip \
  && echo "${MIXMHC2PRED_SHA256}  /tmp/mixmhc2pred.zip" | sha256sum -c - \
  && mkdir -p app/tools/MixMHC2pred-2.0 \
  && unzip -oq /tmp/mixmhc2pred.zip -d app/tools/MixMHC2pred-2.0 \
  && chmod +x app/tools/MixMHC2pred-2.0/MixMHC2pred_unix \
  && rm /tmp/mixmhc2pred.zip

# ── Layer 7b: MoDec (class II motif deconvolution, optional) ────────────────
# Binary-only, free for ACADEMIC NON-COMMERCIAL use under an LICR licence
# (https://github.com/GfellerLab/MoDec/blob/v1.2.1/LICENSE). We must not
# redistribute it, so it is fetched here at build time rather than committed, and
# the licence requires every user of a server running it to have agreed to the
# terms — hence the required tick box on the submit form when MoDec is chosen.
#
# Unlike the licensed DTU tools this download needs no manual step, but a failed
# download (offline build, GitHub outage) must not break the whole image: it only
# warns, and the app then hides the MoDec option (config.MODEC_BIN is absent).
# A CHECKSUM MISMATCH, on the other hand, fails the build — a fetched binary that
# differs from the one we reviewed should never silently end up in the image.
# The Linux binary is byte-identical in v1.2, v1.2.1 and master; the tag is pinned
# so a future upstream change cannot slip in. Update MODEC_SHA256 deliberately.
ARG MODEC_TAG=v1.2.1
ARG MODEC_SHA256=b647bd192b8513f22c78eea1da25b9af5196f1f4b63ccd69a130962f3983031a
RUN mkdir -p app/tools/MoDec \
  && if wget -q "https://raw.githubusercontent.com/GfellerLab/MoDec/${MODEC_TAG}/MoDec_unix" -O app/tools/MoDec/MoDec_unix; then \
       echo "${MODEC_SHA256}  app/tools/MoDec/MoDec_unix" | sha256sum -c - \
       && chmod +x app/tools/MoDec/MoDec_unix; \
     else \
       rm -f app/tools/MoDec/MoDec_unix; \
       echo "WARNING: could not download MoDec ${MODEC_TAG} — the MoDec option will be unavailable"; \
     fi

# ── Layer 8: HLA-PepClust (MHC-TP) ─────────────────────────────────────────
# Clone into HLA-PepClust/ — this directory name is hardcoded in app/utils.py
RUN git clone --depth 1 --branch immunolyser/class2-mhctp \
      https://github.com/PurcellLab/MHC-TP.git app/tools/HLA-PepClust \
  && cd app/tools/HLA-PepClust \
  && python3 -m venv hlapepclust-env \
  && hlapepclust-env/bin/pip install --quiet -e .

# data/ref_data ships committed content (mouse Class I motifs, human/mouse .db,
# warmed numba_cache for those) but at runtime a HOST volume is bind-mounted
# over this exact path (see docker-compose.yml / HLA_PEPCLUST_REF_DATA) so the
# large human/Class II downloads persist across rebuilds. A bind mount hides
# whatever's here, so stash a copy for entrypoint.sh to seed the mount with.
RUN cp -r app/tools/HLA-PepClust/data/ref_data /app/.hlapepclust-ref-data-seed

# ── Layer 9: Python 3 virtualenv + app dependencies ────────────────────────
# mhcflurry still imports pkg_resources, which setuptools 82.0.0 (Feb 2026)
# removed entirely — pin to a version that still has it.
#
# IMPORTANT: do NOT set MHCFLURRY_DOWNLOADS_DIR before this fetch runs.
# Read mhcflurry's own source (mhcflurry/downloads.py configure()): the
# release-lookup logic that `mhcflurry-downloads fetch` depends on
# (get_current_release() / get_downloads_metadata()) only runs inside an
# `if not MHCFLURRY_DOWNLOADS_DIR:` branch — if that var is set at all,
# release detection is skipped entirely and fetch fails with a `KeyError:
# None`, regardless of --release or MHCFLURRY_DOWNLOADS_CURRENT_RELEASE
# (neither is consulted on this path — verified by reading the source
# after both appeared to have no effect). MHCFLURRY_DOWNLOADS_DIR is only
# meant for "data already exists here", not "download fresh data here".
# So: fetch to the default location (as root) first, then relocate.
RUN python3 -m venv lenv \
  && lenv/bin/pip install --quiet --upgrade pip \
  && lenv/bin/pip install --quiet "setuptools==75.6.0" \
  && lenv/bin/pip install --quiet -r requirements_python3.txt \
  && lenv/bin/mhcflurry-downloads fetch

# Downloads fetched as root (above) land in /root/.local/share/mhcflurry/...
# — but the app runs as non-root appuser at runtime (Layer 10) and would
# look in appuser's own home dir instead, finding nothing. Relocate to a
# fixed path under /app (owned by appuser via the chown below) and point
# MHCFLURRY_DOWNLOADS_DIR there for runtime lookups only — safe here since
# by this point the fetch is already done and this env var only affects
# get_path()-based lookups (used by the app), not the release-detection
# logic in `fetch` above.
RUN mkdir -p /app/.mhcflurry \
  && cp -r /root/.local/share/mhcflurry/. /app/.mhcflurry/
ENV MHCFLURRY_DOWNLOADS_DIR=/app/.mhcflurry/4/2.2.0

# Run any hotfix patching needed after package install
RUN lenv/bin/python hotfix_package_files.py 2>/dev/null || true

# ── Layer 10: Non-root user (fixes IaC/CON findings) ───────────────────────
RUN useradd -m --uid 1000 appuser \
  && chown -R appuser:appuser /app

USER appuser

# ── Healthcheck (fixes IaC finding) ─────────────────────────────────────────
HEALTHCHECK --interval=30s --timeout=10s --start-period=90s --retries=3 \
  CMD wget -q --spider http://localhost:5000/healthz || exit 1

EXPOSE 5000

ENTRYPOINT ["./entrypoint.sh"]
