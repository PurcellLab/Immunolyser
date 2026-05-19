#!/bin/bash
# Immunolyser Production Deployment — security/hardening branch
# Run this on the new prod VM after:
# 1. New 200GB volume created + attached at /dev/vdb
# 2. Ubuntu 22.04 LTS fresh install
# 3. SSH key configured

set -e

echo "=== Step 1: System Dependencies ==="
sudo apt update && sudo apt install -y \
    python3-pip python3-venv redis-server \
    build-essential wget git sqlite3 \
    r-base tcsh ncompress unzip perl \
    python2 fonts-urw-base35 gsfonts

# Python2 pip + deps
python2 -m ensurepip
python2 -m pip install numpy matplotlib

echo "=== Step 2: Clone Repo ==="
sudo mkdir -p /var/www/Immunolyser
sudo chown ubuntu:ubuntu /var/www/Immunolyser
git clone -b security/hardening https://github.com/PurcellLab/Immunolyser.git /var/www/Immunolyser

echo "=== Step 3: Install Bioinformatics Tools ==="
cd /var/www/Immunolyser/app/tools

# seq2logo + gibbscluster
tar -xzf seq2logo-2.1.all.tar.gz
tar -xzf gibbscluster-2.0f.Linux.tar.gz

# netMHCpan-4.2 (you must have downloaded this from DTU)
tar -xzf netMHCpan-4.2c.Linux.tar.gz
mkdir -p netMHCpan-4.2/tmp
chmod +x netMHCpan-4.2/netMHCpan
# EDIT: replace setenv NMHOME paths in netMHCpan script:
sed -i 's|setenv  NMHOME.*|setenv  NMHOME  /var/www/Immunolyser/app/tools/netMHCpan-4.2|g' netMHCpan-4.2/netMHCpan

# netMHCIIpan-4.3 (you must have downloaded this from DTU)
tar -xzf netMHCIIpan-4.3j.Linux.tar.gz
mkdir -p netMHCIIpan-4.3/tmp
chmod +x netMHCIIpan-4.3/netMHCIIpan
sed -i 's|setenv  NMHOME.*|setenv  NMHOME /var/www/Immunolyser/app/tools/netMHCIIpan-4.3|g' netMHCIIpan-4.3/netMHCIIpan

# MixMHCpred
wget https://github.com/GfellerLab/MixMHCpred/archive/refs/tags/v3.0.tar.gz -O MixMHCpred.tar.gz
tar -xzf MixMHCpred.tar.gz
mv MixMHCpred-3.0 MixMHCpred
chmod +x MixMHCpred/MixMHCpred

# MixMHC2pred
wget https://github.com/GfellerLab/MixMHC2pred/archive/refs/tags/v2.0.2.2.tar.gz -O MixMHC2pred.tar.gz
tar -xzf MixMHC2pred.tar.gz
mv MixMHC2pred-2.0.2.2 MixMHC2pred-2.0
chmod +x MixMHC2pred-2.0/MixMHC2pred_unix

# Ghostscript 9.53.3 (not 9.55.0 — font bug)
wget https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs9533/ghostscript-9.53.3-linux-x86_64.tgz
tar -xzf ghostscript-9.53.3-linux-x86_64.tgz
chmod +x ghostscript-9.53.3-linux-x86_64/gs-9533-linux-x86_64
# Point seq2logo to it:
sed -i "s|gsPath='gs'|gsPath='/var/www/Immunolyser/app/tools/ghostscript-9.53.3-linux-x86_64/gs-9533-linux-x86_64'|" \
    /var/www/Immunolyser/app/tools/seq2logo-2.1/Seq2Logo.py

# HLA-PepClust (MUST clone as HLA-PepClust — hardcoded in code)
git clone https://github.com/PurcellLab/MHC-TP.git HLA-PepClust
cd HLA-PepClust
git fetch origin netmhcpan-data-update-2025
git checkout netmhcpan-data-update-2025
# Fix scipy version pin bug
sed -i 's/scipy==1.16.1/scipy>=1.13/' setup.py
# Create venv + install
python3 -m venv hlapepclust-env
source hlapepclust-env/bin/activate
pip install -e .
deactivate
# Download ref_data
python3 -m pip install gdown
mkdir -p data/ref_data
cd data/ref_data
python3 -m gdown 'https://drive.google.com/uc?id=1iAAvir1woMOnURkP46zr_ETqpW2oUgGD'
unzip Gibbs_motifs_human.zip
rm Gibbs_motifs_human.zip

echo "=== Step 4: Python Virtualenv + Dependencies ==="
cd /var/www/Immunolyser
python3 -m venv lenv
source lenv/bin/activate
pip install -r requirements_python3.txt
pip install mhcflurry
mhcflurry-downloads fetch   # ~2GB, takes 10-15 min

echo "=== Step 5: Mount Data Volume ==="
sudo mkdir -p /pvol
sudo mount /dev/vdb /pvol
sudo chown ubuntu:ubuntu /pvol
# Add to fstab for persistence:
echo "/dev/vdb /pvol ext4 defaults,nofail 0 2" | sudo tee -a /etc/fstab

echo "=== Step 6: Create Symlink for Job Output ==="
ln -s /pvol /var/www/Immunolyser/app/static/images

echo "=== Step 7: Environment File ==="
cd /var/www/Immunolyser
cp .env.example .env
# EDIT .env with your values:
# SECRET_KEY=$(python3 -c "import secrets; print(secrets.token_hex(32))")
# IMMUNOLYSER_DATA=/pvol
# IS_DOCKER=false
# DEBUG=False
# EMAIL_ADDRESS=noreply@immunolyser.erc.monash.edu
# BASE_URL=https://immunolyser.erc.monash.edu
# DEMO_TASK_ID=<uuid from demo job, or leave empty>
# DATA_RETENTION_DAYS=30

echo "!!! PAUSE: Edit .env manually with SECRET_KEY and other values !!!"
nano .env

echo "=== Step 8: Initialize Database ==="
source lenv/bin/activate
python3 -c "from app.job_registry import init_job_registry; from app.email_registry import init_email_registry; init_job_registry(); init_email_registry(); print('DB initialised')"
deactivate

echo "=== Step 9: Systemd Services ==="

# Flask app service
sudo tee /etc/systemd/system/immunolyser.service > /dev/null <<EOF
[Unit]
Description=Immunolyser Flask App
After=network.target redis.service

[Service]
User=ubuntu
WorkingDirectory=/var/www/Immunolyser
EnvironmentFile=/var/www/Immunolyser/.env
ExecStart=/var/www/Immunolyser/lenv/bin/gunicorn --workers 2 --bind 0.0.0.0:5000 --timeout 120 firstdemo:app
Restart=always

[Install]
WantedBy=multi-user.target
EOF

# Celery worker service
sudo tee /etc/systemd/system/immunolyser-celery.service > /dev/null <<EOF
[Unit]
Description=Immunolyser Celery Worker
After=network.target redis.service

[Service]
User=ubuntu
WorkingDirectory=/var/www/Immunolyser
EnvironmentFile=/var/www/Immunolyser/.env
ExecStart=/var/www/Immunolyser/lenv/bin/celery -A app.celery worker --loglevel=info --concurrency=1
Restart=always

[Install]
WantedBy=multi-user.target
EOF

sudo systemctl daemon-reload
sudo systemctl enable immunolyser immunolyser-celery
sudo systemctl start immunolyser immunolyser-celery

echo "=== Deployment Complete ==="
echo "Flask app: http://<prod-ip>:5000"
echo "Check status: sudo systemctl status immunolyser immunolyser-celery"
echo "Logs: sudo journalctl -u immunolyser -f"
