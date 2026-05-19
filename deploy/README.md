# Deployment

## Install systemd services (one-time)

```bash
sudo cp deploy/immunolyser.service /etc/systemd/system/
sudo cp deploy/immunolyser-celery.service /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl enable immunolyser immunolyser-celery
sudo systemctl start immunolyser immunolyser-celery
```

## Redeploy after code changes

```bash
cd /var/www/Immunolyser
git pull
sudo systemctl restart immunolyser immunolyser-celery
```

## Check status / logs

```bash
sudo systemctl status immunolyser
sudo journalctl -u immunolyser -f
sudo journalctl -u immunolyser-celery -f
```
