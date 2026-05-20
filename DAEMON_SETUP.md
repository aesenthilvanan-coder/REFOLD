# REFOLD Pipeline Daemon — Setup

The daemon runs the full REFOLD pipeline continuously, processing one ClinVar pathogenic
missense variant at a time and injecting results into `PCD_global_atlas.json`.

## Install (macOS LaunchAgent — survives reboots, restarts on crash)

```bash
# Copy the plist to LaunchAgents
cp /Users/TheAES/Desktop/REFOLD/com.syalis.refold.pcd.plist \
   ~/Library/LaunchAgents/com.syalis.refold.pcd.plist

# Load and start immediately
launchctl load ~/Library/LaunchAgents/com.syalis.refold.pcd.plist
```

## Monitor

```bash
# Live log
tail -f /Users/TheAES/Desktop/REFOLD/logs/pipeline_daemon.log

# Check progress
cat /Users/TheAES/Desktop/REFOLD/data/pipeline_state.json | python3 -m json.tool | head -20

# Watch atlas entry count grow
watch -n 10 "cat /Users/TheAES/Desktop/REFOLD/data/results/PCD_global_atlas.json | python3 -c \"import json,sys; d=json.load(sys.stdin); print(d['total_entries'], 'entries')\""
```

## Stop / Restart

```bash
launchctl unload ~/Library/LaunchAgents/com.syalis.refold.pcd.plist
launchctl load ~/Library/LaunchAgents/com.syalis.refold.pcd.plist
```

## Run manually (foreground, for debugging)

```bash
cd /Users/TheAES/Desktop/REFOLD
python3 scripts/pipeline_daemon.py
```

## Architecture

The daemon pulls from a 46-variant initial queue covering:
- Lysosomal Storage Diseases (Gaucher, Tay-Sachs, Fabry, Krabbe...)
- Amino Acid Metabolism (PKU, Tyrosinemia, Galactosemia, MSUD...)
- Cystic Fibrosis (3 CFTR variants)
- Neurological (ALS-SOD1, Parkinson-LRRK2/PINK1, TTR amyloidosis...)
- Cardiac (HCM-MYH7, Long QT, Brugada, Familial Hypercholesterolemia)
- Hematological (Sickle Cell, Hemophilia A, G6PD)
- Oncological (BRCA1, TP53, PTEN, VHL)
- Other (Wilson Disease, Alpha-1 Antitrypsin, OI, Marfan)

Each variant goes through Stage 1 rescue classification → Stage 2 ANM+fpocket
pocket detection → Stage 3 chaperone generation. Results with druggability < 0.70
are logged as skipped. All complete entries are injected into the atlas and the
website's public JSON is synced automatically.

The daemon never stops — when the queue empties it sleeps and waits for new entries.
Add new variants to `INITIAL_QUEUE` in `scripts/pipeline_daemon.py` at any time.
