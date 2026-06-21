"""
Full Meltome Atlas scoring run — all 9,160 human proteins.

Features
--------
- Resumes automatically from data/all_9160_scores.csv if interrupted
- Parallel AF2 downloads (8 workers)
- Incremental CSV flush every SAVE_EVERY proteins
- Skips proteins that previously failed AF2/length checks (failed.txt)
- Progress log with ETA to stdout (tee'd to logs/run_all_9160.log by launcher)

Usage
-----
  python3 scripts/run_all_9160.py          # normal run / resume
  python3 scripts/run_all_9160.py --reset  # wipe checkpoint and restart
"""

import csv, json, re, sys, time, traceback
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
import numpy as np

# ── Config ────────────────────────────────────────────────────────────────────
ROOT        = Path(__file__).resolve().parent.parent
OUT_CSV     = ROOT / 'data' / 'all_9160_scores.csv'
FAILED_FILE = ROOT / 'data' / 'all_9160_failed.txt'
TM_SRC      = ROOT / 'data' / 'human_meltome_tms.csv'
TMP_TM      = Path('/tmp/human_meltome_tms.csv')
MIN_LEN     = 150
SAVE_EVERY  = 25
DL_WORKERS  = 8

canonical_re = re.compile(
    r'^([A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]|[OPQ][0-9][A-Z0-9]{3}[0-9])$'
)

# ── Helpers ───────────────────────────────────────────────────────────────────
def log(msg):
    ts = time.strftime('%H:%M:%S')
    print(f'[{ts}] {msg}', flush=True)


def load_pool():
    """Return list of {uniprot, tm, n} dicts for all canonical IDs with n>=3."""
    src = TMP_TM if TMP_TM.exists() else TM_SRC
    parsed = {}
    with open(src) as f:
        for r in csv.DictReader(f):
            pid = r['uniprot_id']
            uid = pid.split('_')[0].split('-')[0]
            if not canonical_re.match(uid):
                continue
            n  = int(r['n_datasets'])
            tm = float(r['median_tm'])
            if n >= 3 and (uid not in parsed or n > parsed[uid]['n']):
                parsed[uid] = {'uniprot': uid, 'tm': tm, 'n': n}
    return list(parsed.values())


def load_done():
    """Return set of UniProt IDs already in the output CSV."""
    done = set()
    if OUT_CSV.exists():
        with open(OUT_CSV) as f:
            for r in csv.DictReader(f):
                done.add(r['uniprot'])
    # also absorb existing pilot runs
    for path in [ROOT / 'data' / 'random_1000_scores.csv',
                 ROOT / 'data' / 'pilot_200_scores.csv']:
        if path.exists():
            with open(path) as f:
                for r in csv.DictReader(f):
                    done.add(r['uniprot'])
    return done


def load_failed():
    """Return set of UniProt IDs that previously failed AF2/length checks."""
    if FAILED_FILE.exists():
        return set(FAILED_FILE.read_text().splitlines())
    return set()


def probe_one(p):
    """Return (uid, pdb_url_or_None, length_or_0)."""
    uid = p['uniprot']
    try:
        r1 = requests.get(
            f'https://alphafold.ebi.ac.uk/api/prediction/{uid}', timeout=20)
        if r1.status_code != 200 or not r1.json():
            return uid, None, 0
        pdb_url = r1.json()[0]['pdbUrl']
        r2 = requests.get(
            f'https://rest.uniprot.org/uniprotkb/{uid}.json',
            params={'fields': 'sequence'}, timeout=20)
        length = r2.json().get('sequence', {}).get('length', 0) if r2.ok else 0
        return uid, pdb_url, length
    except Exception:
        return uid, None, 0


def fetch_pdb(uid, pdb_url, data_dir):
    """Download PDB if not cached; return local path."""
    from prosurf.io.fetch import fetch_af2
    return fetch_af2(uid, data_dir)


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--reset', action='store_true',
                    help='Delete checkpoint and restart from scratch')
    args = ap.parse_args()

    if args.reset:
        OUT_CSV.unlink(missing_ok=True)
        FAILED_FILE.unlink(missing_ok=True)
        log('Checkpoint reset — starting from scratch')

    from prosurf.config import MetricConfig, PathsConfig
    from prosurf.pipeline import analyze_structure

    cfg   = MetricConfig()
    paths = PathsConfig()

    # Load existing scores from prior runs into the checkpoint CSV
    # (so resume works across all score files)
    done    = load_done()
    failed  = load_failed()
    pool    = load_pool()

    log(f'Pool: {len(pool)} proteins | Already done: {len(done)} | '
        f'Previously failed: {len(failed)}')

    todo = [p for p in pool if p['uniprot'] not in done
                             and p['uniprot'] not in failed]
    log(f'Remaining: {len(todo)} proteins to process')

    if not todo:
        log('Nothing left to do — all proteins scored.')
        return

    # If output CSV doesn't exist yet, seed it from existing score files
    if not OUT_CSV.exists():
        log('Seeding checkpoint from existing score files...')
        field_names = ['uniprot', 'tm', 'z_mean', 'z_max', 'z_frac', 'n_patches']
        with open(OUT_CSV, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=field_names)
            w.writeheader()
            for src in [ROOT / 'data' / 'random_1000_scores.csv',
                        ROOT / 'data' / 'pilot_200_scores.csv']:
                if not src.exists():
                    continue
                with open(src) as sf:
                    for r in csv.DictReader(sf):
                        w.writerow({k: r[k] for k in field_names if k in r})
        log(f'Checkpoint seeded with {len(done)} existing scores')

    # ── Process in batches: probe → download → score ──────────────────────────
    BATCH = 200
    total_new  = 0
    total_fail = 0
    new_failed = []
    t0 = time.time()
    field_names = ['uniprot', 'tm', 'z_mean', 'z_max', 'z_frac', 'n_patches']

    for batch_start in range(0, len(todo), BATCH):
        batch = todo[batch_start: batch_start + BATCH]

        # Phase 1: probe AF2 + length in parallel
        probe_results = {}
        with ThreadPoolExecutor(max_workers=DL_WORKERS) as ex:
            futs = {ex.submit(probe_one, p): p for p in batch}
            for fut in as_completed(futs):
                uid, pdb_url, length = fut.result()
                probe_results[uid] = (pdb_url, length)

        confirmed = [(p, probe_results[p['uniprot']][0])
                     for p in batch
                     if probe_results[p['uniprot']][0] is not None
                     and probe_results[p['uniprot']][1] >= MIN_LEN]

        batch_failed = [p['uniprot'] for p in batch
                        if p['uniprot'] not in {c[0]['uniprot'] for c in confirmed}]
        new_failed.extend(batch_failed)
        total_fail += len(batch_failed)

        log(f'Batch {batch_start//BATCH + 1}: '
            f'{len(confirmed)}/{len(batch)} confirmed, '
            f'{len(batch_failed)} failed AF2/length')

        # Phase 2: score each confirmed protein
        buffer = []
        for i, (p, pdb_url) in enumerate(confirmed):
            uid = p['uniprot']
            try:
                pdb_path = fetch_pdb(uid, pdb_url, paths.data_dir)
                _, ps    = analyze_structure(pdb_path, uid, cfg)
                row = {
                    'uniprot':   uid,
                    'tm':        p['tm'],
                    'z_mean':    ps.z_mean,
                    'z_max':     ps.z_max,
                    'z_frac':    ps.z_frac,
                    'n_patches': ps.n_patches,
                }
                buffer.append(row)
                total_new += 1

                # Progress
                done_total = len(done) + total_new
                elapsed    = time.time() - t0
                rate       = total_new / elapsed if elapsed > 0 else 1
                remaining  = len(todo) - total_new - total_fail
                eta_min    = remaining / rate / 60 if rate > 0 else 0
                pct        = done_total / len(pool) * 100

                if total_new % 10 == 0 or i == len(confirmed) - 1:
                    log(f'  [{done_total}/{len(pool)} {pct:.1f}%] '
                        f'{uid}  z_mean={ps.z_mean:.4f}  '
                        f'ETA {eta_min:.0f}min')

            except Exception as e:
                new_failed.append(uid)
                total_fail += 1
                log(f'  FAIL {uid}: {e}')

            # Flush buffer to CSV
            if len(buffer) >= SAVE_EVERY:
                with open(OUT_CSV, 'a', newline='') as f:
                    w = csv.DictWriter(f, fieldnames=field_names)
                    w.writerows(buffer)
                buffer.clear()

        # Flush remainder
        if buffer:
            with open(OUT_CSV, 'a', newline='') as f:
                w = csv.DictWriter(f, fieldnames=field_names)
                w.writerows(buffer)
            buffer.clear()

        # Persist failed list
        if new_failed:
            with open(FAILED_FILE, 'a') as f:
                f.write('\n'.join(new_failed) + '\n')
            new_failed.clear()

    # ── Final summary ─────────────────────────────────────────────────────────
    elapsed_h = (time.time() - t0) / 3600
    log(f'\n{"="*60}')
    log(f'DONE in {elapsed_h:.2f}h')
    log(f'  New proteins scored: {total_new}')
    log(f'  Failures (no AF2 / too short / error): {total_fail}')
    log(f'  Total rows in {OUT_CSV.name}: '
        f'{sum(1 for _ in open(OUT_CSV)) - 1}')

    # Quick Spearman on final combined set
    from scipy.stats import spearmanr
    rows_out = []
    with open(OUT_CSV) as f:
        for r in csv.DictReader(f):
            rows_out.append((float(r['tm']), float(r['z_mean'])))
    tms_all = np.array([x[0] for x in rows_out])
    z_all   = np.array([x[1] for x in rows_out])
    rho, p  = spearmanr(z_all, tms_all)
    log(f'\nFinal Spearman ρ(z_mean, Tm) = {rho:+.4f}   p = {p:.2e}   '
        f'n = {len(rows_out)}')
    log('='*60)


if __name__ == '__main__':
    main()
