"""
Uporedni eksperiment: nas GVNS vs objavljeni METS/GRASP/BKS na PRAVIM instancama.

Ovo je ispravna metodologija (za razliku od starog main.py koji je koristio
IZMISLJENE instance):
    1. ucitaj PRAVE benchmark instance iz METS repozitorijuma (.mat)
    2. pusti nas GVNS n puta po instanci
    3. uporedi nas rezultat sa objavljenim BKS / GRASP / METS (Xu et al. 2025)
    4. izvezi tabelu (tekst + CSV) sa gap-om do BKS i do METS

Koriscenje:
    python run_comparison.py --instances-dir /tmp/METS-Algorithm/Instances \\
        --set S-Central --n-runs 10 --time-limit 30

    --set:        S-Central | M-Central25 | M-Central50 | M-Central100 | all
    --n-runs:     broj nezavisnih pokretanja po instanci (default 10)
    --time-limit: sekundi po jednom pokretanju (default zavisi od velicine seta)
"""

import os
import csv
import time
import argparse
import numpy as np

from mets_loader import load_set
from gvns import gvns
from solution import evaluate_solution
from reference_results import REFERENCE

# prefiks fajla -> (naziv seta, podrazumevano vreme po pokretanju u s)
SETS = {
    'S-Central':    ('15',  30),
    'M-Central25':  ('25',  60),
    'M-Central50':  ('50',  120),
    'M-Central100': ('100', 240),
}


def run_instance(inst, n_runs, time_limit, max_no_improve=400, k_max=6, base_seed=7):
    """Pusti GVNS n_runs puta na jednoj instanci; vrati statistiku po izvodljivim resenjima."""
    feas_dists = []
    all_dists = []
    best_eval = None
    t0 = time.time()
    for run in range(n_runs):
        sol, _ = gvns(inst, time_limit=time_limit, max_no_improve=max_no_improve,
                      k_max=k_max, seed=base_seed + run * 41, verbose=False)
        ev = evaluate_solution(sol)
        all_dists.append(ev['total_distance'])
        if ev['feasible']:
            feas_dists.append(ev['total_distance'])
            if best_eval is None or ev['total_distance'] < best_eval['total_distance']:
                best_eval = ev
    elapsed = time.time() - t0

    return {
        'name': inst.name,
        'feasible_count': len(feas_dists),
        'n_runs': n_runs,
        'best': min(feas_dists) if feas_dists else None,
        'avg': float(np.mean(feas_dists)) if feas_dists else None,
        'std': float(np.std(feas_dists)) if feas_dists else None,
        'worst': max(feas_dists) if feas_dists else None,
        'n_routes': best_eval['n_routes'] if best_eval else None,
        'total_time': elapsed,
    }


def build_rows(results):
    """Spoji nase rezultate sa referentnim (BKS/GRASP/METS) i izracunaj gap-ove."""
    rows = []
    for r in results:
        ref = REFERENCE.get(r['name'])
        bks, grasp, mets = (ref if ref else (None, None, None))
        gap_bks = gap_mets = None
        if r['best'] is not None and bks:
            gap_bks = 100.0 * (r['best'] - bks) / bks
        if r['best'] is not None and mets:
            gap_mets = 100.0 * (r['best'] - mets) / mets
        rows.append({**r, 'bks': bks, 'grasp': grasp, 'mets': mets,
                     'gap_bks': gap_bks, 'gap_mets': gap_mets})
    return rows


def format_table(rows, set_name):
    """Napravi citljivu tekstualnu uporednu tabelu."""
    L = []
    L.append("=" * 104)
    L.append(f"  Uporedni rezultati: NAS GVNS vs METS/GRASP/BKS  |  set: {set_name}")
    L.append("  (PRAVE instance iz Xu et al. 2025; gap>0 znaci da smo losiji)")
    L.append("=" * 104)
    L.append(f"{'Instanca':<16}{'BKS':>9}{'GRASP':>9}{'METS':>9}{'OurBest':>9}{'OurAvg±Std':>16}"
             f"{'gapBKS%':>9}{'gapMETS%':>9}{'feas':>6}")
    L.append("-" * 104)
    g_bks, g_mets = [], []
    for r in rows:
        bks = f"{r['bks']:.2f}" if r['bks'] else "-"
        grasp = f"{r['grasp']:.2f}" if r['grasp'] else "-"
        mets = f"{r['mets']:.2f}" if r['mets'] else "-"
        if r['best'] is None:
            L.append(f"{r['name']:<16}{bks:>9}{grasp:>9}{mets:>9}{'INF':>9}{'-':>16}"
                     f"{'-':>9}{'-':>9}{r['feasible_count']:>3}/{r['n_runs']}")
            continue
        avgstd = f"{r['avg']:.2f}±{r['std']:.2f}"
        gb = f"{r['gap_bks']:+.2f}" if r['gap_bks'] is not None else "-"
        gm = f"{r['gap_mets']:+.2f}" if r['gap_mets'] is not None else "-"
        if r['gap_bks'] is not None: g_bks.append(r['gap_bks'])
        if r['gap_mets'] is not None: g_mets.append(r['gap_mets'])
        L.append(f"{r['name']:<16}{bks:>9}{grasp:>9}{mets:>9}{r['best']:>9.2f}{avgstd:>16}"
                 f"{gb:>9}{gm:>9}{r['feasible_count']:>3}/{r['n_runs']}")
    L.append("-" * 104)
    if g_bks:
        L.append(f"{'PROSEK gap':<16}{'':>27}{'':>9}{'':>16}"
                 f"{np.mean(g_bks):>+9.2f}{np.mean(g_mets):>+9.2f}")
    n_inf = sum(1 for r in rows if r['best'] is None)
    if n_inf:
        L.append(f"  ! {n_inf} instanci bez ijednog izvodljivog resenja (treba popraviti robusnost/vreme)")
    L.append("=" * 104)
    return "\n".join(L)


def save_csv(rows, path):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['instance', 'bks', 'grasp', 'mets', 'our_best', 'our_avg', 'our_std',
                    'our_worst', 'gap_bks_pct', 'gap_mets_pct', 'feasible_count', 'n_runs',
                    'n_routes', 'time_s'])
        for r in rows:
            w.writerow([r['name'], r['bks'], r['grasp'], r['mets'], r['best'], r['avg'],
                        r['std'], r['worst'], r['gap_bks'], r['gap_mets'],
                        r['feasible_count'], r['n_runs'], r['n_routes'], f"{r['total_time']:.1f}"])


def main():
    ap = argparse.ArgumentParser(description="Uporedni eksperiment GVNS vs METS na pravim instancama")
    ap.add_argument('--instances-dir', default='benchmark_instances')
    ap.add_argument('--set', default='S-Central',
                    choices=['S-Central', 'M-Central25', 'M-Central50', 'M-Central100', 'all'])
    ap.add_argument('--n-runs', type=int, default=10)
    ap.add_argument('--time-limit', type=float, default=None, help='s po pokretanju (default po setu)')
    ap.add_argument('--output-dir', default='results_real')
    ap.add_argument('--n-instances', type=int, default=10, help='koliko instanci iz seta (default 10)')
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    set_names = list(SETS.keys()) if args.set == 'all' else [args.set]

    for set_name in set_names:
        prefix, default_tl = SETS[set_name]
        tl = args.time_limit if args.time_limit is not None else default_tl
        instances = load_set(args.instances_dir, prefix)[:args.n_instances]

        print(f"\n>>> {set_name}: {len(instances)} instanci, {args.n_runs} pokretanja x {tl:.0f}s")
        results = []
        for i, inst in enumerate(instances):
            print(f"  [{i+1}/{len(instances)}] {inst.name} ...", flush=True)
            results.append(run_instance(inst, args.n_runs, tl))

        rows = build_rows(results)
        table = format_table(rows, set_name)
        print("\n" + table)

        with open(os.path.join(args.output_dir, f"{set_name}_table.txt"), 'w') as f:
            f.write(table)
        save_csv(rows, os.path.join(args.output_dir, f"{set_name}_results.csv"))
        print(f"\nSnimljeno u {args.output_dir}/{set_name}_*.{{txt,csv}}")


if __name__ == "__main__":
    main()
