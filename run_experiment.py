"""
Brza skripta za eksperiment (demonstracija i provera).
Pokrece GVNS na reprezentativnom podskupu instanci sa umerenim parametrima.
"""

# =============================================================================
# run_experiment.py  -- KRACI "DEMO" POKRETAC (alternativa za main.py)
# -----------------------------------------------------------------------------
# Manji, brzi eksperiment sa fiksnim parametrima (5 instanci, 3 pokretanja, 20s).
# Sluzi za brzu proveru/demonstraciju, dok je main.py pun eksperiment sa opcijama.
# =============================================================================

import os
import sys
import time
import json
import numpy as np
from datetime import datetime

from instance import Instance
from solution import evaluate_solution
from gvns import gvns
from generate_instances import (
    generate_s_central, generate_m_central, generate_triangle, generate_emh_like
)
from results_analysis import format_results_table, save_results_csv, generate_latex_table


def run_single(instance, n_runs=3, time_per_run=20, seed_base=42):
    """Pokreni GVNS nad jednom instancom, vrati recnik rezultata."""
    # ===== Resi JEDNU instancu: pokreni GVNS n_runs puta, vrati sazetak. =====
    distances = []
    times = []
    n_routes_list = []
    feasible_count = 0
    best_dist = float('inf')

    for run in range(n_runs):                           # vise pokretanja (zbog nasumicnosti)
        sol, stats = gvns(
            instance,
            time_limit=time_per_run,
            max_iterations=20000,
            max_no_improve=100,
            seed=seed_base + run,                       # razlicit seed po pokretanju
            verbose=False,
        )
        ed = evaluate_solution(sol)
        distances.append(ed['total_distance'])
        times.append(stats['time_to_best'])
        n_routes_list.append(ed['n_routes'])
        if ed['feasible']:
            feasible_count += 1
        if ed['total_distance'] < best_dist:
            best_dist = ed['total_distance']            # prati najbolju kilometrazu

    return {                                            # sazetak za instancu
        'instance_name': instance.name,
        'best_distance': min(distances),
        'avg_distance': np.mean(distances),
        'std_distance': np.std(distances),
        'worst_distance': max(distances),
        'n_routes': int(np.median(n_routes_list)),
        'avg_time': np.mean(times),
        'feasible': feasible_count > 0,
        'feasible_count': feasible_count,
    }


def main():
    # ===== Pokreni mali eksperiment nad 5 instanci svakog seta i snimi tabele. =====
    output_dir = "results"
    os.makedirs(output_dir, exist_ok=True)

    all_results = {}
    total_start = time.time()

    # Konfiguracija (fiksni parametri ovog demo pokretanja)
    N_INSTANCES = 5
    N_RUNS = 3
    TIME_PER_RUN = 20  # sekundi

    print(f"GrVRP-PCAFS GVNS Experiment")
    print(f"Runs={N_RUNS}, Time/run={TIME_PER_RUN}s, Instances/set={N_INSTANCES}")
    print("=" * 70)

    # === S-Central ===
    print("\n[S-Central] 15 customers, 1 AFS, η=1")
    results = []
    for i in range(1, N_INSTANCES + 1):
        inst = generate_s_central(i, seed=1000 + i)     # napravi instancu (deterministicno)
        r = run_single(inst, N_RUNS, TIME_PER_RUN)      # resi je
        results.append(r)
        print(f"  {r['instance_name']}: best={r['best_distance']:.2f} "
              f"avg={r['avg_distance']:.2f}±{r['std_distance']:.2f} "
              f"routes={r['n_routes']} feas={'Y' if r['feasible'] else 'N'}")
    all_results['S-Central'] = results

    # === M-Central25 ===
    print("\n[M-Central25] 25 customers, 1 AFS, η=2")
    results = []
    for i in range(1, N_INSTANCES + 1):
        inst = generate_m_central(i, 25, seed=2000 + 2500 + i)
        r = run_single(inst, N_RUNS, TIME_PER_RUN)
        results.append(r)
        print(f"  {r['instance_name']}: best={r['best_distance']:.2f} "
              f"avg={r['avg_distance']:.2f}±{r['std_distance']:.2f} "
              f"routes={r['n_routes']} feas={'Y' if r['feasible'] else 'N'}")
    all_results['M-Central25'] = results

    # === M-Central50 ===
    print("\n[M-Central50] 50 customers, 1 AFS, η=3")
    results = []
    for i in range(1, N_INSTANCES + 1):
        inst = generate_m_central(i, 50, seed=2000 + 5000 + i)
        r = run_single(inst, N_RUNS, TIME_PER_RUN * 2)  # vise vremena za veci set
        results.append(r)
        print(f"  {r['instance_name']}: best={r['best_distance']:.2f} "
              f"avg={r['avg_distance']:.2f}±{r['std_distance']:.2f} "
              f"routes={r['n_routes']} feas={'Y' if r['feasible'] else 'N'}")
    all_results['M-Central50'] = results

    # === Triangle ===
    print("\n[Triangle] 15 customers, 3 AFSs, η=1")
    results = []
    for i in range(1, N_INSTANCES + 1):
        inst = generate_triangle(i, seed=3000 + i)
        r = run_single(inst, N_RUNS, TIME_PER_RUN)
        results.append(r)
        print(f"  {r['instance_name']}: best={r['best_distance']:.2f} "
              f"avg={r['avg_distance']:.2f}±{r['std_distance']:.2f} "
              f"routes={r['n_routes']} feas={'Y' if r['feasible'] else 'N'}")
    all_results['Triangle'] = results

    # === EMH ===
    print("\n[EMH] 20 customers, 6 AFSs, η=1")
    results = []
    for i in range(1, N_INSTANCES + 1):
        inst = generate_emh_like(i, seed=4000 + i)
        r = run_single(inst, N_RUNS, TIME_PER_RUN)
        results.append(r)
        print(f"  {r['instance_name']}: best={r['best_distance']:.2f} "
              f"avg={r['avg_distance']:.2f}±{r['std_distance']:.2f} "
              f"routes={r['n_routes']} feas={'Y' if r['feasible'] else 'N'}")
    all_results['EMH'] = results

    # === Ispisi tabele ===
    print("\n" + "=" * 70)
    for set_name, results in all_results.items():       # ispisi i snimi tabele za svaki set
        table = format_results_table(results, set_name)
        print(table)

        # Snimi
        set_dir = os.path.join(output_dir, set_name)
        os.makedirs(set_dir, exist_ok=True)
        save_results_csv(results, os.path.join(set_dir, "results.csv"))

        latex = generate_latex_table(results, set_name)
        with open(os.path.join(set_dir, "results_latex.tex"), 'w') as f:
            f.write(latex)

    # Sazetak
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"{'Set':<20} {'N':>5} {'Avg Best':>10} {'Avg Mean':>10} {'Feas':>6}")
    print("-" * 55)
    for set_name, results in all_results.items():       # objedinjeni pregled
        n = len(results)
        avg_best = np.mean([r['best_distance'] for r in results])
        avg_mean = np.mean([r['avg_distance'] for r in results])
        feas = sum(1 for r in results if r['feasible'])
        print(f"{set_name:<20} {n:>5} {avg_best:>10.2f} {avg_mean:>10.2f} {feas:>5}/{n}")

    total_time = time.time() - total_start
    print(f"\nTotal time: {total_time:.1f}s")

    # Snimi sve
    serializable = {}                                   # snimi sve rezultate u .json
    for sn, res in all_results.items():
        serializable[sn] = res
    with open(os.path.join(output_dir, "all_results.json"), 'w') as f:
        json.dump(serializable, f, indent=2)

    print(f"Results saved to {output_dir}/")


if __name__ == "__main__":
    main()                                              # pokreni demo eksperiment
