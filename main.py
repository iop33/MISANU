"""
Glavni pokretac eksperimenta za GrVRP-PCAFS istrazivanje.

Ova skripta:
1. Generise benchmark instance (S-Central, M-Central, Triangle, EMH)
2. Pokrece GVNS na svim instancama
3. Skuplja rezultate i pravi tabele za rad
4. Snima sve izlaze u results direktorijum

Koriscenje:
    python main.py                    # Pokreni sve eksperimente
    python main.py --quick            # Brzi test
    python main.py --set S-Central    # Pokreni samo S-Central set
"""

# =============================================================================
# main.py  -- POKRETAC EKSPERIMENTA (ulazna tacka programa)
# -----------------------------------------------------------------------------
# Ovo se pokrece sa "python main.py". Redom:
#   1) procita parametre iz komandne linije,
#   2) za svaki set: napravi 10 instanci i resi ih GVNS-om (5 puta po instanci),
#   3) snimi rezultate (tabela, CSV, LaTeX, .json resenja, sazetak).
# =============================================================================

import os
import sys
import time
import json
import argparse                                         # citanje opcija iz komandne linije
import numpy as np
from datetime import datetime

from instance import Instance
from solution import Solution, evaluate_solution, compute_total_distance
from construction import greedy_construction
from neighborhoods import fix_stations
from gvns import gvns, run_multiple, solution_cost
from generate_instances import (
    generate_s_central, generate_m_central, generate_triangle,
    generate_emh_like, generate_all_instances
)
from results_analysis import (
    format_results_table, save_results_csv, compute_summary_stats,
    save_solution_details, generate_latex_table
)


def run_experiment_on_instance(instance: Instance, n_runs: int = 5,
                                time_limit_per_run: float = 120.0,
                                verbose: bool = True) -> dict:
    """
    Pokreni GVNS eksperiment nad jednom instancom.

    Vraca recnik rezultata sa statistikom.
    """
    # ===== Resi JEDNU instancu: pokreni GVNS n_runs puta i skupi statistiku. =====
    distances = []
    times_to_best = []
    all_stats = []
    best_solution = None
    best_dist = float('inf')
    n_routes_list = []
    feasible_count = 0

    for run in range(n_runs):                           # n_runs nezavisnih pokretanja
        if verbose:
            print(f"\n--- {instance.name} - Run {run+1}/{n_runs} ---")

        sol, stats = gvns(                              # pokreni GVNS
            instance,
            time_limit=time_limit_per_run,              # vreme po pokretanju (npr. 120s)
            max_iterations=50000,
            max_no_improve=200,
            k_max=6,
            seed=run * 42 + 7,                          # razlicit seed svaki put (7, 49, 91, ...)
            verbose=verbose,
        )

        eval_data = evaluate_solution(sol)
        dist = eval_data['total_distance']
        distances.append(dist)
        times_to_best.append(stats['time_to_best'])
        n_routes_list.append(eval_data['n_routes'])
        all_stats.append(stats)

        if eval_data['feasible']:
            feasible_count += 1                         # broj izvodljivih pokretanja

        if dist < best_dist:                            # zapamti NAJBOLJE od svih pokretanja
            best_dist = dist
            best_solution = sol

    result = {                                          # sazetak za ovu instancu
        'instance_name': instance.name,
        'best_distance': min(distances),
        'avg_distance': np.mean(distances),
        'std_distance': np.std(distances),
        'worst_distance': max(distances),
        'n_routes': int(np.median(n_routes_list)),
        'avg_time': np.mean(times_to_best),
        'feasible': feasible_count > 0,
        'feasible_count': feasible_count,
        'n_runs': n_runs,
        'best_solution': best_solution,
        'all_distances': distances,
    }

    return result


def run_experiment_set(set_name: str, instances: list,
                        n_runs: int = 5, time_limit: float = 120.0,
                        output_dir: str = "results",
                        verbose: bool = True) -> list:
    """Pokreni eksperimente nad celim setom instanci i snimi rezultate."""
    # ===== Resi CEO set (10 instanci) i snimi tabele/CSV/LaTeX/resenja. =====
    set_dir = os.path.join(output_dir, set_name)        # npr. results/S-Central
    os.makedirs(set_dir, exist_ok=True)

    results = []

    for idx, instance in enumerate(instances):          # za svaku instancu u setu...
        print(f"\n{'#'*70}")
        print(f"# Instance {idx+1}/{len(instances)}: {instance.name}")
        print(f"# Set: {set_name}")
        print(f"{'#'*70}")

        result = run_experiment_on_instance(            # ...resi je
            instance, n_runs=n_runs,
            time_limit_per_run=time_limit,
            verbose=verbose
        )
        results.append(result)

        # Snimi detalje najboljeg resenja
        if result['best_solution'] is not None:         # snimi najbolje resenje te instance u .json
            save_solution_details(
                result['best_solution'],
                os.path.join(set_dir, f"{instance.name}_solution.json")
            )

    # Ispisi i snimi tabelu rezultata
    table = format_results_table(results, set_name)     # lepa tekstualna tabela
    print(table)

    with open(os.path.join(set_dir, "results_table.txt"), 'w') as f:
        f.write(table)

    # Snimi CSV
    save_results_csv(results, os.path.join(set_dir, "results.csv"))  # CSV tabela

    # Snimi LaTeX tabelu
    latex = generate_latex_table(results, set_name)     # LaTeX tabela (za rad)
    with open(os.path.join(set_dir, "results_latex.tex"), 'w') as f:
        f.write(latex)

    return results


def main():
    # ===== GLAVNA FUNKCIJA: citanje parametara + petlja po svim setovima. =====
    parser = argparse.ArgumentParser(
        description="GrVRP-PCAFS Experiment Runner - GVNS Metaheuristic"
    )
    parser.add_argument('--quick', action='store_true',          # brzi test (manji parametri)
                       help='Quick test with reduced parameters')
    parser.add_argument('--set', type=str, default='all',        # koji set pokrenuti
                       choices=['all', 'S-Central', 'M-Central25', 'M-Central50',
                               'M-Central100', 'Triangle', 'EMH'],
                       help='Which instance set to run')
    parser.add_argument('--n-runs', type=int, default=5,         # broj pokretanja po instanci
                       help='Number of independent runs per instance')
    parser.add_argument('--time-limit', type=float, default=120.0,  # vreme po pokretanju (s)
                       help='Time limit per run in seconds')
    parser.add_argument('--n-instances', type=int, default=10,   # broj instanci po setu
                       help='Number of instances per set')
    parser.add_argument('--output-dir', type=str, default='results',  # gde se snimaju rezultati
                       help='Output directory for results')
    parser.add_argument('--verbose', action='store_true', default=True,
                       help='Verbose output')
    parser.add_argument('--quiet', action='store_true',
                       help='Minimal output')

    args = parser.parse_args()                          # procitaj sve opcije

    if args.quiet:
        args.verbose = False

    # Brzi rezim za testiranje
    if args.quick:                                      # --quick -> smanji sve za brzi test
        args.n_runs = 2
        args.time_limit = 30.0
        args.n_instances = 3
        print("=== QUICK TEST MODE ===")

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    # Zapisi konfiguraciju
    config = {                                          # zapis parametara ovog pokretanja
        'timestamp': datetime.now().isoformat(),
        'n_runs': args.n_runs,
        'time_limit_per_run': args.time_limit,
        'n_instances': args.n_instances,
        'set': args.set,
        'algorithm': 'GVNS',
    }
    with open(os.path.join(output_dir, "config.json"), 'w') as f:
        json.dump(config, f, indent=2)                  # snimi results/config.json

    print(f"\n{'='*70}")
    print(f"GrVRP-PCAFS Experiment Runner")
    print(f"Algorithm: General Variable Neighborhood Search (GVNS)")
    print(f"{'='*70}")
    print(f"Configuration:")
    print(f"  Runs per instance: {args.n_runs}")
    print(f"  Time limit per run: {args.time_limit}s")
    print(f"  Instances per set: {args.n_instances}")
    print(f"  Instance set: {args.set}")
    print(f"  Output: {output_dir}/")
    print(f"{'='*70}\n")

    all_results = {}
    total_start = time.time()

    # ---- Generisi i pokreni S-Central ----
    if args.set in ['all', 'S-Central']:                # set S-Central (15 musterija)
        instances = [generate_s_central(i, seed=1000+i)  # napravi 10 instanci (seed 1001..1010)
                    for i in range(1, args.n_instances + 1)]
        results = run_experiment_set(
            'S-Central', instances,
            n_runs=args.n_runs, time_limit=args.time_limit,
            output_dir=output_dir, verbose=args.verbose
        )
        all_results['S-Central'] = results

    # ---- Generisi i pokreni M-Central25 ----
    if args.set in ['all', 'M-Central25']:              # set M-Central25 (25 musterija)
        instances = [generate_m_central(i, n_customers=25, seed=2000+25*100+i)
                    for i in range(1, args.n_instances + 1)]
        results = run_experiment_set(
            'M-Central25', instances,
            n_runs=args.n_runs, time_limit=args.time_limit,
            output_dir=output_dir, verbose=args.verbose
        )
        all_results['M-Central25'] = results

    # ---- Generisi i pokreni M-Central50 ----
    if args.set in ['all', 'M-Central50']:              # set M-Central50 (50 musterija)
        instances = [generate_m_central(i, n_customers=50, seed=2000+50*100+i)
                    for i in range(1, args.n_instances + 1)]
        results = run_experiment_set(
            'M-Central50', instances,
            n_runs=args.n_runs, time_limit=args.time_limit,
            output_dir=output_dir, verbose=args.verbose
        )
        all_results['M-Central50'] = results

    # ---- Generisi i pokreni M-Central100 ----
    if args.set in ['all', 'M-Central100']:             # set M-Central100 (100 musterija, najteze)
        instances = [generate_m_central(i, n_customers=100, seed=2000+100*100+i)
                    for i in range(1, args.n_instances + 1)]
        results = run_experiment_set(
            'M-Central100', instances,
            n_runs=args.n_runs, time_limit=args.time_limit,
            output_dir=output_dir, verbose=args.verbose
        )
        all_results['M-Central100'] = results

    # ---- Generisi i pokreni Triangle ----
    if args.set in ['all', 'Triangle']:                 # set Triangle (3 punionice na sredini)
        instances = [generate_triangle(i, seed=3000+i)
                    for i in range(1, args.n_instances + 1)]
        results = run_experiment_set(
            'Triangle', instances,
            n_runs=args.n_runs, time_limit=args.time_limit,
            output_dir=output_dir, verbose=args.verbose
        )
        all_results['Triangle'] = results

    # ---- Generisi i pokreni EMH ----
    if args.set in ['all', 'EMH']:                      # set EMH (sve izmesano, 6 punionica)
        instances = [generate_emh_like(i, seed=4000+i)
                    for i in range(1, args.n_instances + 1)]
        results = run_experiment_set(
            'EMH', instances,
            n_runs=args.n_runs, time_limit=args.time_limit,
            output_dir=output_dir, verbose=args.verbose
        )
        all_results['EMH'] = results

    # ---- Sazetak ----
    total_time = time.time() - total_start

    summary = compute_summary_stats(all_results)        # ukupan pregled svih setova
    print(summary)

    with open(os.path.join(output_dir, "summary.txt"), 'w') as f:
        f.write(summary)
        f.write(f"\n\nTotal experiment time: {total_time:.1f}s")

    # Snimi sve rezultate kao JSON
    serializable_results = {}                           # spakuj sve rezultate u .json (bez objekata resenja)
    for set_name, results in all_results.items():
        serializable_results[set_name] = [
            {k: v for k, v in r.items()
             if k not in ['best_solution', 'all_distances']}  # izbaci polja nepodesna za json
            for r in results
        ]

    with open(os.path.join(output_dir, "all_results.json"), 'w') as f:
        json.dump(serializable_results, f, indent=2)

    print(f"\nTotal experiment time: {total_time:.1f}s")
    print(f"All results saved to: {output_dir}/")

    return all_results


if __name__ == "__main__":
    main()                                              # "prekidac za paljenje" -- pokrece ceo eksperiment
