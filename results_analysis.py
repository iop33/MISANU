"""
Analiza rezultata i izvestavanje za GrVRP-PCAFS eksperimente.

Pravi tabele i statistiku za naucni rad.
"""

# =============================================================================
# results_analysis.py  -- IZVESTAJI I TABELE (za naucni rad)
# -----------------------------------------------------------------------------
# Pretvara rezultate u: tekstualnu tabelu, CSV, LaTeX tabelu, .json resenja
# i objedinjeni sazetak po svim setovima. Ovde NEMA algoritma -- samo formatiranje.
# =============================================================================

import os
import json
import numpy as np
import csv
from typing import List, Dict, Any
from datetime import datetime


def format_results_table(results: List[Dict], set_name: str) -> str:
    """
    Formatiraj rezultate kao tabelu spremnu za rad.

    Parametri:
    -----------
    results : lista recnika sa kljucevima:
        - instance_name
        - best_distance
        - avg_distance
        - std_distance
        - worst_distance
        - avg_time
        - n_routes
        - feasible
    """
    # ===== Napravi lepu TEKSTUALNU tabelu (poravnate kolone) za prikaz/snimanje. =====
    header = (f"{'Instance':<25} {'Best':>10} {'Avg±Std':>18} "
              f"{'Worst':>10} {'Routes':>7} {'Time(s)':>8} {'Feas':>5}")
    separator = "-" * len(header)

    lines = [
        f"\n{'='*len(header)}",
        f"Results for {set_name}",
        f"{'='*len(header)}",
        header,
        separator,
    ]

    best_vals = []
    avg_vals = []

    for r in results:                                   # jedan red po instanci
        avg_std = f"{r['avg_distance']:.2f}±{r['std_distance']:.2f}"
        feas = "Yes" if r['feasible'] else "No"
        line = (f"{r['instance_name']:<25} {r['best_distance']:>10.2f} "
                f"{avg_std:>18} {r['worst_distance']:>10.2f} "
                f"{r['n_routes']:>7} {r['avg_time']:>8.1f} {feas:>5}")
        lines.append(line)
        best_vals.append(r['best_distance'])
        avg_vals.append(r['avg_distance'])

    lines.append(separator)

    # Red sa prosecima
    if best_vals:                                       # poslednji red = proseci po setu
        lines.append(
            f"{'Mean':<25} {np.mean(best_vals):>10.2f} "
            f"{np.mean(avg_vals):>10.2f}{'':>8} "
            f"{np.mean([r['worst_distance'] for r in results]):>10.2f} "
            f"{'':>7} {np.mean([r['avg_time'] for r in results]):>8.1f} "
            f"{sum(1 for r in results if r['feasible']):>5}"
        )

    lines.append(f"{'='*len(header)}\n")

    return "\n".join(lines)


def save_results_csv(results: List[Dict], filepath: str):
    """Snimi rezultate u CSV fajl."""
    # ===== Snimi rezultate kao CSV (tabela koju otvara Excel). =====
    if not results:
        return

    keys = ['instance_name', 'best_distance', 'avg_distance', 'std_distance',
            'worst_distance', 'n_routes', 'avg_time', 'feasible']  # kolone CSV-a

    with open(filepath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()                            # red sa nazivima kolona
        for r in results:
            row = {k: r.get(k, '') for k in keys}
            writer.writerow(row)                        # jedan red po instanci


def save_comparison_csv(gvns_results: List[Dict],
                         reference_results: List[Dict],
                         filepath: str):
    """Snimi poredjenje izmedju GVNS-a i referentnog algoritma."""
    # ===== POREDJENJE sa drugim algoritmom (npr. METS/MILP): racuna GAP u %.
    #       (Trenutno se ne poziva iz main.py -- to je sledeci korak.) =====
    keys = ['instance_name', 'gvns_best', 'ref_best', 'gap_%',
            'gvns_avg', 'ref_avg', 'gap_avg_%']

    rows = []
    for g, r in zip(gvns_results, reference_results):   # uparuj nase i referentne rezultate
        gap = 100 * (g['best_distance'] - r['best_distance']) / r['best_distance'] if r['best_distance'] > 0 else 0
        gap_avg = 100 * (g['avg_distance'] - r['avg_distance']) / r['avg_distance'] if r['avg_distance'] > 0 else 0
        rows.append({
            'instance_name': g['instance_name'],
            'gvns_best': g['best_distance'],
            'ref_best': r['best_distance'],
            'gap_%': f"{gap:.2f}",                       # zaostatak nas u odnosu na referencu
            'gvns_avg': g['avg_distance'],
            'ref_avg': r['avg_distance'],
            'gap_avg_%': f"{gap_avg:.2f}",
        })

    with open(filepath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def compute_summary_stats(all_results: Dict[str, List[Dict]]) -> str:
    """Izracunaj i formatiraj objedinjenu statistiku po svim setovima."""
    # ===== Objedinjeni pregled SVIH setova (jedan red po setu). =====
    lines = [
        "\n" + "=" * 70,
        "SUMMARY STATISTICS",
        "=" * 70,
        f"{'Set':<20} {'Instances':>10} {'Avg Best':>10} {'Avg Avg':>10} "
        f"{'Feasible':>10} {'Avg Time':>10}",
        "-" * 70,
    ]

    total_instances = 0
    total_feasible = 0

    for set_name, results in all_results.items():       # za svaki set izracunaj proseke
        n = len(results)
        total_instances += n
        feasible = sum(1 for r in results if r['feasible'])
        total_feasible += feasible

        lines.append(
            f"{set_name:<20} {n:>10} "
            f"{np.mean([r['best_distance'] for r in results]):>10.2f} "
            f"{np.mean([r['avg_distance'] for r in results]):>10.2f} "
            f"{feasible:>10} "
            f"{np.mean([r['avg_time'] for r in results]):>10.1f}"
        )

    lines.append("-" * 70)
    lines.append(f"{'Total':<20} {total_instances:>10} {'':>10} {'':>10} "
                 f"{total_feasible:>10}")
    lines.append("=" * 70)

    return "\n".join(lines)


def save_solution_details(solution, filepath: str):
    """Snimi detaljno resenje u fajl."""
    # ===== Snimi JEDNO resenje u .json: rute + kilometraze + trajanja.
    #       (Ovo pravi one *_solution.json fajlove koje vidis u results/.) =====
    from solution import evaluate_solution

    eval_data = evaluate_solution(solution)             # oceni resenje

    data = {
        'total_distance': eval_data['total_distance'],          # ukupna kilometraza
        'n_routes': eval_data['n_routes'],                      # broj ruta
        'feasible': eval_data['feasible'],                      # da li je izvodljivo
        'routes': [r.nodes for r in solution.routes],           # putevi (nizovi cvorova)
        'route_distances': [e['distance'] for e in eval_data['route_evaluations']],   # km po ruti
        'route_durations': [e['duration'] for e in eval_data['route_evaluations']],   # sati po ruti
    }

    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)


def generate_latex_table(results: List[Dict], set_name: str) -> str:
    """Napravi LaTeX tabelu za rad."""
    # ===== Napravi LaTeX tabelu spremnu za copy-paste u naucni rad. =====
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        f"\\caption{{Results on {set_name} instances}}",
        f"\\label{{tab:{set_name.lower().replace(' ', '_')}}}",
        r"\begin{tabular}{lrrrrr}",
        r"\toprule",
        r"Instance & Best & Avg$\pm$Std & Worst & Routes & Time (s) \\",
        r"\midrule",
    ]

    for r in results:                                   # jedan red po instanci
        line = (f"{r['instance_name']} & {r['best_distance']:.2f} & "
                f"{r['avg_distance']:.2f}$\\pm${r['std_distance']:.2f} & "
                f"{r['worst_distance']:.2f} & {r['n_routes']} & "
                f"{r['avg_time']:.1f} \\\\")
        lines.append(line)

    lines.append(r"\midrule")

    # Red sa prosecima
    if results:                                         # poslednji red = proseci
        lines.append(
            f"Mean & {np.mean([r['best_distance'] for r in results]):.2f} & "
            f"{np.mean([r['avg_distance'] for r in results]):.2f} & "
            f"{np.mean([r['worst_distance'] for r in results]):.2f} & "
            f"-- & {np.mean([r['avg_time'] for r in results]):.1f} \\\\"
        )

    lines.extend([
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ])

    return "\n".join(lines)
