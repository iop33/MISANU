"""
Main experiment runner for GrVRP-PCAFS research.

This script:
1. Generates benchmark instances (S-Central, M-Central, Triangle, EMH)
2. Runs GVNS on all instances
3. Collects results and generates tables for the paper
4. Saves all outputs to the results directory

Usage:
    python main.py                    # Run all experiments
    python main.py --quick            # Quick test run
    python main.py --set S-Central    # Run only S-Central set
"""

import os
import sys
import time
import json
import argparse
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
    Run GVNS experiment on a single instance.
    
    Returns result dict with statistics.
    """
    distances = []
    times_to_best = []
    all_stats = []
    best_solution = None
    best_dist = float('inf')
    n_routes_list = []
    feasible_count = 0
    
    for run in range(n_runs):
        if verbose:
            print(f"\n--- {instance.name} - Run {run+1}/{n_runs} ---")
        
        sol, stats = gvns(
            instance, 
            time_limit=time_limit_per_run,
            max_iterations=50000,
            max_no_improve=200,
            k_max=6,
            seed=run * 42 + 7,
            verbose=verbose,
        )
        
        eval_data = evaluate_solution(sol)
        dist = eval_data['total_distance']
        distances.append(dist)
        times_to_best.append(stats['time_to_best'])
        n_routes_list.append(eval_data['n_routes'])
        all_stats.append(stats)
        
        if eval_data['feasible']:
            feasible_count += 1
        
        if dist < best_dist:
            best_dist = dist
            best_solution = sol
    
    result = {
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
    """Run experiments on a set of instances and save results."""
    set_dir = os.path.join(output_dir, set_name)
    os.makedirs(set_dir, exist_ok=True)
    
    results = []
    
    for idx, instance in enumerate(instances):
        print(f"\n{'#'*70}")
        print(f"# Instance {idx+1}/{len(instances)}: {instance.name}")
        print(f"# Set: {set_name}")
        print(f"{'#'*70}")
        
        result = run_experiment_on_instance(
            instance, n_runs=n_runs, 
            time_limit_per_run=time_limit,
            verbose=verbose
        )
        results.append(result)
        
        # Save solution details
        if result['best_solution'] is not None:
            save_solution_details(
                result['best_solution'],
                os.path.join(set_dir, f"{instance.name}_solution.json")
            )
    
    # Print and save results table
    table = format_results_table(results, set_name)
    print(table)
    
    with open(os.path.join(set_dir, "results_table.txt"), 'w') as f:
        f.write(table)
    
    # Save CSV
    save_results_csv(results, os.path.join(set_dir, "results.csv"))
    
    # Save LaTeX table
    latex = generate_latex_table(results, set_name)
    with open(os.path.join(set_dir, "results_latex.tex"), 'w') as f:
        f.write(latex)
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description="GrVRP-PCAFS Experiment Runner - GVNS Metaheuristic"
    )
    parser.add_argument('--quick', action='store_true',
                       help='Quick test with reduced parameters')
    parser.add_argument('--set', type=str, default='all',
                       choices=['all', 'S-Central', 'M-Central25', 'M-Central50', 
                               'M-Central100', 'Triangle', 'EMH'],
                       help='Which instance set to run')
    parser.add_argument('--n-runs', type=int, default=5,
                       help='Number of independent runs per instance')
    parser.add_argument('--time-limit', type=float, default=120.0,
                       help='Time limit per run in seconds')
    parser.add_argument('--n-instances', type=int, default=10,
                       help='Number of instances per set')
    parser.add_argument('--output-dir', type=str, default='results',
                       help='Output directory for results')
    parser.add_argument('--verbose', action='store_true', default=True,
                       help='Verbose output')
    parser.add_argument('--quiet', action='store_true',
                       help='Minimal output')
    
    args = parser.parse_args()
    
    if args.quiet:
        args.verbose = False
    
    # Quick mode for testing
    if args.quick:
        args.n_runs = 2
        args.time_limit = 30.0
        args.n_instances = 3
        print("=== QUICK TEST MODE ===")
    
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    # Log configuration
    config = {
        'timestamp': datetime.now().isoformat(),
        'n_runs': args.n_runs,
        'time_limit_per_run': args.time_limit,
        'n_instances': args.n_instances,
        'set': args.set,
        'algorithm': 'GVNS',
    }
    with open(os.path.join(output_dir, "config.json"), 'w') as f:
        json.dump(config, f, indent=2)
    
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
    
    # ---- Generate and run S-Central ----
    if args.set in ['all', 'S-Central']:
        instances = [generate_s_central(i, seed=1000+i) 
                    for i in range(1, args.n_instances + 1)]
        results = run_experiment_set(
            'S-Central', instances, 
            n_runs=args.n_runs, time_limit=args.time_limit,
            output_dir=output_dir, verbose=args.verbose
        )
        all_results['S-Central'] = results
    
    # ---- Generate and run M-Central25 ----
    if args.set in ['all', 'M-Central25']:
        instances = [generate_m_central(i, n_customers=25, seed=2000+25*100+i) 
                    for i in range(1, args.n_instances + 1)]
        results = run_experiment_set(
            'M-Central25', instances,
            n_runs=args.n_runs, time_limit=args.time_limit,
            output_dir=output_dir, verbose=args.verbose
        )
        all_results['M-Central25'] = results
    
    # ---- Generate and run M-Central50 ----
    if args.set in ['all', 'M-Central50']:
        instances = [generate_m_central(i, n_customers=50, seed=2000+50*100+i) 
                    for i in range(1, args.n_instances + 1)]
        results = run_experiment_set(
            'M-Central50', instances,
            n_runs=args.n_runs, time_limit=args.time_limit,
            output_dir=output_dir, verbose=args.verbose
        )
        all_results['M-Central50'] = results
    
    # ---- Generate and run M-Central100 ----
    if args.set in ['all', 'M-Central100']:
        instances = [generate_m_central(i, n_customers=100, seed=2000+100*100+i) 
                    for i in range(1, args.n_instances + 1)]
        results = run_experiment_set(
            'M-Central100', instances,
            n_runs=args.n_runs, time_limit=args.time_limit,
            output_dir=output_dir, verbose=args.verbose
        )
        all_results['M-Central100'] = results
    
    # ---- Generate and run Triangle ----
    if args.set in ['all', 'Triangle']:
        instances = [generate_triangle(i, seed=3000+i) 
                    for i in range(1, args.n_instances + 1)]
        results = run_experiment_set(
            'Triangle', instances,
            n_runs=args.n_runs, time_limit=args.time_limit,
            output_dir=output_dir, verbose=args.verbose
        )
        all_results['Triangle'] = results
    
    # ---- Generate and run EMH ----
    if args.set in ['all', 'EMH']:
        instances = [generate_emh_like(i, seed=4000+i) 
                    for i in range(1, args.n_instances + 1)]
        results = run_experiment_set(
            'EMH', instances,
            n_runs=args.n_runs, time_limit=args.time_limit,
            output_dir=output_dir, verbose=args.verbose
        )
        all_results['EMH'] = results
    
    # ---- Summary ----
    total_time = time.time() - total_start
    
    summary = compute_summary_stats(all_results)
    print(summary)
    
    with open(os.path.join(output_dir, "summary.txt"), 'w') as f:
        f.write(summary)
        f.write(f"\n\nTotal experiment time: {total_time:.1f}s")
    
    # Save all results as JSON
    serializable_results = {}
    for set_name, results in all_results.items():
        serializable_results[set_name] = [
            {k: v for k, v in r.items() 
             if k not in ['best_solution', 'all_distances']}
            for r in results
        ]
    
    with open(os.path.join(output_dir, "all_results.json"), 'w') as f:
        json.dump(serializable_results, f, indent=2)
    
    print(f"\nTotal experiment time: {total_time:.1f}s")
    print(f"All results saved to: {output_dir}/")
    
    return all_results


if __name__ == "__main__":
    main()
