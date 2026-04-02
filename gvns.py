"""
General Variable Neighborhood Search (GVNS) for GrVRP-PCAFS.

References:
- Hansen, P. and Mladenovic, N. (2001). Variable neighborhood search.
  European Journal of Operational Research, 130(3), 449-467.
- Mladenovic, N. and Hansen, P. (1997). Variable neighborhood search.
  Computers & Operations Research, 24(11), 1097-1100.

GVNS Structure:
1. Generate initial solution
2. Repeat until stopping criterion:
   a. Set k = 1
   b. While k <= k_max:
      - Shaking: Generate random neighbor x' in N_k(x)
      - Local search: Apply VND to x' → x''
      - Move or not: If x'' improves x, set x = x'' and k = 1
        Otherwise, k = k + 1
3. Return best solution found
"""

import time
import random
import copy
import numpy as np
from typing import List, Tuple, Optional, Callable
from instance import Instance
from solution import (Solution, Route, compute_total_distance, 
                       evaluate_solution, compute_penalty)
from construction import greedy_construction, savings_construction, insert_station_if_needed
from neighborhoods import (
    SHAKE_NEIGHBORHOODS, LOCAL_SEARCH_NEIGHBORHOODS,
    fix_stations, relocate_best, swap_best, two_opt_intra_best,
)


def solution_cost(solution: Solution) -> float:
    """Compute solution cost (distance + penalty for infeasibility)."""
    dist = compute_total_distance(solution)
    penalty = compute_penalty(solution)
    return dist + 1000 * penalty


def is_feasible(solution: Solution) -> bool:
    """Check if solution is feasible."""
    eval_data = evaluate_solution(solution)
    return eval_data['feasible']


def vnd(solution: Solution, neighborhoods: List[Callable], 
        time_limit: float = None, start_time: float = None) -> Solution:
    """
    Variable Neighborhood Descent (VND).
    
    Apply local search neighborhoods sequentially.
    When improvement found, restart from first neighborhood.
    """
    current = solution
    current_cost = solution_cost(current)
    k = 0
    
    while k < len(neighborhoods):
        if time_limit and start_time:
            if time.time() - start_time >= time_limit:
                break
        
        # Apply local search in neighborhood k
        neighbor = neighborhoods[k](current)
        neighbor_cost = solution_cost(neighbor)
        
        if neighbor_cost < current_cost - 1e-6:
            current = neighbor
            current_cost = neighbor_cost
            k = 0  # Restart from first neighborhood
        else:
            k += 1
    
    return current


def gvns(instance: Instance, 
         time_limit: float = 600.0,
         max_iterations: int = 10000,
         max_no_improve: int = 500,
         k_max: int = 6,
         seed: int = None,
         verbose: bool = True) -> Tuple[Solution, dict]:
    """
    General Variable Neighborhood Search for GrVRP-PCAFS.
    
    Parameters:
    -----------
    instance : Instance
        Problem instance
    time_limit : float
        Maximum time in seconds
    max_iterations : int
        Maximum number of iterations
    max_no_improve : int
        Max iterations without improvement before restart
    k_max : int
        Maximum neighborhood index for shaking
    seed : int
        Random seed for reproducibility
    verbose : bool
        Print progress
    
    Returns:
    --------
    (best_solution, stats): Best solution found and algorithm statistics
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
    start_time = time.time()
    
    # Statistics tracking
    stats = {
        'iterations': 0,
        'improvements': 0,
        'time_to_best': 0.0,
        'best_costs': [],
        'feasible_found': False,
    }
    
    # Generate initial solution
    if verbose:
        print("Generating initial solution...")
    
    best_solution = greedy_construction(instance)
    best_solution = fix_stations(best_solution)
    
    # Also try savings construction
    savings_sol = savings_construction(instance)
    savings_sol = fix_stations(savings_sol)
    
    if solution_cost(savings_sol) < solution_cost(best_solution):
        best_solution = savings_sol
    
    best_cost = solution_cost(best_solution)
    current = best_solution.copy()
    current_cost = best_cost
    
    if verbose:
        eval_data = evaluate_solution(best_solution)
        print(f"Initial solution: distance={eval_data['total_distance']:.2f}, "
              f"routes={eval_data['n_routes']}, feasible={eval_data['feasible']}")
    
    stats['best_costs'].append(best_cost)
    
    no_improve_count = 0
    iteration = 0
    
    shake_neighborhoods = SHAKE_NEIGHBORHOODS[:k_max]
    ls_neighborhoods = LOCAL_SEARCH_NEIGHBORHOODS
    
    while iteration < max_iterations:
        elapsed = time.time() - start_time
        if elapsed >= time_limit:
            break
        
        k = 0
        while k < len(shake_neighborhoods):
            elapsed = time.time() - start_time
            if elapsed >= time_limit:
                break
            
            iteration += 1
            stats['iterations'] = iteration
            
            # ---- SHAKING ----
            neighbor = shake_neighborhoods[k](current)
            
            # Fix stations after shaking
            neighbor = fix_stations(neighbor)
            
            # ---- LOCAL SEARCH (VND) ----
            remaining_time = time_limit - elapsed
            neighbor = vnd(neighbor, ls_neighborhoods, 
                          time_limit=remaining_time, start_time=time.time())
            
            neighbor_cost = solution_cost(neighbor)
            
            # ---- MOVE OR NOT ----
            if neighbor_cost < current_cost - 1e-6:
                current = neighbor
                current_cost = neighbor_cost
                k = 0  # Reset neighborhood
                no_improve_count = 0
                
                # Update best
                if current_cost < best_cost - 1e-6:
                    best_solution = current.copy()
                    best_cost = current_cost
                    stats['improvements'] += 1
                    stats['time_to_best'] = time.time() - start_time
                    
                    if is_feasible(best_solution):
                        stats['feasible_found'] = True
                    
                    if verbose and iteration % 10 == 0:
                        eval_data = evaluate_solution(best_solution)
                        print(f"  Iter {iteration}: cost={best_cost:.2f}, "
                              f"dist={eval_data['total_distance']:.2f}, "
                              f"feasible={eval_data['feasible']}, "
                              f"time={time.time()-start_time:.1f}s")
            else:
                k += 1
                no_improve_count += 1
            
            stats['best_costs'].append(best_cost)
        
        # Restart if stuck
        if no_improve_count >= max_no_improve:
            if verbose:
                print(f"  Restarting at iter {iteration} (no improvement for {no_improve_count} iters)")
            
            # Perturbation: generate new solution and restart
            current = _perturb(best_solution, instance)
            current = fix_stations(current)
            current_cost = solution_cost(current)
            no_improve_count = 0
    
    elapsed = time.time() - start_time
    stats['total_time'] = elapsed
    
    eval_data = evaluate_solution(best_solution)
    stats['final_distance'] = eval_data['total_distance']
    stats['final_feasible'] = eval_data['feasible']
    stats['final_n_routes'] = eval_data['n_routes']
    
    if verbose:
        print(f"\nGVNS completed:")
        print(f"  Total time: {elapsed:.2f}s")
        print(f"  Iterations: {stats['iterations']}")
        print(f"  Improvements: {stats['improvements']}")
        print(f"  Best distance: {stats['final_distance']:.2f}")
        print(f"  Feasible: {stats['final_feasible']}")
        print(f"  Routes: {stats['final_n_routes']}")
        print(f"  Time to best: {stats['time_to_best']:.2f}s")
    
    return best_solution, stats


def _perturb(solution: Solution, instance: Instance) -> Solution:
    """Strong perturbation for restart: randomly remove and reinsert customers."""
    sol = solution.copy()
    
    # Remove 30-50% of customers randomly
    all_customers = list(sol.get_all_customers())
    n_remove = max(1, int(len(all_customers) * random.uniform(0.3, 0.5)))
    to_remove = random.sample(all_customers, min(n_remove, len(all_customers)))
    
    for c in to_remove:
        for route in sol.routes:
            if c in route.nodes:
                route.nodes.remove(c)
                break
    
    sol.remove_empty_routes()
    
    # Reinsert removed customers greedily
    for c in to_remove:
        best_cost_increase = float('inf')
        best_route_idx = -1
        best_pos = -1
        
        for r_idx, route in enumerate(sol.routes):
            for pos in range(1, len(route.nodes)):
                prev = route.nodes[pos - 1]
                next_n = route.nodes[pos]
                cost_increase = (instance.dist(prev, c) + instance.dist(c, next_n) 
                                - instance.dist(prev, next_n))
                
                if cost_increase < best_cost_increase:
                    best_cost_increase = cost_increase
                    best_route_idx = r_idx
                    best_pos = pos
        
        # Also try new route
        if len(sol.routes) < instance.n_vehicles:
            new_route_cost = instance.dist(0, c) + instance.dist(c, 0)
            if new_route_cost < best_cost_increase:
                sol.routes.append(Route([0, c, 0]))
                continue
        
        if best_route_idx >= 0:
            sol.routes[best_route_idx].insert(best_pos, c)
    
    sol.invalidate_cache()
    return sol


def run_multiple(instance: Instance, n_runs: int = 5, 
                 time_limit: float = 600.0, **kwargs) -> Tuple[Solution, dict]:
    """
    Run GVNS multiple times and return best solution.
    
    Returns: (best_solution, aggregated_stats)
    """
    best_solution = None
    best_cost = float('inf')
    all_stats = []
    
    per_run_time = time_limit / n_runs
    
    for run in range(n_runs):
        print(f"\n{'='*60}")
        print(f"Run {run + 1}/{n_runs} (time limit: {per_run_time:.0f}s)")
        print(f"{'='*60}")
        
        sol, stats = gvns(instance, time_limit=per_run_time, 
                         seed=run * 42 + 7, **kwargs)
        all_stats.append(stats)
        
        cost = solution_cost(sol)
        if cost < best_cost:
            best_cost = cost
            best_solution = sol
    
    # Aggregate stats
    agg_stats = {
        'n_runs': n_runs,
        'best_distance': min(s['final_distance'] for s in all_stats),
        'avg_distance': np.mean([s['final_distance'] for s in all_stats]),
        'std_distance': np.std([s['final_distance'] for s in all_stats]),
        'worst_distance': max(s['final_distance'] for s in all_stats),
        'avg_time_to_best': np.mean([s['time_to_best'] for s in all_stats]),
        'avg_iterations': np.mean([s['iterations'] for s in all_stats]),
        'feasible_runs': sum(1 for s in all_stats if s['final_feasible']),
        'all_stats': all_stats,
    }
    
    print(f"\n{'='*60}")
    print(f"Aggregated Results ({n_runs} runs):")
    print(f"  Best distance: {agg_stats['best_distance']:.2f}")
    print(f"  Avg distance:  {agg_stats['avg_distance']:.2f} ± {agg_stats['std_distance']:.2f}")
    print(f"  Worst distance: {agg_stats['worst_distance']:.2f}")
    print(f"  Feasible runs: {agg_stats['feasible_runs']}/{n_runs}")
    print(f"{'='*60}")
    
    return best_solution, agg_stats
