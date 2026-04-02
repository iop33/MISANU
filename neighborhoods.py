"""
Neighborhood structures for GVNS metaheuristic applied to GrVRP-PCAFS.

Neighborhoods (ordered by complexity):
1. Relocate: Move one customer to another position
2. Swap: Exchange two customers between routes
3. Or-opt: Move a sequence of 2-3 consecutive customers
4. 2-opt intra-route: Reverse a segment within a route
5. 2-opt* inter-route: Exchange tails between two routes
6. Station relocate: Change AFS visit position or station
7. Station insert/remove: Add or remove AFS visits

Each neighborhood has:
- shake(): Generate random neighbor (for diversification)
- local_search(): Find best neighbor (for intensification)
"""

import random
import copy
import math
from typing import List, Tuple, Optional
from instance import Instance
from solution import (Route, Solution, evaluate_route, compute_total_distance, 
                       evaluate_solution, compute_penalty)
from construction import find_nearest_feasible_station, insert_station_if_needed


def _route_distance(route: Route, instance: Instance) -> float:
    """Quick distance computation for a route."""
    d = 0.0
    for i in range(len(route.nodes) - 1):
        d += instance.dist(route.nodes[i], route.nodes[i + 1])
    return d


def _check_fuel_feasibility(route: Route, instance: Instance) -> bool:
    """Quick fuel feasibility check."""
    fuel = instance.tank_capacity
    for i in range(1, len(route.nodes)):
        d = instance.dist(route.nodes[i-1], route.nodes[i])
        fuel -= instance.consumption_rate * d
        if fuel < -1e-6:
            return False
        if instance.is_station(route.nodes[i]):
            fuel = instance.tank_capacity
    return True


def _check_duration_feasibility(route: Route, instance: Instance) -> bool:
    """Quick duration feasibility check."""
    t = instance.p_start
    for i in range(1, len(route.nodes)):
        t += instance.travel_time(route.nodes[i-1], route.nodes[i])
        if instance.is_customer(route.nodes[i]):
            t += instance.service_time_customer
        elif instance.is_station(route.nodes[i]):
            t += instance.refueling_time
    return t <= instance.t_max + 1e-6


def _is_route_feasible(route: Route, instance: Instance) -> bool:
    """Check if a route is feasible (fuel + duration)."""
    return _check_fuel_feasibility(route, instance) and _check_duration_feasibility(route, instance)


# ==================== NEIGHBORHOOD 1: RELOCATE ====================

def relocate_shake(solution: Solution) -> Solution:
    """Randomly relocate a customer to a different position."""
    sol = solution.copy()
    instance = sol.instance
    
    # Find all customer positions
    positions = []
    for r_idx, route in enumerate(sol.routes):
        for pos in range(1, len(route.nodes) - 1):
            if instance.is_customer(route.nodes[pos]):
                positions.append((r_idx, pos))
    
    if not positions:
        return sol
    
    # Pick random customer
    r_idx, pos = random.choice(positions)
    customer = sol.routes[r_idx].nodes[pos]
    sol.routes[r_idx].remove_at(pos)
    
    # Remove empty stations that are no longer needed
    _clean_stations(sol.routes[r_idx], instance)
    
    # Pick random insertion position (could be same or different route)
    possible_routes = list(range(len(sol.routes)))
    if len(sol.routes) < instance.n_vehicles:
        possible_routes.append(-1)  # new route option
    
    target_r = random.choice(possible_routes)
    
    if target_r == -1:
        # New route
        new_route = Route([0, customer, 0])
        new_route = insert_station_if_needed(new_route, instance)
        sol.routes.append(new_route)
    else:
        route = sol.routes[target_r]
        # Random position (between depot nodes)
        insert_pos = random.randint(1, len(route.nodes) - 1)
        route.insert(insert_pos, customer)
    
    sol.remove_empty_routes()
    sol.invalidate_cache()
    return sol


def relocate_best(solution: Solution) -> Solution:
    """Find best relocate move."""
    instance = solution.instance
    best_sol = solution
    best_dist = compute_total_distance(solution)
    best_penalty = compute_penalty(solution)
    best_cost = best_dist + 1000 * best_penalty
    
    # Try all customer relocations
    for r_idx, route in enumerate(solution.routes):
        for pos in range(1, len(route.nodes) - 1):
            if not instance.is_customer(route.nodes[pos]):
                continue
            
            customer = route.nodes[pos]
            
            # Try all insertion positions
            for t_idx, t_route in enumerate(solution.routes):
                for t_pos in range(1, len(t_route.nodes)):
                    if t_idx == r_idx and (t_pos == pos or t_pos == pos + 1):
                        continue
                    
                    # Create candidate
                    sol = solution.copy()
                    sol.routes[r_idx].remove_at(pos)
                    
                    # Adjust target position if same route and pos shifted
                    actual_t_pos = t_pos
                    if t_idx == r_idx and t_pos > pos:
                        actual_t_pos -= 1
                    
                    sol.routes[t_idx].insert(actual_t_pos, customer)
                    sol.remove_empty_routes()
                    
                    dist = compute_total_distance(sol)
                    penalty = compute_penalty(sol)
                    cost = dist + 1000 * penalty
                    
                    if cost < best_cost - 1e-6:
                        best_cost = cost
                        best_sol = sol
    
    best_sol.invalidate_cache()
    return best_sol


# ==================== NEIGHBORHOOD 2: SWAP ====================

def swap_shake(solution: Solution) -> Solution:
    """Randomly swap two customers."""
    sol = solution.copy()
    instance = sol.instance
    
    positions = []
    for r_idx, route in enumerate(sol.routes):
        for pos in range(1, len(route.nodes) - 1):
            if instance.is_customer(route.nodes[pos]):
                positions.append((r_idx, pos))
    
    if len(positions) < 2:
        return sol
    
    (r1, p1), (r2, p2) = random.sample(positions, 2)
    
    # Swap
    sol.routes[r1].nodes[p1], sol.routes[r2].nodes[p2] = \
        sol.routes[r2].nodes[p2], sol.routes[r1].nodes[p1]
    
    sol.invalidate_cache()
    return sol


def swap_best(solution: Solution) -> Solution:
    """Find best swap move."""
    instance = solution.instance
    best_sol = solution
    best_cost = compute_total_distance(solution) + 1000 * compute_penalty(solution)
    
    positions = []
    for r_idx, route in enumerate(solution.routes):
        for pos in range(1, len(route.nodes) - 1):
            if instance.is_customer(route.nodes[pos]):
                positions.append((r_idx, pos))
    
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            r1, p1 = positions[i]
            r2, p2 = positions[j]
            
            sol = solution.copy()
            sol.routes[r1].nodes[p1], sol.routes[r2].nodes[p2] = \
                sol.routes[r2].nodes[p2], sol.routes[r1].nodes[p1]
            
            dist = compute_total_distance(sol)
            penalty = compute_penalty(sol)
            cost = dist + 1000 * penalty
            
            if cost < best_cost - 1e-6:
                best_cost = cost
                best_sol = sol
    
    best_sol.invalidate_cache()
    return best_sol


# ==================== NEIGHBORHOOD 3: OR-OPT ====================

def or_opt_shake(solution: Solution) -> Solution:
    """Randomly move a segment of 1-3 customers to another position."""
    sol = solution.copy()
    instance = sol.instance
    
    # Pick a random route with customers
    routes_with_customers = [
        (r_idx, r) for r_idx, r in enumerate(sol.routes) 
        if len(r.customers(instance)) >= 1
    ]
    if not routes_with_customers:
        return sol
    
    r_idx, route = random.choice(routes_with_customers)
    
    # Find customer positions
    cust_positions = [
        p for p in range(1, len(route.nodes) - 1) 
        if instance.is_customer(route.nodes[p])
    ]
    if not cust_positions:
        return sol
    
    seg_len = random.choice([1, 2, 3])
    seg_len = min(seg_len, len(cust_positions))
    
    start_idx = random.randint(0, len(cust_positions) - seg_len)
    segment_positions = cust_positions[start_idx:start_idx + seg_len]
    
    # Extract segment nodes
    segment = [route.nodes[p] for p in segment_positions]
    
    # Remove segment from route (reverse order to maintain indices)
    for p in sorted(segment_positions, reverse=True):
        sol.routes[r_idx].remove_at(p)
    
    # Insert segment at random position in random route
    possible_routes = list(range(len(sol.routes)))
    target_r = random.choice(possible_routes)
    target_route = sol.routes[target_r]
    insert_pos = random.randint(1, max(1, len(target_route.nodes) - 1))
    
    for k, node in enumerate(segment):
        target_route.insert(insert_pos + k, node)
    
    sol.remove_empty_routes()
    sol.invalidate_cache()
    return sol


# ==================== NEIGHBORHOOD 4: 2-OPT INTRA-ROUTE ====================

def two_opt_intra_shake(solution: Solution) -> Solution:
    """Randomly apply 2-opt within a route (reverse a segment)."""
    sol = solution.copy()
    instance = sol.instance
    
    non_empty = [r_idx for r_idx, r in enumerate(sol.routes) 
                 if len(r.nodes) > 3]
    if not non_empty:
        return sol
    
    r_idx = random.choice(non_empty)
    route = sol.routes[r_idx]
    
    # Pick two positions (not depot)
    n = len(route.nodes)
    if n <= 3:
        return sol
    
    i = random.randint(1, n - 3)
    j = random.randint(i + 1, n - 2)
    
    # Reverse segment between i and j
    route.nodes[i:j+1] = route.nodes[i:j+1][::-1]
    
    sol.invalidate_cache()
    return sol


def two_opt_intra_best(solution: Solution) -> Solution:
    """Find best 2-opt intra-route move."""
    instance = solution.instance
    best_sol = solution
    best_cost = compute_total_distance(solution) + 1000 * compute_penalty(solution)
    
    for r_idx, route in enumerate(solution.routes):
        n = len(route.nodes)
        if n <= 3:
            continue
        
        for i in range(1, n - 2):
            for j in range(i + 1, n - 1):
                sol = solution.copy()
                sol.routes[r_idx].nodes[i:j+1] = sol.routes[r_idx].nodes[i:j+1][::-1]
                
                dist = compute_total_distance(sol)
                penalty = compute_penalty(sol)
                cost = dist + 1000 * penalty
                
                if cost < best_cost - 1e-6:
                    best_cost = cost
                    best_sol = sol
    
    best_sol.invalidate_cache()
    return best_sol


# ==================== NEIGHBORHOOD 5: 2-OPT* INTER-ROUTE ====================

def two_opt_star_shake(solution: Solution) -> Solution:
    """Randomly exchange tails between two routes."""
    sol = solution.copy()
    instance = sol.instance
    
    if len(sol.routes) < 2:
        return sol
    
    r1_idx, r2_idx = random.sample(range(len(sol.routes)), 2)
    route1 = sol.routes[r1_idx]
    route2 = sol.routes[r2_idx]
    
    if len(route1.nodes) <= 2 or len(route2.nodes) <= 2:
        return sol
    
    # Pick cut points
    cut1 = random.randint(1, len(route1.nodes) - 2)
    cut2 = random.randint(1, len(route2.nodes) - 2)
    
    # Exchange tails
    tail1 = route1.nodes[cut1:]
    tail2 = route2.nodes[cut2:]
    
    sol.routes[r1_idx].nodes = route1.nodes[:cut1] + tail2
    sol.routes[r2_idx].nodes = route2.nodes[:cut2] + tail1
    
    sol.remove_empty_routes()
    sol.invalidate_cache()
    return sol


# ==================== NEIGHBORHOOD 6: STATION INSERT/REMOVE ====================

def station_change_shake(solution: Solution) -> Solution:
    """Randomly insert, remove, or change a station visit."""
    sol = solution.copy()
    instance = sol.instance
    
    action = random.choice(['insert', 'remove', 'change'])
    
    if action == 'insert':
        # Insert a random station at a random position
        non_empty = [r_idx for r_idx, r in enumerate(sol.routes) if len(r.nodes) > 2]
        if not non_empty:
            return sol
        r_idx = random.choice(non_empty)
        route = sol.routes[r_idx]
        station = random.choice(instance.station_indices)
        pos = random.randint(1, len(route.nodes) - 1)
        route.insert(pos, station)
    
    elif action == 'remove':
        # Remove a random station
        station_positions = []
        for r_idx, route in enumerate(sol.routes):
            for pos in range(1, len(route.nodes) - 1):
                if instance.is_station(route.nodes[pos]):
                    station_positions.append((r_idx, pos))
        
        if station_positions:
            r_idx, pos = random.choice(station_positions)
            sol.routes[r_idx].remove_at(pos)
    
    elif action == 'change':
        # Change a station to a different one
        station_positions = []
        for r_idx, route in enumerate(sol.routes):
            for pos in range(1, len(route.nodes) - 1):
                if instance.is_station(route.nodes[pos]):
                    station_positions.append((r_idx, pos))
        
        if station_positions:
            r_idx, pos = random.choice(station_positions)
            current_station = sol.routes[r_idx].nodes[pos]
            other_stations = [s for s in instance.station_indices if s != current_station]
            if other_stations:
                new_station = random.choice(other_stations)
                sol.routes[r_idx].nodes[pos] = new_station
    
    sol.invalidate_cache()
    return sol


# ==================== UTILITY ====================

def _clean_stations(route: Route, instance: Instance):
    """Remove unnecessary station visits from a route."""
    changed = True
    while changed:
        changed = False
        for i in range(len(route.nodes) - 2, 0, -1):
            if instance.is_station(route.nodes[i]):
                # Try removing this station
                test_nodes = route.nodes[:i] + route.nodes[i+1:]
                test_route = Route(test_nodes)
                if _check_fuel_feasibility(test_route, instance):
                    route.nodes = test_nodes
                    changed = True
                    break


def fix_stations(solution: Solution) -> Solution:
    """Fix station visits in all routes: remove unnecessary, add where needed."""
    sol = solution.copy()
    instance = sol.instance
    
    for r_idx in range(len(sol.routes)):
        route = sol.routes[r_idx]
        
        # First remove all stations
        customer_nodes = [0] + [n for n in route.nodes[1:-1] if instance.is_customer(n)] + [0]
        route = Route(customer_nodes)
        
        # Then insert stations where needed
        route = insert_station_if_needed(route, instance)
        sol.routes[r_idx] = route
    
    sol.invalidate_cache()
    return sol


# ==================== NEIGHBORHOOD COLLECTIONS ====================

# Shaking neighborhoods (for GVNS diversification)
SHAKE_NEIGHBORHOODS = [
    relocate_shake,
    swap_shake,
    or_opt_shake,
    two_opt_intra_shake,
    two_opt_star_shake,
    station_change_shake,
]

# Local search neighborhoods (for GVNS intensification)
LOCAL_SEARCH_NEIGHBORHOODS = [
    relocate_best,
    swap_best,
    two_opt_intra_best,
]
