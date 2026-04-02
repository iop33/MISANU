"""
Solution representation and evaluation for GrVRP-PCAFS.

A solution is a set of routes. Each route is a list of node indices
starting and ending at depot (0). Routes may contain AFS visits.

Includes the Reschedule procedure from Bruglieri et al. (2022) / Xu et al. (2025)
for managing limited AFS capacity.
"""

import copy
import math
import numpy as np
from typing import List, Tuple, Optional, Dict
from instance import Instance


class Route:
    """A single vehicle route: sequence of nodes starting and ending at depot."""
    
    def __init__(self, nodes: List[int] = None):
        """
        nodes: list of node indices, e.g. [0, 3, 5, 11, 2, 0]
        where 0=depot, 3,5,2=customers, 11=AFS
        """
        if nodes is None:
            self.nodes = [0, 0]  # empty route
        else:
            self.nodes = list(nodes)
    
    @property
    def is_empty(self) -> bool:
        """Route is empty if it only contains depot nodes."""
        return len(self.nodes) <= 2
    
    def customers(self, instance: Instance) -> List[int]:
        """Return list of customer nodes in this route."""
        return [n for n in self.nodes if instance.is_customer(n)]
    
    def stations(self, instance: Instance) -> List[int]:
        """Return list of AFS nodes in this route."""
        return [n for n in self.nodes if instance.is_station(n)]
    
    def insert(self, position: int, node: int):
        """Insert a node at given position (1-indexed within route body)."""
        self.nodes.insert(position, node)
    
    def remove_node(self, node: int):
        """Remove first occurrence of node from route."""
        self.nodes.remove(node)
    
    def remove_at(self, position: int):
        """Remove node at position."""
        self.nodes.pop(position)
    
    def copy(self) -> 'Route':
        return Route(list(self.nodes))
    
    def __repr__(self):
        return f"Route({self.nodes})"
    
    def __len__(self):
        return len(self.nodes)


class Solution:
    """A complete solution: set of routes covering all customers."""
    
    def __init__(self, instance: Instance, routes: List[Route] = None):
        self.instance = instance
        self.routes = routes if routes is not None else []
        self._total_distance = None
        self._feasible = None
        self._evaluation = None
    
    def invalidate_cache(self):
        self._total_distance = None
        self._feasible = None
        self._evaluation = None
    
    @property
    def n_routes(self) -> int:
        return len(self.routes)
    
    def add_route(self, route: Route):
        self.routes.append(route)
        self.invalidate_cache()
    
    def remove_empty_routes(self):
        self.routes = [r for r in self.routes if not r.is_empty]
        self.invalidate_cache()
    
    def copy(self) -> 'Solution':
        new_sol = Solution(self.instance, [r.copy() for r in self.routes])
        return new_sol
    
    def get_all_customers(self) -> set:
        """Get set of all customers in the solution."""
        customers = set()
        for route in self.routes:
            for n in route.nodes:
                if self.instance.is_customer(n):
                    customers.add(n)
        return customers
    
    def get_unserved_customers(self) -> set:
        """Get customers not yet in any route."""
        served = self.get_all_customers()
        all_customers = set(self.instance.customer_indices)
        return all_customers - served
    
    def __repr__(self):
        eval_info = self.evaluate()
        return (f"Solution(routes={self.n_routes}, "
                f"dist={eval_info['total_distance']:.2f}, "
                f"feasible={eval_info['feasible']})")


def evaluate_route(route: Route, instance: Instance) -> dict:
    """
    Evaluate a single route: compute distance, duration, fuel feasibility.
    Does NOT check AFS capacity (that's a global constraint).
    
    Returns dict with:
    - distance: total travel distance
    - duration: total route duration (travel + service + refueling)
    - fuel_feasible: True if fuel constraints are satisfied
    - duration_feasible: True if duration <= Tmax
    - arrival_times: list of arrival times at each node
    - fuel_levels: list of fuel levels at arrival at each node
    """
    nodes = route.nodes
    n = len(nodes)
    
    if n <= 2:
        return {
            'distance': 0.0,
            'duration': 0.0,
            'fuel_feasible': True,
            'duration_feasible': True,
            'arrival_times': [0.0, 0.0],
            'fuel_levels': [instance.tank_capacity, instance.tank_capacity],
            'departure_times': [0.0, 0.0],
        }
    
    total_distance = 0.0
    arrival_times = [0.0] * n
    departure_times = [0.0] * n
    fuel_levels = [0.0] * n
    
    # Start at depot with full tank
    fuel_levels[0] = instance.tank_capacity
    departure_times[0] = instance.p_start  # initial refuel time
    
    fuel_feasible = True
    
    for i in range(1, n):
        prev = nodes[i - 1]
        curr = nodes[i]
        
        d = instance.dist(prev, curr)
        t = instance.travel_time(prev, curr)
        total_distance += d
        
        # Arrival time
        arrival_times[i] = departure_times[i - 1] + t
        
        # Fuel level on arrival
        fuel_consumed = instance.consumption_rate * d
        fuel_levels[i] = fuel_levels[i - 1] - fuel_consumed
        
        if fuel_levels[i] < -1e-6:
            fuel_feasible = False
        
        # Departure time
        if instance.is_customer(curr):
            departure_times[i] = arrival_times[i] + instance.get_service_time(curr)
        elif instance.is_station(curr):
            # Full refuel
            departure_times[i] = arrival_times[i] + instance.refueling_time
            fuel_levels[i] = instance.tank_capacity  # refuel to full
        elif instance.is_depot(curr):
            departure_times[i] = arrival_times[i]
        else:
            departure_times[i] = arrival_times[i]
    
    duration = arrival_times[-1]  # arrival time at final depot
    
    return {
        'distance': total_distance,
        'duration': duration,
        'fuel_feasible': fuel_feasible,
        'duration_feasible': duration <= instance.t_max + 1e-6,
        'arrival_times': arrival_times,
        'departure_times': departure_times,
        'fuel_levels': fuel_levels,
    }


def reschedule(instance: Instance, solution: 'Solution') -> Tuple[bool, Dict[int, float]]:
    """
    Reschedule procedure: check and resolve AFS capacity conflicts.
    
    Based on Algorithm A1 from Xu et al. (2025) supplementary.
    
    For each AFS, collect all vehicles that visit it with their arrival times.
    If more than eta_s vehicles overlap, delay some vehicles.
    
    Returns:
        (feasible, waiting_times): 
        - feasible: True if scheduling is possible within Tmax
        - waiting_times: dict mapping (route_idx, position) -> waiting_time
    """
    waiting_times = {}
    
    # Collect all AFS visits across all routes, grouped by station
    # First evaluate all routes without waiting
    route_evals = []
    for route in solution.routes:
        route_evals.append(evaluate_route(route, instance))
    
    # Group AFS visits by station
    station_visits = {}  # station_index -> list of (route_idx, pos_in_route, arrival_time)
    
    for r_idx, route in enumerate(solution.routes):
        eval_data = route_evals[r_idx]
        for pos, node in enumerate(route.nodes):
            if instance.is_station(node):
                # Map to actual station
                station_idx = node
                if station_idx not in station_visits:
                    station_visits[station_idx] = []
                station_visits[station_idx].append(
                    (r_idx, pos, eval_data['arrival_times'][pos])
                )
    
    # For each station, check capacity and resolve conflicts
    feasible = True
    
    for station_idx, visits in station_visits.items():
        capacity = instance.get_station_capacity(station_idx)
        refuel_time = instance.refueling_time
        
        if len(visits) <= capacity:
            continue  # No conflict
        
        # Sort by arrival time
        visits_sorted = sorted(visits, key=lambda x: x[2])
        
        # Apply reschedule algorithm
        arrival_times = [v[2] for v in visits_sorted]
        n_vehicles = len(arrival_times)
        
        # Reschedule: ensure at most capacity vehicles refuel simultaneously
        adjusted_arrivals = list(arrival_times)
        changed = True
        max_iterations = 1000
        iteration = 0
        
        while changed and iteration < max_iterations:
            changed = False
            iteration += 1
            for i in range(n_vehicles - 1):
                num_overlap = 0
                start_i = adjusted_arrivals[i]
                end_i = start_i + refuel_time
                
                for j in range(i + 1, n_vehicles):
                    start_j = adjusted_arrivals[j]
                    end_j = start_j + refuel_time
                    
                    if start_i < end_j and start_j < end_i:
                        num_overlap += 1
                        
                        if num_overlap > capacity - 1:
                            # Conflict: delay one vehicle
                            if end_j <= end_i:
                                adjusted_arrivals[i] = end_j
                            else:
                                adjusted_arrivals[j] = end_i
                            changed = True
                            break
                if changed:
                    break
        
        # Compute waiting times
        for k, (r_idx, pos, orig_arrival) in enumerate(visits_sorted):
            wt = adjusted_arrivals[k] - orig_arrival
            if wt > 1e-6:
                waiting_times[(r_idx, pos)] = wt
    
    # Check if waiting times cause any route to exceed Tmax
    for r_idx, route in enumerate(solution.routes):
        total_wait = sum(
            wt for (ri, pos), wt in waiting_times.items() if ri == r_idx
        )
        eval_data = route_evals[r_idx]
        if eval_data['duration'] + total_wait > instance.t_max + 1e-6:
            feasible = False
    
    return feasible, waiting_times


def evaluate_solution(solution: 'Solution') -> dict:
    """
    Full evaluation of a solution including AFS capacity constraints.
    
    Returns dict with:
    - total_distance: sum of all route distances
    - route_evaluations: list of per-route evaluations
    - all_customers_served: True if all customers are served exactly once
    - vehicle_count_feasible: True if n_routes <= n_vehicles
    - fuel_feasible: True if all routes satisfy fuel constraints
    - duration_feasible: True if all routes satisfy duration constraints
    - capacity_feasible: True if AFS capacity constraints can be satisfied
    - feasible: True if ALL constraints are satisfied
    - waiting_times: dict of waiting times from Reschedule
    """
    instance = solution.instance
    
    # Evaluate each route
    route_evals = []
    total_distance = 0.0
    fuel_feasible = True
    duration_feasible = True
    
    for route in solution.routes:
        eval_data = evaluate_route(route, instance)
        route_evals.append(eval_data)
        total_distance += eval_data['distance']
        if not eval_data['fuel_feasible']:
            fuel_feasible = False
        if not eval_data['duration_feasible']:
            duration_feasible = False
    
    # Check all customers served exactly once
    customer_count = {}
    for route in solution.routes:
        for node in route.nodes:
            if instance.is_customer(node):
                customer_count[node] = customer_count.get(node, 0) + 1
    
    all_served = (set(customer_count.keys()) == set(instance.customer_indices))
    no_duplicates = all(v == 1 for v in customer_count.values())
    customers_ok = all_served and no_duplicates
    
    # Check vehicle count
    vehicle_ok = len(solution.routes) <= instance.n_vehicles
    
    # Check AFS capacity via Reschedule
    capacity_feasible, waiting_times = reschedule(instance, solution)
    
    # Recheck duration with waiting times
    if waiting_times:
        for r_idx, route in enumerate(solution.routes):
            total_wait = sum(
                wt for (ri, pos), wt in waiting_times.items() if ri == r_idx
            )
            if route_evals[r_idx]['duration'] + total_wait > instance.t_max + 1e-6:
                duration_feasible = False
    
    feasible = (customers_ok and vehicle_ok and fuel_feasible 
                and duration_feasible and capacity_feasible)
    
    return {
        'total_distance': total_distance,
        'route_evaluations': route_evals,
        'all_customers_served': customers_ok,
        'vehicle_count_feasible': vehicle_ok,
        'fuel_feasible': fuel_feasible,
        'duration_feasible': duration_feasible,
        'capacity_feasible': capacity_feasible,
        'feasible': feasible,
        'waiting_times': waiting_times,
        'n_routes': len(solution.routes),
    }


# Attach evaluate method to Solution class
def _solution_evaluate(self) -> dict:
    if self._evaluation is None:
        self._evaluation = evaluate_solution(self)
    return self._evaluation

Solution.evaluate = _solution_evaluate


def compute_total_distance(solution: Solution) -> float:
    """Quick computation of total distance without full evaluation."""
    td = 0.0
    for route in solution.routes:
        for i in range(len(route.nodes) - 1):
            td += solution.instance.dist(route.nodes[i], route.nodes[i + 1])
    return td


def compute_penalty(solution: Solution, w_duration: float = 1.0, 
                     w_fuel: float = 1.0, w_capacity: float = 1.0) -> float:
    """
    Compute penalty for constraint violations.
    Used in metaheuristic to guide search through infeasible regions.
    
    Returns total penalty value (0 if feasible).
    """
    instance = solution.instance
    penalty = 0.0
    
    for route in solution.routes:
        eval_data = evaluate_route(route, instance)
        
        # Duration violation
        if eval_data['duration'] > instance.t_max:
            penalty += w_duration * (eval_data['duration'] - instance.t_max)
        
        # Fuel violation: check each segment
        for fl in eval_data['fuel_levels']:
            if fl < -1e-6:
                penalty += w_fuel * abs(fl)
    
    # Capacity violation
    capacity_feasible, waiting_times = reschedule(instance, solution)
    if not capacity_feasible:
        # Sum of excess waiting
        total_excess_wait = sum(waiting_times.values()) if waiting_times else 0
        penalty += w_capacity * (total_excess_wait + 1.0)
    
    return penalty
