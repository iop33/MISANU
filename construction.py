"""
Construction heuristics for GrVRP-PCAFS.

Implements:
1. Greedy nearest-neighbor construction
2. Savings-based construction (Clarke & Wright)
"""

import math
import numpy as np
from typing import List, Tuple, Optional, Set
from instance import Instance
from solution import Route, Solution, evaluate_route


def find_nearest_feasible_station(instance: Instance, from_node: int, 
                                   fuel_remaining: float) -> Optional[int]:
    """Find nearest AFS reachable with current fuel from given node."""
    best_station = None
    best_dist = float('inf')
    
    for s in instance.station_indices:
        d = instance.dist(from_node, s)
        fuel_needed = instance.consumption_rate * d
        if fuel_needed <= fuel_remaining + 1e-6 and d < best_dist:
            best_dist = d
            best_station = s
    
    return best_station


def can_return_to_depot(instance: Instance, from_node: int, 
                        fuel_remaining: float) -> bool:
    """Check if vehicle can return to depot from current node."""
    d = instance.dist(from_node, instance.depot)
    fuel_needed = instance.consumption_rate * d
    return fuel_needed <= fuel_remaining + 1e-6


def greedy_construction(instance: Instance) -> Solution:
    """
    Greedy nearest-neighbor construction heuristic.
    
    Build routes one by one:
    1. Start from depot with full tank
    2. Go to nearest unvisited customer that is reachable
    3. If fuel is low, insert AFS visit
    4. If no customer can be added, close route and start new one
    """
    solution = Solution(instance)
    unserved = set(instance.customer_indices)
    
    while unserved and len(solution.routes) < instance.n_vehicles:
        route_nodes = [0]  # start at depot
        current_node = 0
        fuel = instance.tank_capacity
        current_time = instance.p_start
        
        while unserved:
            # Find nearest unserved customer
            best_customer = None
            best_dist = float('inf')
            
            for c in unserved:
                d = instance.dist(current_node, c)
                fuel_needed = instance.consumption_rate * d
                time_after = current_time + instance.travel_time(current_node, c) + instance.service_time_customer
                
                # Check if we can reach customer and then return to depot
                # (possibly via AFS)
                d_to_depot = instance.dist(c, 0)
                fuel_after = fuel - fuel_needed
                fuel_to_depot = instance.consumption_rate * d_to_depot
                time_to_depot = time_after + instance.travel_time(c, 0)
                
                can_reach = fuel_needed <= fuel + 1e-6
                
                # Check if we can return (directly or via AFS)
                can_close = False
                if fuel_after >= fuel_to_depot - 1e-6 and time_to_depot <= instance.t_max + 1e-6:
                    can_close = True
                else:
                    # Check via nearest AFS
                    for s in instance.station_indices:
                        d_c_s = instance.dist(c, s)
                        d_s_0 = instance.dist(s, 0)
                        fuel_c_s = instance.consumption_rate * d_c_s
                        if fuel_after >= fuel_c_s - 1e-6:
                            time_via_s = (time_after + instance.travel_time(c, s) 
                                         + instance.refueling_time 
                                         + instance.travel_time(s, 0))
                            if time_via_s <= instance.t_max + 1e-6:
                                # After refueling, enough fuel to reach depot?
                                fuel_after_refuel = instance.tank_capacity
                                if fuel_after_refuel >= instance.consumption_rate * d_s_0 - 1e-6:
                                    can_close = True
                                    break
                
                if can_reach and can_close and d < best_dist:
                    best_dist = d
                    best_customer = c
            
            if best_customer is None:
                break  # No more customers can be added
            
            # Check if we need to refuel before visiting customer
            d_to_cust = instance.dist(current_node, best_customer)
            fuel_needed = instance.consumption_rate * d_to_cust
            
            if fuel_needed > fuel + 1e-6:
                # Need refueling first - try all stations, pick best
                best_station = None
                best_station_cost = float('inf')
                for s in instance.station_indices:
                    d_cs = instance.dist(current_node, s)
                    d_sc = instance.dist(s, best_customer)
                    if instance.consumption_rate * d_cs <= fuel + 1e-6:
                        if instance.consumption_rate * d_sc <= instance.tank_capacity + 1e-6:
                            cost = d_cs + d_sc
                            if cost < best_station_cost:
                                best_station_cost = cost
                                best_station = s
                
                if best_station is None:
                    break
                
                route_nodes.append(best_station)
                current_time += instance.travel_time(current_node, best_station) + instance.refueling_time
                fuel = instance.tank_capacity
                current_node = best_station
                
                # Recheck if customer is reachable
                d_to_cust = instance.dist(current_node, best_customer)
                fuel_needed = instance.consumption_rate * d_to_cust
                if fuel_needed > fuel + 1e-6:
                    break
            
            # Visit customer
            route_nodes.append(best_customer)
            fuel -= instance.consumption_rate * d_to_cust
            current_time += instance.travel_time(current_node, best_customer) + instance.service_time_customer
            current_node = best_customer
            unserved.remove(best_customer)
            
            # Check if we need refueling to return to depot
            d_to_depot = instance.dist(current_node, 0)
            fuel_to_depot = instance.consumption_rate * d_to_depot
            
            if fuel < fuel_to_depot - 1e-6:
                # Need refueling before returning
                station = find_nearest_feasible_station(instance, current_node, fuel)
                if station is not None:
                    d_to_s = instance.dist(current_node, station)
                    d_s_0 = instance.dist(station, 0)
                    time_check = current_time + instance.travel_time(current_node, station) + instance.refueling_time + instance.travel_time(station, 0)
                    if time_check <= instance.t_max + 1e-6:
                        route_nodes.append(station)
                        fuel = instance.tank_capacity
                        current_time += instance.travel_time(current_node, station) + instance.refueling_time
                        current_node = station
        
        # Close route
        route_nodes.append(0)
        route = Route(route_nodes)
        
        if not route.is_empty and len(route.customers(instance)) > 0:
            solution.add_route(route)
        else:
            # If route has no customers, try to force one
            break
    
    solution.remove_empty_routes()
    
    # Fallback: create individual routes for remaining unserved customers
    while unserved and len(solution.routes) < instance.n_vehicles:
        c = unserved.pop()
        d_to = instance.dist(0, c)
        d_from = instance.dist(c, 0)
        if instance.consumption_rate * (d_to + d_from) <= instance.tank_capacity + 1e-6:
            t = instance.p_start + instance.travel_time(0, c) + instance.service_time_customer + instance.travel_time(c, 0)
            if t <= instance.t_max + 1e-6:
                solution.add_route(Route([0, c, 0]))
                continue
        placed = False
        for s in instance.station_indices:
            d_0s = instance.dist(0, s)
            d_sc = instance.dist(s, c)
            d_cs = instance.dist(c, s)
            d_s0 = instance.dist(s, 0)
            if (instance.consumption_rate * d_0s <= instance.tank_capacity + 1e-6 and
                instance.consumption_rate * (d_sc + d_cs) <= instance.tank_capacity + 1e-6 and
                instance.consumption_rate * d_s0 <= instance.tank_capacity + 1e-6):
                t = (instance.p_start + instance.travel_time(0, s) + instance.refueling_time 
                     + instance.travel_time(s, c) + instance.service_time_customer 
                     + instance.travel_time(c, s) + instance.refueling_time 
                     + instance.travel_time(s, 0))
                if t <= instance.t_max + 1e-6:
                    solution.add_route(Route([0, s, c, s, 0]))
                    placed = True
                    break
    
    solution.remove_empty_routes()
    return solution


def savings_construction(instance: Instance) -> Solution:
    """
    Clarke & Wright savings-based construction.
    
    1. Start with one route per customer: depot -> customer -> depot
    2. Compute savings s(i,j) = d(i,0) + d(0,j) - d(i,j)
    3. Merge routes based on highest savings
    """
    # Initialize: one route per customer
    routes = {}
    for c in instance.customer_indices:
        routes[c] = Route([0, c, 0])
    
    # Compute savings
    savings = []
    for i in instance.customer_indices:
        for j in instance.customer_indices:
            if i != j:
                s = (instance.dist(i, 0) + instance.dist(0, j) 
                     - instance.dist(i, j))
                savings.append((s, i, j))
    
    savings.sort(reverse=True, key=lambda x: x[0])
    
    # Route assignment: which route is each customer in?
    customer_route = {c: c for c in instance.customer_indices}
    
    for saving, i, j in savings:
        if saving <= 0:
            break
        
        route_i_key = customer_route[i]
        route_j_key = customer_route[j]
        
        if route_i_key == route_j_key:
            continue  # Same route
        
        route_i = routes[route_i_key]
        route_j = routes[route_j_key]
        
        # Check if i is last customer in its route and j is first in its route
        customers_i = route_i.customers(instance)
        customers_j = route_j.customers(instance)
        
        if not customers_i or not customers_j:
            continue
            
        if customers_i[-1] != i or customers_j[0] != j:
            continue
        
        # Try to merge: route_i + route_j (remove depot in between)
        merged_nodes = route_i.nodes[:-1] + route_j.nodes[1:]
        merged_route = Route(merged_nodes)
        
        # Check feasibility
        eval_data = evaluate_route(merged_route, instance)
        
        if (eval_data['fuel_feasible'] and eval_data['duration_feasible'] 
            and len(merged_route.customers(instance)) <= instance.n_customers):
            
            # Check merged customers count for vehicle limit
            total_routes_after = len(routes) - 1
            
            # Perform merge
            routes[route_i_key] = merged_route
            del routes[route_j_key]
            
            # Update assignments
            for c in customers_j:
                customer_route[c] = route_i_key
    
    solution = Solution(instance, list(routes.values()))
    solution.remove_empty_routes()
    return solution


def insert_station_if_needed(route: Route, instance: Instance) -> Route:
    """
    Post-process a route: insert AFS visits where fuel runs out.
    Uses nearest feasible AFS insertion.
    """
    nodes = route.nodes
    new_nodes = [nodes[0]]
    fuel = instance.tank_capacity
    
    for i in range(1, len(nodes)):
        prev = new_nodes[-1]
        curr = nodes[i]
        d = instance.dist(prev, curr)
        fuel_needed = instance.consumption_rate * d
        
        if fuel_needed > fuel + 1e-6:
            # Need to insert AFS before this node
            station = find_nearest_feasible_station(instance, prev, fuel)
            if station is not None:
                new_nodes.append(station)
                fuel = instance.tank_capacity
                # Recompute fuel needed from station to curr
                fuel_needed = instance.consumption_rate * instance.dist(station, curr)
        
        new_nodes.append(curr)
        
        if instance.is_station(curr):
            fuel = instance.tank_capacity
        else:
            fuel -= fuel_needed
    
    return Route(new_nodes)
