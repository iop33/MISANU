"""
Instance data model for the Green Vehicle Routing Problem 
with Capacitated Alternative Fuel Stations (GrVRP-PCAFS / GVRP-CAFS).

Based on:
- Bruglieri, Mancini & Pisacane (2019), Computers & Operations Research
- Xu, Fan, Liu, Chen & Tang (2025), arXiv:2504.04527
"""

import math
import json
import os
import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional


@dataclass
class Instance:
    """Represents a GrVRP-PCAFS instance."""
    name: str
    
    # Node data: index 0 = depot, 1..n = customers, n+1..n+s = AFSs
    n_customers: int          # number of customers
    n_stations: int           # number of AFSs
    coords: np.ndarray        # (n_total, 2) coordinates of all nodes
    
    # Vehicle parameters
    n_vehicles: int           # max number of available AFVs
    speed: float              # average speed (miles/h or km/h)
    tank_capacity: float      # Q: max fuel/energy capacity
    consumption_rate: float   # r: fuel consumption per unit distance
    
    # Time parameters
    t_max: float              # max route duration
    service_time_customer: float  # service time at each customer
    refueling_time: float     # refueling time at each AFS (full refuel)
    p_start: float            # initial refuel time at depot
    
    # Station capacity
    station_capacity: np.ndarray  # eta_s: number of pumps at each AFS
    
    # Derived
    d_max: float = 0.0        # max driving range = Q / r
    
    # Distance matrix (computed)
    dist_matrix: np.ndarray = field(default=None, repr=False)
    time_matrix: np.ndarray = field(default=None, repr=False)
    
    def __post_init__(self):
        self.d_max = self.tank_capacity / self.consumption_rate
        self._compute_matrices()
    
    @property
    def n_total(self) -> int:
        """Total number of nodes (depot + customers + stations)."""
        return 1 + self.n_customers + self.n_stations
    
    @property
    def depot(self) -> int:
        return 0
    
    @property
    def customer_indices(self) -> List[int]:
        return list(range(1, self.n_customers + 1))
    
    @property
    def station_indices(self) -> List[int]:
        return list(range(self.n_customers + 1, self.n_total))
    
    def is_customer(self, node: int) -> bool:
        return 1 <= node <= self.n_customers
    
    def is_station(self, node: int) -> bool:
        return self.n_customers + 1 <= node < self.n_total
    
    def is_depot(self, node: int) -> bool:
        return node == 0
    
    def get_station_capacity(self, station_node: int) -> int:
        """Get capacity (number of pumps) for a station node."""
        idx = station_node - self.n_customers - 1
        return int(self.station_capacity[idx])
    
    def get_service_time(self, node: int) -> float:
        """Get service time at a node."""
        if self.is_customer(node):
            return self.service_time_customer
        return 0.0
    
    def get_refueling_time(self, node: int) -> float:
        """Get refueling time at a station (full refuel)."""
        if self.is_station(node):
            return self.refueling_time
        return 0.0
    
    def _compute_matrices(self):
        """Compute distance and time matrices."""
        n = self.n_total
        self.dist_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i != j:
                    dx = self.coords[i, 0] - self.coords[j, 0]
                    dy = self.coords[i, 1] - self.coords[j, 1]
                    self.dist_matrix[i, j] = math.sqrt(dx * dx + dy * dy)
        self.time_matrix = self.dist_matrix / self.speed
    
    def dist(self, i: int, j: int) -> float:
        return self.dist_matrix[i, j]
    
    def travel_time(self, i: int, j: int) -> float:
        return self.time_matrix[i, j]
    
    def save(self, filepath: str):
        """Save instance to JSON file."""
        data = {
            'name': self.name,
            'n_customers': self.n_customers,
            'n_stations': self.n_stations,
            'coords': self.coords.tolist(),
            'n_vehicles': self.n_vehicles,
            'speed': self.speed,
            'tank_capacity': self.tank_capacity,
            'consumption_rate': self.consumption_rate,
            't_max': self.t_max,
            'service_time_customer': self.service_time_customer,
            'refueling_time': self.refueling_time,
            'p_start': self.p_start,
            'station_capacity': self.station_capacity.tolist(),
        }
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2)
    
    @classmethod
    def load(cls, filepath: str) -> 'Instance':
        """Load instance from JSON file."""
        with open(filepath, 'r') as f:
            data = json.load(f)
        return cls(
            name=data['name'],
            n_customers=data['n_customers'],
            n_stations=data['n_stations'],
            coords=np.array(data['coords']),
            n_vehicles=data['n_vehicles'],
            speed=data['speed'],
            tank_capacity=data['tank_capacity'],
            consumption_rate=data['consumption_rate'],
            t_max=data['t_max'],
            service_time_customer=data['service_time_customer'],
            refueling_time=data['refueling_time'],
            p_start=data['p_start'],
            station_capacity=np.array(data['station_capacity']),
        )
    
    def __str__(self):
        return (f"Instance(name={self.name}, customers={self.n_customers}, "
                f"stations={self.n_stations}, vehicles={self.n_vehicles}, "
                f"Dmax={self.d_max:.1f}, Tmax={self.t_max})")
