"""
Instance generator for GrVRP-PCAFS benchmark instances.

Generates instances matching the benchmark sets described in:
- Bruglieri et al. (2019): CENTRAL set (15 customers)
- Bruglieri et al. (2022): CENTRAL set with GRASP
- Xu et al. (2025): S-Central, M-Central, Beijing sets

Instance layouts:
- S-Central: 15 customers, 1 AFS at center, depot far away
- M-Central: 25/50/100 customers, 1 AFS, similar layout
- EMH: Based on Erdogan & Miller-Hooks instances
- TRIANGLE: AFSs between depot and customer area

All instances use unit AFS capacity (eta_s = 1) by default.
"""

import os
import numpy as np
import random
from instance import Instance


def generate_s_central(instance_id: int, seed: int = None) -> Instance:
    """
    Generate S-Central instance (15 customers, 1 AFS).
    
    Based on Bruglieri et al. (2019) CENTRAL set and Xu et al. (2025) S-Central.
    
    Layout:
    - Customers randomly placed in a 50x50 area centered at (50, 50)
    - 1 AFS at center of customer area (50, 50)
    - Depot at (50, -30) (about 2 hours away at 40 mph = 80 miles)
    """
    if seed is not None:
        np.random.seed(seed)
    
    n_customers = 15
    n_stations = 1
    n_total = 1 + n_customers + n_stations
    
    coords = np.zeros((n_total, 2))
    
    # Depot at bottom
    coords[0] = [50.0, -30.0]
    
    # Customers in 50x50 area centered at (50, 50)
    for i in range(1, n_customers + 1):
        coords[i] = [
            np.random.uniform(25, 75),
            np.random.uniform(25, 75)
        ]
    
    # AFS at center of customer area
    coords[n_customers + 1] = [50.0, 50.0]
    
    # Station capacity
    station_capacity = np.array([1])
    
    return Instance(
        name=f"S-Central_{instance_id}",
        n_customers=n_customers,
        n_stations=n_stations,
        coords=coords,
        n_vehicles=15,
        speed=40.0,           # 40 mph
        tank_capacity=50.0,   # Q = 50
        consumption_rate=0.2, # r = 0.2 gallons/mile → Dmax = 250
        t_max=7.0,            # 7 hours
        service_time_customer=0.5,  # 0.5 hours
        refueling_time=0.5,         # 0.5 hours
        p_start=0.0,
        station_capacity=station_capacity,
    )


def generate_m_central(instance_id: int, n_customers: int = 25, 
                        seed: int = None) -> Instance:
    """
    Generate M-Central instance (25, 50, or 100 customers).
    
    Based on Xu et al. (2025) M-Central set.
    
    Layout similar to S-Central but scaled.
    AFS capacity varies with size: 2 for 25, 3 for 50, 8 for 100.
    """
    if seed is not None:
        np.random.seed(seed)
    
    assert n_customers in [25, 50, 100], "n_customers must be 25, 50, or 100"
    
    n_stations = 1
    n_total = 1 + n_customers + n_stations
    
    # Scale area with number of customers
    area_size = 50 + n_customers
    center = area_size / 2 + 25
    
    coords = np.zeros((n_total, 2))
    
    # Depot far from customer area
    coords[0] = [center, center - area_size * 0.8]
    
    # Customers randomly in area
    for i in range(1, n_customers + 1):
        coords[i] = [
            np.random.uniform(center - area_size/2, center + area_size/2),
            np.random.uniform(center - area_size/2, center + area_size/2)
        ]
    
    # AFS at center
    coords[n_customers + 1] = [center, center]
    
    # Capacity depends on size
    capacity_map = {25: 2, 50: 3, 100: 8}
    n_vehicles_map = {25: 7, 50: 13, 100: 25}
    
    station_capacity = np.array([capacity_map[n_customers]])
    
    return Instance(
        name=f"M-Central{n_customers}_{instance_id}",
        n_customers=n_customers,
        n_stations=n_stations,
        coords=coords,
        n_vehicles=n_vehicles_map[n_customers],
        speed=40.0,
        tank_capacity=50.0,
        consumption_rate=0.2,
        t_max=7.5,  # 7.5 hours for M-Central
        service_time_customer=0.5,
        refueling_time=0.5,
        p_start=0.0,
        station_capacity=station_capacity,
    )


def generate_triangle(instance_id: int, seed: int = None) -> Instance:
    """
    Generate TRIANGLE instance (15 customers, 3 AFSs).
    
    Based on Bruglieri et al. (2019) TRIANGLE set.
    
    Layout:
    - AFSs in the middle between depot and customer area
    - Every vehicle needs to refuel
    """
    if seed is not None:
        np.random.seed(seed)
    
    n_customers = 15
    n_stations = 3
    n_total = 1 + n_customers + n_stations
    
    coords = np.zeros((n_total, 2))
    
    # Depot at bottom
    coords[0] = [50.0, 0.0]
    
    # Customers in upper area (50x50 centered at 50, 120)
    for i in range(1, n_customers + 1):
        coords[i] = [
            np.random.uniform(25, 75),
            np.random.uniform(95, 145)
        ]
    
    # 3 AFSs in the middle (between depot and customers)
    coords[n_customers + 1] = [30.0, 60.0]
    coords[n_customers + 2] = [50.0, 60.0]
    coords[n_customers + 3] = [70.0, 60.0]
    
    station_capacity = np.array([1, 1, 1])
    
    return Instance(
        name=f"Triangle_{instance_id}",
        n_customers=n_customers,
        n_stations=n_stations,
        coords=coords,
        n_vehicles=10,
        speed=40.0,
        tank_capacity=50.0,
        consumption_rate=0.2,
        t_max=11.0,
        service_time_customer=0.75,
        refueling_time=0.5,
        p_start=0.0,
        station_capacity=station_capacity,
    )


def generate_emh_like(instance_id: int, n_customers: int = 20, 
                       n_stations: int = 6, seed: int = None) -> Instance:
    """
    Generate EMH-like instance (customers and AFSs integrated in same area).
    
    Based on Erdogan & Miller-Hooks (2012) benchmark structure.
    """
    if seed is not None:
        np.random.seed(seed)
    
    n_total = 1 + n_customers + n_stations
    coords = np.zeros((n_total, 2))
    
    # All nodes in a 100x100 area
    coords[0] = [50.0, 50.0]  # Depot at center
    
    for i in range(1, n_total):
        coords[i] = [
            np.random.uniform(0, 100),
            np.random.uniform(0, 100)
        ]
    
    station_capacity = np.ones(n_stations, dtype=int)
    
    return Instance(
        name=f"EMH_{n_customers}c{n_stations}s_{instance_id}",
        n_customers=n_customers,
        n_stations=n_stations,
        coords=coords,
        n_vehicles=max(3, n_customers // 3),
        speed=40.0,
        tank_capacity=60.0,
        consumption_rate=0.2,
        t_max=11.0,
        service_time_customer=0.5,
        refueling_time=0.25,
        p_start=0.25,
        station_capacity=station_capacity,
    )


def generate_all_instances(output_dir: str = "instances", 
                           n_per_set: int = 10):
    """Generate all benchmark instance sets and save to files."""
    os.makedirs(output_dir, exist_ok=True)
    
    instance_sets = {}
    
    # S-Central set (10 instances)
    print("Generating S-Central instances...")
    s_central = []
    for i in range(1, n_per_set + 1):
        inst = generate_s_central(i, seed=1000 + i)
        inst.save(os.path.join(output_dir, f"{inst.name}.json"))
        s_central.append(inst)
    instance_sets['S-Central'] = s_central
    
    # M-Central sets
    for n_cust in [25, 50, 100]:
        print(f"Generating M-Central{n_cust} instances...")
        m_central = []
        for i in range(1, n_per_set + 1):
            inst = generate_m_central(i, n_customers=n_cust, seed=2000 + n_cust * 100 + i)
            inst.save(os.path.join(output_dir, f"{inst.name}.json"))
            m_central.append(inst)
        instance_sets[f'M-Central{n_cust}'] = m_central
    
    # Triangle set
    print("Generating Triangle instances...")
    triangles = []
    for i in range(1, n_per_set + 1):
        inst = generate_triangle(i, seed=3000 + i)
        inst.save(os.path.join(output_dir, f"{inst.name}.json"))
        triangles.append(inst)
    instance_sets['Triangle'] = triangles
    
    # EMH set
    print("Generating EMH instances...")
    emh = []
    for i in range(1, n_per_set + 1):
        inst = generate_emh_like(i, seed=4000 + i)
        inst.save(os.path.join(output_dir, f"{inst.name}.json"))
        emh.append(inst)
    instance_sets['EMH'] = emh
    
    print(f"\nGenerated {sum(len(v) for v in instance_sets.values())} instances "
          f"in {len(instance_sets)} sets")
    print(f"Saved to: {output_dir}/")
    
    return instance_sets


if __name__ == "__main__":
    generate_all_instances()
