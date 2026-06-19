"""
Generator benchmark instanci za GrVRP-PCAFS.

Pravi instance koje odgovaraju benchmark setovima opisanim u:
- Bruglieri et al. (2019): CENTRAL set (15 musterija)
- Bruglieri et al. (2022): CENTRAL set sa GRASP-om
- Xu et al. (2025): S-Central, M-Central, Beijing setovi

Rasporedi instanci:
- S-Central: 15 musterija, 1 punionica u centru, depo daleko
- M-Central: 25/50/100 musterija, 1 punionica, slican raspored
- EMH: zasnovano na Erdogan & Miller-Hooks instancama
- TRIANGLE: punionice izmedju depoa i zone musterija

Sve instance podrazumevano koriste jedinicni kapacitet punionice (eta_s = 1).
"""

# =============================================================================
# generate_instances.py  -- PRAVLJENJE TEST-PRIMERA (instanci)
# -----------------------------------------------------------------------------
# Nemamo pravi skup podataka, pa instance GENERISEMO proceduralno:
# postavimo depo, musterije i punionice na koordinate i zadamo parametre.
# Slucajnost je "zakucana" seed-om -> isti seed = iste koordinate (ponovljivo).
# Svaka funkcija pravi jedan TIP rasporeda (S-Central, M-Central, Triangle, EMH).
# =============================================================================

import os
import numpy as np                                      # nasumicne koordinate
import random                                           # (uvezeno)
from instance import Instance


def generate_s_central(instance_id: int, seed: int = None) -> Instance:
    """
    Napravi S-Central instancu (15 musterija, 1 punionica).

    Zasnovano na Bruglieri et al. (2019) CENTRAL setu i Xu et al. (2025) S-Central.

    Raspored:
    - Musterije nasumicno u zoni 50x50 oko centra (50, 50)
    - 1 punionica u centru zone musterija (50, 50)
    - Depo na (50, -30) (oko 2 sata na 40 mph = 80 milja)
    """
    # ===== S-Central: MALI primer -- 15 musterija, 1 punionica u centru. =====
    if seed is not None:
        np.random.seed(seed)                            # zakucaj slucajnost (npr. seed=1001)

    n_customers = 15
    n_stations = 1
    n_total = 1 + n_customers + n_stations              # = 17 cvorova ukupno

    coords = np.zeros((n_total, 2))                     # tabela (x,y) za svaki cvor

    # Depo dole
    coords[0] = [50.0, -30.0]                           # DEPO -- namerno daleko dole

    # Musterije u zoni 50x50 oko centra (50, 50)
    for i in range(1, n_customers + 1):                 # 15 MUSTERIJA -- nasumicno u kvadratu [25,75]
        coords[i] = [
            np.random.uniform(25, 75),
            np.random.uniform(25, 75)
        ]

    # Punionica u centru zone musterija
    coords[n_customers + 1] = [50.0, 50.0]              # PUNIONICA -- u centru musterija

    # Kapacitet punionice
    station_capacity = np.array([1])                    # eta = 1 (jedna pumpa)

    return Instance(                                    # spakuj sve u objekat Instance
        name=f"S-Central_{instance_id}",
        n_customers=n_customers,
        n_stations=n_stations,
        coords=coords,
        n_vehicles=15,
        speed=40.0,           # 40 milja/h
        tank_capacity=50.0,   # Q = 50
        consumption_rate=0.2, # r = 0.2 -> domet = 50/0.2 = 250
        t_max=7.0,            # 7 sati
        service_time_customer=0.5,  # 0.5 sati
        refueling_time=0.5,         # 0.5 sati
        p_start=0.0,
        station_capacity=station_capacity,
    )


def generate_m_central(instance_id: int, n_customers: int = 25,
                        seed: int = None) -> Instance:
    """
    Napravi M-Central instancu (25, 50 ili 100 musterija).

    Zasnovano na Xu et al. (2025) M-Central setu.

    Raspored slican S-Central, ali skaliran.
    Kapacitet punionice zavisi od velicine: 2 za 25, 3 za 50, 8 za 100.
    """
    # ===== M-Central: SREDNJI/VELIKI primer -- 25/50/100 musterija, 1 punionica. =====
    if seed is not None:
        np.random.seed(seed)

    assert n_customers in [25, 50, 100], "n_customers must be 25, 50, or 100"  # dozvoljene velicine

    n_stations = 1
    n_total = 1 + n_customers + n_stations

    # Zona se skalira sa brojem musterija
    area_size = 50 + n_customers                        # veca zona kad ima vise musterija
    center = area_size / 2 + 25

    coords = np.zeros((n_total, 2))

    # Depo daleko od zone musterija
    coords[0] = [center, center - area_size * 0.8]      # DEPO -- daleko ispod zone

    # Musterije nasumicno u zoni
    for i in range(1, n_customers + 1):                 # MUSTERIJE -- nasumicno u zoni
        coords[i] = [
            np.random.uniform(center - area_size/2, center + area_size/2),
            np.random.uniform(center - area_size/2, center + area_size/2)
        ]

    # Punionica u centru
    coords[n_customers + 1] = [center, center]          # PUNIONICA -- u centru

    # Kapacitet zavisi od velicine
    capacity_map = {25: 2, 50: 3, 100: 8}               # broj pumpi (eta) zavisi od velicine
    n_vehicles_map = {25: 7, 50: 13, 100: 25}           # broj vozila zavisi od velicine

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
        t_max=7.5,  # 7.5 sati za M-Central
        service_time_customer=0.5,
        refueling_time=0.5,
        p_start=0.0,
        station_capacity=station_capacity,
    )


def generate_triangle(instance_id: int, seed: int = None) -> Instance:
    """
    Napravi TRIANGLE instancu (15 musterija, 3 punionice).

    Zasnovano na Bruglieri et al. (2019) TRIANGLE setu.

    Raspored:
    - Punionice na sredini, izmedju depoa i zone musterija
    - Svako vozilo mora da se dopuni
    """
    # ===== Triangle: 3 punionice NA SREDINI izmedju depoa i musterija ->
    #       svako vozilo mora da se dopuni (tu se kapacitet stvarno testira). =====
    if seed is not None:
        np.random.seed(seed)

    n_customers = 15
    n_stations = 3
    n_total = 1 + n_customers + n_stations

    coords = np.zeros((n_total, 2))

    # Depo dole
    coords[0] = [50.0, 0.0]                             # DEPO -- dole

    # Musterije u gornjoj zoni (50x50 oko (50, 120))
    for i in range(1, n_customers + 1):                 # MUSTERIJE -- gore
        coords[i] = [
            np.random.uniform(25, 75),
            np.random.uniform(95, 145)
        ]

    # 3 punionice na sredini (izmedju depoa i musterija)
    coords[n_customers + 1] = [30.0, 60.0]              # 3 PUNIONICE -- na sredini puta
    coords[n_customers + 2] = [50.0, 60.0]
    coords[n_customers + 3] = [70.0, 60.0]

    station_capacity = np.array([1, 1, 1])              # eta = 1 po svakoj punionici

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
    Napravi EMH instancu (musterije i punionice izmesane u istom prostoru).

    Zasnovano na strukturi Erdogan & Miller-Hooks (2012) benchmarka.
    """
    # ===== EMH: sve IZMESANO u istom prostoru -- 20 musterija, 6 punionica. =====
    if seed is not None:
        np.random.seed(seed)

    n_total = 1 + n_customers + n_stations
    coords = np.zeros((n_total, 2))

    # Svi cvorovi u zoni 100x100
    coords[0] = [50.0, 50.0]                            # DEPO -- u centru

    for i in range(1, n_total):                         # SVI ostali (musterije + punionice) -- nasumicno
        coords[i] = [
            np.random.uniform(0, 100),
            np.random.uniform(0, 100)
        ]

    station_capacity = np.ones(n_stations, dtype=int)   # eta = 1 po svakoj punionici

    return Instance(
        name=f"EMH_{n_customers}c{n_stations}s_{instance_id}",
        n_customers=n_customers,
        n_stations=n_stations,
        coords=coords,
        n_vehicles=max(3, n_customers // 3),
        speed=40.0,
        tank_capacity=60.0,                             # veci rezervoar -> domet 300
        consumption_rate=0.2,
        t_max=11.0,
        service_time_customer=0.5,
        refueling_time=0.25,
        p_start=0.25,
        station_capacity=station_capacity,
    )


def generate_all_instances(output_dir: str = "instances",
                           n_per_set: int = 10):
    """Napravi sve benchmark setove instanci i snimi ih u fajlove."""
    # ===== Napravi SVE setove i SNIMI ih kao .json u folder "instances/".
    #       (Napomena: main.py ovo NE poziva -- on instance pravi u memoriji.) =====
    os.makedirs(output_dir, exist_ok=True)

    instance_sets = {}

    # S-Central set (10 instanci)
    print("Generating S-Central instances...")
    s_central = []
    for i in range(1, n_per_set + 1):
        inst = generate_s_central(i, seed=1000 + i)     # seed = 1000+i (deterministicno)
        inst.save(os.path.join(output_dir, f"{inst.name}.json"))
        s_central.append(inst)
    instance_sets['S-Central'] = s_central

    # M-Central setovi
    for n_cust in [25, 50, 100]:                        # tri velicine M-Central
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
    generate_all_instances()                            # ako se fajl pokrene direktno -> napravi i snimi sve instance
