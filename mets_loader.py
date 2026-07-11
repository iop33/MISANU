"""
Ucitavac PRAVIH benchmark instanci iz METS repozitorijuma (Xu et al. 2025).

Instanca je u MATLAB '.mat' formatu (struct 'vrp'). Ovaj modul je cita BEZ
MATLAB-a, preko scipy.io.loadmat, i pretvara u nas Instance model.

KLJUCNO: koristimo ORIGINALNU 'distance_table' iz fajla (NE preracunavamo
euklidsko rastojanje), da bi poredjenje sa METS-om bilo egzaktno na istim
brojevima.

Polja u 'vrp' struct-u (utvrdjeno dekodiranjem 15_1.mat):
    id              - oznake cvorova (D, C1..Cn, S1..Ss)
    type            - tip cvora ('d'=depo, 'c'=musterija, 'f'=punionica)
    longitude/latitude - koordinate (x, y)
    nb_customer     - broj musterija n
    distance_table  - puna matrica rastojanja (n_total x n_total)
    V_fuel          - Q (kapacitet rezervoara)
    V_fuel_rate     - r (potrosnja po jedinici rastojanja)
    V_Dmax          - domet na pun rezervoar (= Q/r)
    V_speed         - prosecna brzina
    V_nb            - broj vozila M
    T_max_V         - T_max (najduze trajanje rute)
    T_Start         - p_start (pocetna priprema u depou)
    T_Afs           - vreme dopune na punionici
    T_Customer      - vreme usluge kod musterije
    C_Afs           - kapacitet punionice (eta_s, broj pumpi)

Konvencija indeksa (poklapa se sa nasim Instance modelom):
    cvor 0       -> depo
    cvor 1..n    -> musterije
    cvor n+1..   -> punionice
"""

import os
import glob
import numpy as np
import scipy.io as sio

from instance import Instance


def _scalar(x):
    """Izvuci skalarnu vrednost iz ugnezdjenog numpy niza (npr. array([[15]]) -> 15.0)."""
    a = np.asarray(x).ravel()
    return a[0]


def _column(x):
    """Izvuci 1D niz iz kolone (npr. shape (n,1) -> shape (n,))."""
    return np.asarray(x).ravel()


def load_mat_instance(path: str, name: str = None) -> Instance:
    """
    Ucitaj jednu '.mat' instancu i vrati nas Instance objekat.

    Distance/time matrice se preuzimaju iz originalne 'distance_table'
    (ne preracunavaju se), a d_max se postavlja na V_Dmax iz fajla.
    """
    mat = sio.loadmat(path)
    vrp = mat['vrp'][0, 0]                       # raspakuj struct 'vrp'

    # ----- tipovi cvorova -> broj musterija i punionica -----
    types = [str(t[0]).strip().lower() for t in _column(vrp['type'])]
    n_customers = int(_scalar(vrp['nb_customer']))
    n_stations = sum(1 for t in types if t == 'f')
    n_total = len(types)
    assert n_total == 1 + n_customers + n_stations, (
        f"Neslaganje broja cvorova u {path}: total={n_total}, "
        f"n_cust={n_customers}, n_stat={n_stations}"
    )

    # ----- koordinate (x=longitude, y=latitude) -----
    lon = _column(vrp['longitude']).astype(float)
    lat = _column(vrp['latitude']).astype(float)
    coords = np.column_stack([lon, lat])

    # ----- parametri vozila / vremena -----
    Q = float(_scalar(vrp['V_fuel']))
    r = float(_scalar(vrp['V_fuel_rate']))
    Dmax = float(_scalar(vrp['V_Dmax']))
    speed = float(_scalar(vrp['V_speed']))
    M = int(_scalar(vrp['V_nb']))
    t_max = float(_scalar(vrp['T_max_V']))
    p_start = float(_scalar(vrp['T_Start']))
    t_afs = float(_scalar(vrp['T_Afs']))
    t_cust = float(_scalar(vrp['T_Customer']))
    eta = int(_scalar(vrp['C_Afs']))

    station_capacity = np.full(n_stations, eta, dtype=int)

    inst = Instance(
        name=name or os.path.splitext(os.path.basename(path))[0],
        n_customers=n_customers,
        n_stations=n_stations,
        coords=coords,
        n_vehicles=M,
        speed=speed,
        tank_capacity=Q,
        consumption_rate=r,
        t_max=t_max,
        service_time_customer=t_cust,
        refueling_time=t_afs,
        p_start=p_start,
        station_capacity=station_capacity,
    )

    # ----- KLJUCNO: zameni preracunate matrice ORIGINALNOM distance_table -----
    # METS interno KLONIRA jedinu punionicu, pa je distance_table prosirena na
    # (2*n_cust+2) cvorova. Nama treba samo gornji-levi blok [D, C1..Cn, S1..Ss].
    dist_table = np.asarray(vrp['distance_table'], dtype=float)
    assert dist_table.shape[0] >= n_total and dist_table.shape[1] >= n_total, (
        f"distance_table {dist_table.shape} manja od ocekivanog {n_total}"
    )
    dist_table = dist_table[:n_total, :n_total]  # uzmi samo blok stvarnih cvorova
    inst.dist_matrix = dist_table
    inst.time_matrix = dist_table / speed
    inst.d_max = Dmax                            # koristi V_Dmax direktno (= Q/r)

    return inst


# Mapiranje prefiksa fajla -> naziv seta (radi citljivih izvestaja)
_SET_OF_PREFIX = {
    '15': 'S-Central',
    '25': 'M-Central25',
    '50': 'M-Central50',
    '100': 'M-Central100',
}


def instance_set_name(filename: str) -> str:
    """Vrati naziv seta na osnovu imena fajla (npr. '50_3.mat' -> 'M-Central50')."""
    base = os.path.basename(filename)
    if base.startswith('jd'):
        return 'Beijing'
    prefix = base.split('_')[0]
    return _SET_OF_PREFIX.get(prefix, prefix)


def load_set(instances_dir: str, prefix: str):
    """
    Ucitaj sve instance jednog seta (npr. prefix='15' za S-Central),
    sortirane po rednom broju. Vraca listu Instance objekata.
    """
    paths = glob.glob(os.path.join(instances_dir, f"{prefix}_*.mat"))

    def _idx(p):
        base = os.path.splitext(os.path.basename(p))[0]
        return int(base.split('_')[1])

    paths = sorted(paths, key=_idx)
    return [load_mat_instance(p, name=f"{instance_set_name(p)}_{_idx(p)}") for p in paths]


if __name__ == "__main__":
    import sys
    d = sys.argv[1] if len(sys.argv) > 1 else "/tmp/METS-Algorithm/Instances"
    inst = load_mat_instance(os.path.join(d, "15_1.mat"))
    print(inst)
    print("  Dmax =", inst.d_max, "| Tmax =", inst.t_max,
          "| eta =", inst.station_capacity.tolist(),
          "| M =", inst.n_vehicles)
