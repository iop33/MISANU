"""
Model podataka jedne instance za zeleni problem rutiranja vozila
sa punionicama ogranicenog kapaciteta (GrVRP-PCAFS / GVRP-CAFS).

Zasnovano na:
- Bruglieri, Mancini & Pisacane (2019), Computers & Operations Research
- Xu, Fan, Liu, Chen & Tang (2025), arXiv:2504.04527
"""

# =============================================================================
# instance.py  -- MODEL PODATAKA JEDNE INSTANCE (tj. "tekst zadatka")
# -----------------------------------------------------------------------------
# Ovde se cuvaju svi STATICKI podaci o jednom problemu:
#   - gde su depo, musterije i punionice (koordinate),
#   - parametri vozila (rezervoar, potrosnja, brzina),
#   - vremenska ogranicenja (T_max, vreme usluge, vreme dopune),
#   - kapacitet punionica (broj pumpi eta).
# Ovaj fajl NE sadrzi algoritam ni racunanje ruta -- samo podatke i pomocne
# upite (rastojanje, vreme, tip cvora). Sve ostalo (resavanje) je u drugim fajlovima.
# =============================================================================

import math                                   # za sqrt (euklidsko rastojanje)
import json                                    # za upis/citanje instance u .json
import os                                      # (uvezeno, ali se ne koristi u ovom fajlu)
import numpy as np                             # matrice i niz koordinata
from dataclasses import dataclass, field       # @dataclass = automatski konstruktor
from typing import List, Tuple, Optional       # samo za oznake tipova (citljivost)


@dataclass                                     # @dataclass sam pravi __init__ od polja ispod
class Instance:
    """Predstavlja jednu GrVRP-PCAFS instancu."""
    name: str                                  # naziv instance, npr. "S-Central_1"

    # Podaci o cvorovima. KONVENCIJA INDEKSA (vazi kroz CEO projekat):
    #   cvor 0            -> depo
    #   cvorovi 1..n      -> musterije
    #   cvorovi n+1..n+s  -> punionice (AFS)
    n_customers: int          # koliko ima musterija
    n_stations: int           # koliko ima punionica
    coords: np.ndarray        # tabela (x,y) za svaki cvor, oblik (n_total, 2)

    # ----- parametri vozila -----
    n_vehicles: int           # najveci broj vozila u floti
    speed: float              # prosecna brzina (milje/h ili km/h)
    tank_capacity: float      # Q = velicina rezervoara (max gorivo/energija)
    consumption_rate: float   # r = potrosnja goriva po jedinici rastojanja

    # ----- vremenski parametri -----
    t_max: float              # T_max = najduze dozvoljeno trajanje rute
    service_time_customer: float  # vreme zadrzavanja (usluge) kod musterije
    refueling_time: float     # vreme dopune na punionici (do punog)
    p_start: float            # pocetno vreme pripreme/dopune u depou

    # ----- kapacitet punionica -----
    station_capacity: np.ndarray  # eta_s = broj pumpi po svakoj punionici

    # ----- izvedeno (racuna se automatski) -----
    d_max: float = 0.0        # domet na pun rezervoar = Q / r (popunjava se nize)

    # ----- matrice (racunaju se na startu) -----
    dist_matrix: np.ndarray = field(default=None, repr=False)  # matrica rastojanja izmedju svih cvorova
    time_matrix: np.ndarray = field(default=None, repr=False)  # matrica vremena putovanja

    def __post_init__(self):
        # __post_init__ se POKRECE AUTOMATSKI cim se napravi Instance(...).
        self.d_max = self.tank_capacity / self.consumption_rate   # domet = Q / r (npr. 50/0.2 = 250)
        self._compute_matrices()                                   # izracunaj matrice rastojanja i vremena

    @property
    def n_total(self) -> int:
        """Ukupan broj cvorova (depo + musterije + punionice)."""
        return 1 + self.n_customers + self.n_stations    # 1 depo + sve musterije + sve punionice

    @property
    def depot(self) -> int:
        return 0                                          # depo je uvek cvor sa indeksom 0

    @property
    def customer_indices(self) -> List[int]:
        return list(range(1, self.n_customers + 1))       # indeksi musterija: 1, 2, ..., n

    @property
    def station_indices(self) -> List[int]:
        return list(range(self.n_customers + 1, self.n_total))  # indeksi punionica: n+1, ..., kraj

    def is_customer(self, node: int) -> bool:
        return 1 <= node <= self.n_customers              # da li je dati cvor musterija?

    def is_station(self, node: int) -> bool:
        return self.n_customers + 1 <= node < self.n_total  # da li je dati cvor punionica?

    def is_depot(self, node: int) -> bool:
        return node == 0                                  # da li je dati cvor depo?

    def get_station_capacity(self, station_node: int) -> int:
        """Vrati kapacitet (broj pumpi) za cvor punionice."""
        # cvor punionice -> indeks u nizu station_capacity (oduzmemo depo i musterije)
        idx = station_node - self.n_customers - 1
        return int(self.station_capacity[idx])            # vrati broj pumpi (eta) te stanice

    def get_service_time(self, node: int) -> float:
        """Vrati vreme usluge u datom cvoru."""
        if self.is_customer(node):
            return self.service_time_customer             # musterija -> vreme usluge
        return 0.0                                        # depo i punionica -> nema vremena usluge

    def get_refueling_time(self, node: int) -> float:
        """Vrati vreme dopune na punionici (dopuna do punog)."""
        if self.is_station(node):
            return self.refueling_time                    # punionica -> vreme dopune
        return 0.0                                        # ostalo -> nema dopune

    def _compute_matrices(self):
        """Izracunaj matrice rastojanja i vremena."""
        n = self.n_total                                  # ukupan broj cvorova
        self.dist_matrix = np.zeros((n, n))               # prazna n x n matrica rastojanja
        for i in range(n):                                # za svaki par cvorova (i, j)...
            for j in range(n):
                if i != j:                                # rastojanje cvora do samog sebe ostaje 0
                    dx = self.coords[i, 0] - self.coords[j, 0]   # razlika po x
                    dy = self.coords[i, 1] - self.coords[j, 1]   # razlika po y
                    self.dist_matrix[i, j] = math.sqrt(dx * dx + dy * dy)  # EUKLIDSKO (vazdusno) rastojanje
        self.time_matrix = self.dist_matrix / self.speed  # vreme = rastojanje / brzina (cela matrica odjednom)

    def dist(self, i: int, j: int) -> float:
        return self.dist_matrix[i, j]                     # brzo CITANJE rastojanja (matrica vec izracunata)

    def travel_time(self, i: int, j: int) -> float:
        return self.time_matrix[i, j]                     # brzo citanje vremena putovanja

    def save(self, filepath: str):
        """Snimi instancu u JSON fajl."""
        # Pretvori instancu u obican recnik i upisi ga kao .json.
        # NAPOMENA: matrice i d_max se NE cuvaju -- ponovo se izracunaju pri ucitavanju.
        data = {
            'name': self.name,
            'n_customers': self.n_customers,
            'n_stations': self.n_stations,
            'coords': self.coords.tolist(),               # numpy niz -> obicna lista (da moze u json)
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
            json.dump(data, f, indent=2)                  # upis u fajl, lepo uvuceno (indent=2)

    @classmethod
    def load(cls, filepath: str) -> 'Instance':
        """Ucitaj instancu iz JSON fajla."""
        with open(filepath, 'r') as f:
            data = json.load(f)                           # procitaj recnik iz .json
        return cls(                                       # napravi Instance od procitanih podataka
            name=data['name'],
            n_customers=data['n_customers'],
            n_stations=data['n_stations'],
            coords=np.array(data['coords']),              # lista -> numpy niz
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
        # Citljiv tekstualni opis instance (za print()).
        return (f"Instance(name={self.name}, customers={self.n_customers}, "
                f"stations={self.n_stations}, vehicles={self.n_vehicles}, "
                f"Dmax={self.d_max:.1f}, Tmax={self.t_max})")
