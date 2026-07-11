"""
Predstavljanje i ocena resenja za GrVRP-PCAFS.

Resenje je skup ruta. Svaka ruta je lista indeksa cvorova koja
pocinje i zavrsava u depou (0). Rute mogu sadrzati i posete punionicama.

Ukljucuje Reschedule proceduru iz Bruglieri et al. (2022) / Xu et al. (2025)
za upravljanje ogranicenim kapacitetom punionica.
"""

# =============================================================================
# solution.py  -- RESENJE + OCENA RESENJA (srce projekta)
# -----------------------------------------------------------------------------
# Sadrzi:
#   - klasu Route    (jedna ruta = niz cvorova [0, ..., 0])
#   - klasu Solution (skup ruta koji pokriva sve musterije)
#   - evaluate_route   : oceni jednu rutu (kilometraza, vreme, gorivo)
#   - reschedule       : resi RED NA PUMPI (kapacitet punionica) -- iz Xu et al. 2025
#   - evaluate_solution: puna ocena celog resenja + provera svih ogranicenja
#   - compute_penalty  : KAZNA za prekrsaje (vodi pretragu kroz nedozvoljeni prostor)
# =============================================================================

import copy                                    # (uvezeno; duboko kopiranje)
import math                                     # (uvezeno)
import numpy as np                              # numericke operacije
from typing import List, Tuple, Optional, Dict  # oznake tipova
from instance import Instance                   # koristimo podatke instance (rastojanja, parametri)


class Route:
    """Jedna ruta vozila: niz cvorova koji pocinje i zavrsava u depou."""
    # RUTA = put jednog vozila. Npr. [0, 3, 5, 11, 2, 0]:
    #   0 = depo (pocetak i kraj), 3/5/2 = musterije, 11 = punionica.

    def __init__(self, nodes: List[int] = None):
        """
        nodes: lista indeksa cvorova, npr. [0, 3, 5, 11, 2, 0]
        gde je 0=depo, 3,5,2=musterije, 11=punionica
        """
        if nodes is None:
            self.nodes = [0, 0]                       # prazna ruta = samo depo->depo
        else:
            self.nodes = list(nodes)                  # kopija liste cvorova (da se ne deli referenca)

    @property
    def is_empty(self) -> bool:
        """Ruta je prazna ako sadrzi samo depo cvorove."""
        return len(self.nodes) <= 2                   # [0,0] -> prazna (nema nijedne musterije)

    def customers(self, instance: Instance) -> List[int]:
        """Vrati listu cvorova-musterija u ovoj ruti."""
        return [n for n in self.nodes if instance.is_customer(n)]  # samo cvorovi-musterije iz rute

    def stations(self, instance: Instance) -> List[int]:
        """Vrati listu cvorova-punionica u ovoj ruti."""
        return [n for n in self.nodes if instance.is_station(n)]   # samo cvorovi-punionice iz rute

    def insert(self, position: int, node: int):
        """Ubaci cvor na datu poziciju u ruti."""
        self.nodes.insert(position, node)             # ubaci cvor na datu poziciju

    def remove_node(self, node: int):
        """Izbaci prvo pojavljivanje datog cvora iz rute."""
        self.nodes.remove(node)                       # izbaci prvo pojavljivanje datog cvora

    def remove_at(self, position: int):
        """Izbaci cvor sa date pozicije."""
        self.nodes.pop(position)                      # izbaci cvor sa date pozicije

    def copy(self) -> 'Route':
        return Route(list(self.nodes))                # napravi nezavisnu kopiju rute

    def __repr__(self):
        return f"Route({self.nodes})"                 # tekstualni prikaz, npr. Route([0, 3, 0])

    def __len__(self):
        return len(self.nodes)                        # broj cvorova u ruti


class Solution:
    """Kompletno resenje: skup ruta koji pokriva sve musterije."""
    # RESENJE = lista ruta. Cuva i kes (cache) ocene da se ne racuna stalno iznova.

    def __init__(self, instance: Instance, routes: List[Route] = None):
        self.instance = instance                      # na koju instancu se odnosi
        self.routes = routes if routes is not None else []  # lista ruta (podrazumevano prazna)
        self._total_distance = None                   # kes: ukupna kilometraza (None = jos neizracunato)
        self._feasible = None                         # kes: izvodljivost
        self._evaluation = None                       # kes: cela ocena

    def invalidate_cache(self):
        # Pozvati kad se resenje promeni -> brise kes da se ocena ponovo izracuna.
        self._total_distance = None
        self._feasible = None
        self._evaluation = None

    @property
    def n_routes(self) -> int:
        return len(self.routes)                       # broj ruta (= broj koriscenih vozila)

    def add_route(self, route: Route):
        self.routes.append(route)                     # dodaj rutu...
        self.invalidate_cache()                       # ...i ponisti kes

    def remove_empty_routes(self):
        self.routes = [r for r in self.routes if not r.is_empty]  # izbaci prazne rute [0,0]
        self.invalidate_cache()

    def copy(self) -> 'Solution':
        new_sol = Solution(self.instance, [r.copy() for r in self.routes])  # duboka kopija svih ruta
        return new_sol

    def get_all_customers(self) -> set:
        """Vrati skup svih musterija koje su u resenju."""
        customers = set()
        for route in self.routes:                     # prodji kroz sve rute...
            for n in route.nodes:                     # i sve cvorove u njima...
                if self.instance.is_customer(n):
                    customers.add(n)                  # skupi sve opsluzene musterije
        return customers

    def get_unserved_customers(self) -> set:
        """Vrati musterije koje jos nisu ni u jednoj ruti."""
        served = self.get_all_customers()             # opsluzene
        all_customers = set(self.instance.customer_indices)  # sve musterije
        return all_customers - served                 # razlika = jos NEopsluzene

    def __repr__(self):
        eval_info = self.evaluate()                   # uzmi ocenu (kilometraza, izvodljivost)
        return (f"Solution(routes={self.n_routes}, "
                f"dist={eval_info['total_distance']:.2f}, "
                f"feasible={eval_info['feasible']})")


def evaluate_route(route: Route, instance: Instance) -> dict:
    """
    Oceni jednu rutu: izracunaj kilometrazu, trajanje i izvodljivost po gorivu.
    NE proverava kapacitet punionica (to je globalno ogranicenje, vidi reschedule).

    Vraca recnik sa:
    - distance: ukupna predjena kilometraza
    - duration: ukupno trajanje rute (voznja + usluga + dopuna)
    - fuel_feasible: True ako su uslovi goriva ispostovani
    - duration_feasible: True ako je trajanje <= Tmax
    - arrival_times: vremena dolaska u svaki cvor
    - fuel_levels: nivoi goriva pri dolasku u svaki cvor
    """
    # ===== OCENA JEDNE RUTE: prolazak cvor-po-cvor uz pracenje vremena i goriva =====
    nodes = route.nodes
    n = len(nodes)

    if n <= 2:
        # Prazna ruta ([0,0]) -> sve nule, pun rezervoar, izvodljivo.
        return {
            'distance': 0.0,
            'duration': 0.0,
            'fuel_feasible': True,
            'duration_feasible': True,
            'arrival_times': [0.0, 0.0],
            'fuel_levels': [instance.tank_capacity, instance.tank_capacity],
            'departure_times': [0.0, 0.0],
        }

    total_distance = 0.0                              # zbir svih deonica
    arrival_times = [0.0] * n                         # vreme DOLASKA u svaki cvor
    departure_times = [0.0] * n                       # vreme POLASKA iz svakog cvora
    fuel_levels = [0.0] * n                            # nivo goriva pri dolasku u svaki cvor

    # U depou krecemo sa PUNIM rezervoarom (Q)
    fuel_levels[0] = instance.tank_capacity
    departure_times[0] = instance.p_start              # pocetno vreme priprema u depou

    fuel_feasible = True                               # pretpostavka: gorivo je OK dok se ne dokaze suprotno

    for i in range(1, n):                              # za svaki sledeci cvor u ruti...
        prev = nodes[i - 1]                            # prethodni cvor
        curr = nodes[i]                                # trenutni cvor

        d = instance.dist(prev, curr)                  # rastojanje deonice
        t = instance.travel_time(prev, curr)           # vreme putovanja deonice
        total_distance += d                            # dodaj na ukupnu kilometrazu

        # Vreme dolaska = polazak iz prethodnog cvora + voznja
        arrival_times[i] = departure_times[i - 1] + t

        # Nivo goriva pri dolasku
        fuel_consumed = instance.consumption_rate * d  # potroseno gorivo = r * rastojanje
        fuel_levels[i] = fuel_levels[i - 1] - fuel_consumed  # gorivo se SMANJI za potroseno

        if fuel_levels[i] < -1e-6:                     # ako je gorivo PALO ISPOD NULE...
            fuel_feasible = False                      # ...ruta je nedozvoljena (vozilo bi stalo)

        # Vreme polaska iz trenutnog cvora
        if instance.is_customer(curr):
            departure_times[i] = arrival_times[i] + instance.get_service_time(curr)  # musterija: + usluga
        elif instance.is_station(curr):
            # Dopuna do punog
            departure_times[i] = arrival_times[i] + instance.refueling_time  # punionica: + vreme dopune
            fuel_levels[i] = instance.tank_capacity    # ...i rezervoar se NAPUNI na Q
        elif instance.is_depot(curr):
            departure_times[i] = arrival_times[i]      # depo: nema zadrzavanja
        else:
            departure_times[i] = arrival_times[i]

    duration = arrival_times[-1]                       # trajanje rute = vreme dolaska u krajnji depo

    return {
        'distance': total_distance,
        'duration': duration,
        'fuel_feasible': fuel_feasible,
        'duration_feasible': duration <= instance.t_max + 1e-6,  # da li je trajanje <= T_max
        'arrival_times': arrival_times,
        'departure_times': departure_times,
        'fuel_levels': fuel_levels,
    }


def reschedule(instance: Instance, solution: 'Solution') -> Tuple[bool, Dict[int, float]]:
    """
    Reschedule procedura: proveri i resi konflikte kapaciteta punionica.

    Zasnovano na Algoritmu A1 iz Xu et al. (2025), dodatak.

    Za svaku punionicu skupi sva vozila koja je posecuju sa njihovim vremenima dolaska.
    Ako se preklapa vise od eta_s vozila, odlozi (zadrzi) neka vozila.

    Vraca:
        (feasible, waiting_times):
        - feasible: True ako je raspored moguc unutar Tmax
        - waiting_times: recnik (route_idx, position) -> vreme cekanja
    """
    # ===== RED NA PUMPI: ako vise vozila dodje na istu stanicu nego sto ima pumpi,
    #       neka vozila CEKAJU. Ovo je deo specifican za nas problem (kapacitet punionica). =====
    waiting_times = {}                                 # rezultat: (ruta, pozicija) -> koliko ceka

    # Prvo oceni svaku rutu (bez cekanja), pa skupi posete punionicama po stanici
    route_evals = []
    for route in solution.routes:
        route_evals.append(evaluate_route(route, instance))  # prvo oceni svaku rutu (bez cekanja)

    # Grupisi posete punionicama po stanici
    station_visits = {}  # indeks_stanice -> lista (route_idx, pos_in_route, arrival_time)

    for r_idx, route in enumerate(solution.routes):    # prodji kroz sve rute...
        eval_data = route_evals[r_idx]
        for pos, node in enumerate(route.nodes):       # i sve cvorove...
            if instance.is_station(node):              # ako je cvor punionica...
                station_idx = node
                if station_idx not in station_visits:
                    station_visits[station_idx] = []
                station_visits[station_idx].append(    # zabelezi posetu: (ruta, pozicija, vreme dolaska)
                    (r_idx, pos, eval_data['arrival_times'][pos])
                )

    # Za svaku stanicu proveri kapacitet i resi konflikte
    feasible = True

    for station_idx, visits in station_visits.items():  # za svaku punionicu posebno...
        capacity = instance.get_station_capacity(station_idx)  # broj pumpi (eta)
        refuel_time = instance.refueling_time          # koliko traje jedna dopuna

        if len(visits) <= capacity:
            continue                                   # ako je poseta <= broj pumpi -> nema reda

        # Poredaj posete po vremenu dolaska
        visits_sorted = sorted(visits, key=lambda x: x[2])  # poredaj posete po vremenu dolaska

        # Primeni algoritam preraspodele
        arrival_times = [v[2] for v in visits_sorted]  # samo vremena dolaska
        n_vehicles = len(arrival_times)

        # Preraspodela: najvise "capacity" vozila se puni istovremeno
        adjusted_arrivals = list(arrival_times)        # kopija koju cemo "gurati" da resimo preklapanja
        changed = True
        max_iterations = 1000                          # zastita od beskonacne petlje
        iteration = 0

        while changed and iteration < max_iterations:  # ponavljaj dok ima preklapanja...
            changed = False
            iteration += 1
            for i in range(n_vehicles - 1):
                num_overlap = 0
                start_i = adjusted_arrivals[i]         # pocetak dopune vozila i
                end_i = start_i + refuel_time          # kraj dopune vozila i

                for j in range(i + 1, n_vehicles):
                    start_j = adjusted_arrivals[j]
                    end_j = start_j + refuel_time

                    if start_i < end_j and start_j < end_i:  # da li se intervali i i j PREKLAPAJU?
                        num_overlap += 1

                        if num_overlap > capacity - 1:  # vise preklapanja nego sto ima pumpi -> KONFLIKT
                            # Konflikt: odlozi jedno vozilo na kraj drugog
                            if end_j <= end_i:
                                adjusted_arrivals[i] = end_j   # pomeri i da pocne kad j zavrsi
                            else:
                                adjusted_arrivals[j] = end_i   # ili pomeri j da pocne kad i zavrsi
                            changed = True
                            break
                if changed:
                    break

        # Izracunaj vremena cekanja
        for k, (r_idx, pos, orig_arrival) in enumerate(visits_sorted):
            wt = adjusted_arrivals[k] - orig_arrival   # cekanje = novo vreme - prvobitno vreme dolaska
            if wt > 1e-6:
                waiting_times[(r_idx, pos)] = wt        # zapamti koliko ko ceka

    # Proveri da li cekanje gura neku rutu preko Tmax
    for r_idx, route in enumerate(solution.routes):    # da li cekanje gura neku rutu preko T_max?
        total_wait = sum(
            wt for (ri, pos), wt in waiting_times.items() if ri == r_idx  # ukupno cekanje te rute
        )
        eval_data = route_evals[r_idx]
        if eval_data['duration'] + total_wait > instance.t_max + 1e-6:
            feasible = False                           # ako da -> resenje NIJE izvodljivo po kapacitetu

    return feasible, waiting_times


def evaluate_solution(solution: 'Solution') -> dict:
    """
    Puna ocena resenja, ukljucujuci i ogranicenja kapaciteta punionica.

    Vraca recnik sa:
    - total_distance: zbir kilometraza svih ruta
    - route_evaluations: lista ocena po ruti
    - all_customers_served: True ako je svaka musterija opsluzena tacno jednom
    - vehicle_count_feasible: True ako je broj ruta <= broj vozila
    - fuel_feasible: True ako sve rute postuju gorivo
    - duration_feasible: True ako sve rute postuju trajanje
    - capacity_feasible: True ako se kapacitet punionica moze ispostovati
    - feasible: True ako su SVA ogranicenja zadovoljena
    - waiting_times: recnik vremena cekanja iz Reschedule
    """
    # ===== PUNA OCENA: kilometraza + SVE provere izvodljivosti =====
    instance = solution.instance

    # Oceni svaku rutu
    route_evals = []
    total_distance = 0.0
    fuel_feasible = True
    duration_feasible = True

    for route in solution.routes:                      # oceni svaku rutu posebno...
        eval_data = evaluate_route(route, instance)
        route_evals.append(eval_data)
        total_distance += eval_data['distance']        # saberi kilometrazu
        if not eval_data['fuel_feasible']:
            fuel_feasible = False                       # ako bar jedna ruta ima problem s gorivom...
        if not eval_data['duration_feasible']:
            duration_feasible = False                   # ...ili s trajanjem -> obelezi

    # Proveri da li je svaka musterija opsluzena tacno jednom
    customer_count = {}                                # broji koliko puta je svaka musterija opsluzena
    for route in solution.routes:
        for node in route.nodes:
            if instance.is_customer(node):
                customer_count[node] = customer_count.get(node, 0) + 1

    all_served = (set(customer_count.keys()) == set(instance.customer_indices))  # da li su SVE opsluzene?
    no_duplicates = all(v == 1 for v in customer_count.values())  # da li je svaka tacno JEDNOM?
    customers_ok = all_served and no_duplicates

    # Proveri broj vozila
    vehicle_ok = len(solution.routes) <= instance.n_vehicles  # broj ruta <= broj vozila?

    # Proveri kapacitet punionica preko Reschedule
    capacity_feasible, waiting_times = reschedule(instance, solution)  # provera reda na pumpi

    # Ponovo proveri trajanje sa vremenima cekanja
    if waiting_times:                                  # ako ima cekanja, ponovo proveri trajanje...
        for r_idx, route in enumerate(solution.routes):
            total_wait = sum(
                wt for (ri, pos), wt in waiting_times.items() if ri == r_idx
            )
            if route_evals[r_idx]['duration'] + total_wait > instance.t_max + 1e-6:
                duration_feasible = False              # cekanje moze da obori trajanje preko T_max

    feasible = (customers_ok and vehicle_ok and fuel_feasible  # IZVODLJIVO = svi uslovi zadovoljeni
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


# Zakaci metodu evaluate() na klasu Solution
# (tehnicki trik: dodajemo metodu evaluate() klasi Solution spolja, sa kesiranjem)
def _solution_evaluate(self) -> dict:
    if self._evaluation is None:                       # ako ocena jos nije izracunata...
        self._evaluation = evaluate_solution(self)     # ...izracunaj je i zapamti (kes)
    return self._evaluation                            # vrati (mozda kesiranu) ocenu

Solution.evaluate = _solution_evaluate                 # zakaci funkciju kao metodu Solution.evaluate


def compute_total_distance(solution: Solution) -> float:
    """Brzo racunanje samo ukupne kilometraze (bez pune ocene)."""
    # Brzo sabiranje SAMO kilometraze (bez svih provera) -- koristi se cesto u pretrazi.
    td = 0.0
    for route in solution.routes:
        for i in range(len(route.nodes) - 1):
            td += solution.instance.dist(route.nodes[i], route.nodes[i + 1])  # saberi sve deonice
    return td


def compute_penalty(solution: Solution, w_duration: float = 1.0,
                     w_fuel: float = 1.0, w_capacity: float = 1.0) -> float:
    """
    Izracunaj kaznu za prekrsaje ogranicenja.
    Koristi se u metaheuristici da vodi pretragu kroz nedozvoljena podrucja.

    Vraca ukupnu vrednost kazne (0 ako je resenje izvodljivo).
    """
    # ===== KAZNA ZA PREKRSAJE: 0 ako je sve OK; inace srazmerna velicini prekrsaja.
    #       Omogucava algoritmu da PRIVREMENO prodje kroz nedozvoljena resenja. =====
    instance = solution.instance
    penalty = 0.0

    for route in solution.routes:
        eval_data = evaluate_route(route, instance)

        # Prekrsaj trajanja
        if eval_data['duration'] > instance.t_max:     # ako ruta traje duze od T_max...
            penalty += w_duration * (eval_data['duration'] - instance.t_max)  # kazni za visak vremena

        # Prekrsaj goriva: proveri svaki nivo goriva
        for fl in eval_data['fuel_levels']:            # za svaki nivo goriva u ruti...
            if fl < -1e-6:
                penalty += w_fuel * abs(fl)            # kazni za koliko je gorivo otislo u minus

    # Prekrsaj kapaciteta
    capacity_feasible, waiting_times = reschedule(instance, solution)  # provera reda na pumpi
    if not capacity_feasible:
        # Zbir viska cekanja
        total_excess_wait = sum(waiting_times.values()) if waiting_times else 0
        penalty += w_capacity * (total_excess_wait + 1.0)  # kazni za prekoracenje kapaciteta

    return penalty


# =============================================================================
# ADAPTIVNI KAZNENI FAKTOR -- deli ga CELA pretraga (cena + lokalna pretraga).
# GVNS ga povecava kad pretraga dugo luta po neizvodljivom prostoru (da je gurne
# nazad u izvodljivo), a smanjuje kad je vecina resenja izvodljiva.
# =============================================================================
PENALTY_WEIGHT = 1000.0          # mnozilac kazne u ceni resenja (pocetna vrednost)


def set_penalty_weight(w: float):
    """Postavi globalni kazneni faktor (koristi ga GVNS za adaptaciju)."""
    global PENALTY_WEIGHT
    PENALTY_WEIGHT = float(w)


def penalized_cost(solution: Solution) -> float:
    """
    Jedinstvena CENA resenja = kilometraza + PENALTY_WEIGHT * kazna.
    Sve (i GVNS i lokalna pretraga) racunaju cenu OVDE, da adaptacija kazne
    deluje dosledno na celu pretragu.
    """
    return compute_total_distance(solution) + PENALTY_WEIGHT * compute_penalty(solution)
