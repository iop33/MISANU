"""
Strukture okolina (poteza) za GVNS metaheuristiku primenjenu na GrVRP-PCAFS.

Okoline (poredjane po slozenosti):
1. Relocate: premesti jednu musteriju na drugu poziciju
2. Swap: zameni dve musterije izmedju ruta
3. Or-opt: premesti niz od 2-3 uzastopne musterije
4. 2-opt unutar rute: obrni segment unutar rute
5. 2-opt* izmedju ruta: zameni "repove" dve rute
6. Premestanje punionice: promeni poziciju ili izbor punionice
7. Ubacivanje/izbacivanje punionice: dodaj ili izbaci posetu punionici

Svaka okolina ima:
- shake(): nasumican sused (za diversifikaciju)
- local_search(): najbolji sused (za intenzifikaciju)
"""

# =============================================================================
# neighborhoods.py  -- "POTEZI" (nacini da se ruta malo izmeni)
# -----------------------------------------------------------------------------
# Svaki potez pravi "suseda" trenutnog resenja. Dve uloge:
#   - *_shake : NASUMICAN potez (za bekstvo iz lokalnog minimuma -- diversifikacija)
#   - *_best  : NAJBOLJI potez te vrste (za doterivanje -- intenzifikacija)
# Na dnu fajla su dve liste: SHAKE_NEIGHBORHOODS i LOCAL_SEARCH_NEIGHBORHOODS.
# =============================================================================

import random                                          # nasumicni izbor (za shake poteze)
import copy                                             # (uvezeno)
import math                                             # (uvezeno)
from typing import List, Tuple, Optional                # oznake tipova
from instance import Instance
from solution import (Route, Solution, evaluate_route, compute_total_distance,
                      evaluate_solution, compute_penalty, penalized_cost)
from construction import find_nearest_feasible_station, insert_station_if_needed


def _route_distance(route: Route, instance: Instance) -> float:
    """Brzo racunanje kilometraze jedne rute."""
    d = 0.0
    for i in range(len(route.nodes) - 1):
        d += instance.dist(route.nodes[i], route.nodes[i + 1])
    return d


def _check_fuel_feasibility(route: Route, instance: Instance) -> bool:
    """Brza provera izvodljivosti po gorivu."""
    # Da li gorivo nigde ne padne ispod nule duz rute?
    fuel = instance.tank_capacity
    for i in range(1, len(route.nodes)):
        d = instance.dist(route.nodes[i-1], route.nodes[i])
        fuel -= instance.consumption_rate * d           # potrosi gorivo
        if fuel < -1e-6:
            return False                                # ostalo bez goriva -> nedozvoljeno
        if instance.is_station(route.nodes[i]):
            fuel = instance.tank_capacity               # punionica -> dopuna
    return True


def _check_duration_feasibility(route: Route, instance: Instance) -> bool:
    """Brza provera izvodljivosti po trajanju."""
    # Da li ruta staje u T_max?
    t = instance.p_start
    for i in range(1, len(route.nodes)):
        t += instance.travel_time(route.nodes[i-1], route.nodes[i])  # voznja
        if instance.is_customer(route.nodes[i]):
            t += instance.service_time_customer         # + usluga kod musterije
        elif instance.is_station(route.nodes[i]):
            t += instance.refueling_time                # + dopuna na punionici
    return t <= instance.t_max + 1e-6


def _is_route_feasible(route: Route, instance: Instance) -> bool:
    """Da li je ruta izvodljiva (gorivo + trajanje)?"""
    return _check_fuel_feasibility(route, instance) and _check_duration_feasibility(route, instance)


# ==================== OKOLINA 1: RELOCATE ====================
# RELOCATE = premesti JEDNU musteriju na drugo mesto.

def relocate_shake(solution: Solution) -> Solution:
    """Nasumicno premesti jednu musteriju na drugu poziciju."""
    sol = solution.copy()
    instance = sol.instance

    # Popisi sve pozicije musterija
    positions = []
    for r_idx, route in enumerate(sol.routes):
        for pos in range(1, len(route.nodes) - 1):
            if instance.is_customer(route.nodes[pos]):
                positions.append((r_idx, pos))

    if not positions:
        return sol

    # Izaberi nasumicnu musteriju
    r_idx, pos = random.choice(positions)
    customer = sol.routes[r_idx].nodes[pos]
    sol.routes[r_idx].remove_at(pos)                    # izvadi je iz rute

    # Ocisti punionice koje vise ne trebaju
    _clean_stations(sol.routes[r_idx], instance)

    # Izaberi nasumicnu poziciju ubacivanja (ista ili druga ruta)
    possible_routes = list(range(len(sol.routes)))
    if len(sol.routes) < instance.n_vehicles:
        possible_routes.append(-1)                      # opcija: otvori novu rutu

    target_r = random.choice(possible_routes)

    if target_r == -1:
        # Nova ruta
        new_route = Route([0, customer, 0])             # nova ruta samo sa tom musterijom
        new_route = insert_station_if_needed(new_route, instance)
        sol.routes.append(new_route)
    else:
        route = sol.routes[target_r]
        # Nasumicna pozicija (izmedju depo cvorova)
        insert_pos = random.randint(1, len(route.nodes) - 1)  # nasumicna pozicija ubacivanja
        route.insert(insert_pos, customer)

    sol.remove_empty_routes()
    sol.invalidate_cache()
    return sol


def relocate_best(solution: Solution) -> Solution:
    """Nadji najbolji relocate potez."""
    # NAJBOLJI relocate: probaj SVAKU musteriju na SVAKO mesto, vrati najjeftinije.
    instance = solution.instance
    best_sol = solution
    best_cost = penalized_cost(solution)                # cena = kilometraza + PENALTY_WEIGHT*kazna

    # Probaj sva premestanja musterija
    for r_idx, route in enumerate(solution.routes):
        for pos in range(1, len(route.nodes) - 1):
            if not instance.is_customer(route.nodes[pos]):
                continue

            customer = route.nodes[pos]

            # Probaj sve pozicije ubacivanja
            for t_idx, t_route in enumerate(solution.routes):
                for t_pos in range(1, len(t_route.nodes)):
                    if t_idx == r_idx and (t_pos == pos or t_pos == pos + 1):
                        continue                        # premestanje na isto mesto -> preskoci

                    # Napravi kandidata
                    sol = solution.copy()
                    sol.routes[r_idx].remove_at(pos)    # izvadi musteriju

                    # Ispravi ciljnu poziciju ako je ista ruta i pos se pomerio
                    actual_t_pos = t_pos
                    if t_idx == r_idx and t_pos > pos:
                        actual_t_pos -= 1               # ispravi indeks (jer se lista pomerila)

                    sol.routes[t_idx].insert(actual_t_pos, customer)  # ubaci na novo mesto
                    sol.remove_empty_routes()

                    cost = penalized_cost(sol)

                    if cost < best_cost - 1e-6:         # ako je jeftinije -> zapamti
                        best_cost = cost
                        best_sol = sol

    best_sol.invalidate_cache()
    return best_sol


# ==================== OKOLINA 2: SWAP ====================
# SWAP = zameni mesta dvema musterijama.

def swap_shake(solution: Solution) -> Solution:
    """Nasumicno zameni dve musterije."""
    sol = solution.copy()
    instance = sol.instance

    positions = []
    for r_idx, route in enumerate(sol.routes):          # popisi sve pozicije musterija
        for pos in range(1, len(route.nodes) - 1):
            if instance.is_customer(route.nodes[pos]):
                positions.append((r_idx, pos))

    if len(positions) < 2:
        return sol

    (r1, p1), (r2, p2) = random.sample(positions, 2)    # izaberi dve razlicite musterije

    # Zameni im mesta
    sol.routes[r1].nodes[p1], sol.routes[r2].nodes[p2] = \
        sol.routes[r2].nodes[p2], sol.routes[r1].nodes[p1]

    sol.invalidate_cache()
    return sol


def swap_best(solution: Solution) -> Solution:
    """Nadji najbolju zamenu (swap)."""
    # NAJBOLJA zamena: probaj svaki par musterija, vrati najjeftiniju zamenu.
    instance = solution.instance
    best_sol = solution
    best_cost = penalized_cost(solution)

    positions = []
    for r_idx, route in enumerate(solution.routes):
        for pos in range(1, len(route.nodes) - 1):
            if instance.is_customer(route.nodes[pos]):
                positions.append((r_idx, pos))

    for i in range(len(positions)):                     # za svaki par (i, j) musterija...
        for j in range(i + 1, len(positions)):
            r1, p1 = positions[i]
            r2, p2 = positions[j]

            sol = solution.copy()
            sol.routes[r1].nodes[p1], sol.routes[r2].nodes[p2] = \
                sol.routes[r2].nodes[p2], sol.routes[r1].nodes[p1]  # probna zamena

            cost = penalized_cost(sol)

            if cost < best_cost - 1e-6:
                best_cost = cost
                best_sol = sol

    best_sol.invalidate_cache()
    return best_sol


# ==================== OKOLINA 3: OR-OPT ====================
# OR-OPT = premesti SEGMENT od 1-3 uzastopne musterije na drugo mesto.

def or_opt_shake(solution: Solution) -> Solution:
    """Nasumicno premesti segment od 1-3 musterije na drugu poziciju."""
    sol = solution.copy()
    instance = sol.instance

    # Izaberi nasumicnu rutu koja ima musterije
    routes_with_customers = [
        (r_idx, r) for r_idx, r in enumerate(sol.routes)
        if len(r.customers(instance)) >= 1
    ]
    if not routes_with_customers:
        return sol

    r_idx, route = random.choice(routes_with_customers)  # nasumicna ruta sa musterijama

    # Nadji pozicije musterija
    cust_positions = [
        p for p in range(1, len(route.nodes) - 1)
        if instance.is_customer(route.nodes[p])
    ]
    if not cust_positions:
        return sol

    seg_len = random.choice([1, 2, 3])                  # duzina segmenta: 1, 2 ili 3
    seg_len = min(seg_len, len(cust_positions))

    start_idx = random.randint(0, len(cust_positions) - seg_len)
    segment_positions = cust_positions[start_idx:start_idx + seg_len]

    # Izvadi cvorove segmenta
    segment = [route.nodes[p] for p in segment_positions]

    # Izbaci segment iz rute (unazad, da indeksi ostanu tacni)
    for p in sorted(segment_positions, reverse=True):
        sol.routes[r_idx].remove_at(p)

    # Ubaci segment na nasumicnu poziciju u nasumicnoj ruti
    possible_routes = list(range(len(sol.routes)))
    target_r = random.choice(possible_routes)
    target_route = sol.routes[target_r]
    insert_pos = random.randint(1, max(1, len(target_route.nodes) - 1))

    for k, node in enumerate(segment):                  # ubaci ceo segment na novo mesto
        target_route.insert(insert_pos + k, node)

    sol.remove_empty_routes()
    sol.invalidate_cache()
    return sol


# ==================== OKOLINA 4: 2-OPT UNUTAR RUTE ====================
# 2-OPT (unutar rute) = OBRNI deo rute (raspetljava ukrstanja).

def two_opt_intra_shake(solution: Solution) -> Solution:
    """Nasumicno primeni 2-opt unutar jedne rute (obrni segment)."""
    sol = solution.copy()
    instance = sol.instance

    non_empty = [r_idx for r_idx, r in enumerate(sol.routes)
                 if len(r.nodes) > 3]                   # rute dovoljno duge da ima sta da se obrne
    if not non_empty:
        return sol

    r_idx = random.choice(non_empty)
    route = sol.routes[r_idx]

    # Izaberi dve tacke (ne depo)
    n = len(route.nodes)
    if n <= 3:
        return sol

    i = random.randint(1, n - 3)                        # dve nasumicne tacke unutar rute
    j = random.randint(i + 1, n - 2)

    # Obrni segment izmedju i i j
    route.nodes[i:j+1] = route.nodes[i:j+1][::-1]

    sol.invalidate_cache()
    return sol


def two_opt_intra_best(solution: Solution) -> Solution:
    """Nadji najbolji 2-opt potez unutar rute."""
    # NAJBOLJI 2-opt: probaj sve parove (i,j) u svim rutama, vrati najbolju obrnutu varijantu.
    instance = solution.instance
    best_sol = solution
    best_cost = penalized_cost(solution)

    for r_idx, route in enumerate(solution.routes):
        n = len(route.nodes)
        if n <= 3:
            continue

        for i in range(1, n - 2):
            for j in range(i + 1, n - 1):
                sol = solution.copy()
                sol.routes[r_idx].nodes[i:j+1] = sol.routes[r_idx].nodes[i:j+1][::-1]  # probno obrtanje

                cost = penalized_cost(sol)

                if cost < best_cost - 1e-6:
                    best_cost = cost
                    best_sol = sol

    best_sol.invalidate_cache()
    return best_sol


# ==================== OKOLINA 5: 2-OPT* IZMEDJU RUTA ====================
# 2-OPT* (izmedju ruta) = razmeni "repove" dve rute.

def two_opt_star_shake(solution: Solution) -> Solution:
    """Nasumicno zameni "repove" dve rute."""
    sol = solution.copy()
    instance = sol.instance

    if len(sol.routes) < 2:
        return sol                                      # treba bar dve rute

    r1_idx, r2_idx = random.sample(range(len(sol.routes)), 2)  # dve razlicite rute
    route1 = sol.routes[r1_idx]
    route2 = sol.routes[r2_idx]

    if len(route1.nodes) <= 2 or len(route2.nodes) <= 2:
        return sol

    # Izaberi mesta secenja
    cut1 = random.randint(1, len(route1.nodes) - 2)     # mesto secenja u ruti 1
    cut2 = random.randint(1, len(route2.nodes) - 2)     # mesto secenja u ruti 2

    # Zameni repove
    tail1 = route1.nodes[cut1:]                         # "rep" rute 1
    tail2 = route2.nodes[cut2:]                         # "rep" rute 2

    sol.routes[r1_idx].nodes = route1.nodes[:cut1] + tail2  # zameni repove
    sol.routes[r2_idx].nodes = route2.nodes[:cut2] + tail1

    sol.remove_empty_routes()
    sol.invalidate_cache()
    return sol


# ==================== OKOLINA 6: UBACI/IZBACI PUNIONICU ====================
# Potez nad PUNIONICAMA: ubaci / izbaci / promeni punionicu.

def station_change_shake(solution: Solution) -> Solution:
    """Nasumicno ubaci, izbaci ili promeni posetu punionici."""
    sol = solution.copy()
    instance = sol.instance

    action = random.choice(['insert', 'remove', 'change'])  # nasumicno izaberi akciju

    if action == 'insert':
        # Ubaci nasumicnu punionicu na nasumicnu poziciju
        non_empty = [r_idx for r_idx, r in enumerate(sol.routes) if len(r.nodes) > 2]
        if not non_empty:
            return sol
        r_idx = random.choice(non_empty)
        route = sol.routes[r_idx]
        station = random.choice(instance.station_indices)  # nasumicna punionica
        pos = random.randint(1, len(route.nodes) - 1)
        route.insert(pos, station)                      # ubaci je na nasumicno mesto

    elif action == 'remove':
        # Izbaci nasumicnu punionicu
        station_positions = []
        for r_idx, route in enumerate(sol.routes):      # popisi sve punionice u resenju
            for pos in range(1, len(route.nodes) - 1):
                if instance.is_station(route.nodes[pos]):
                    station_positions.append((r_idx, pos))

        if station_positions:
            r_idx, pos = random.choice(station_positions)
            sol.routes[r_idx].remove_at(pos)            # izbaci nasumicnu punionicu

    elif action == 'change':
        # Zameni jednu punionicu drugom
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
                sol.routes[r_idx].nodes[pos] = new_station  # zameni jednu punionicu drugom

    sol.invalidate_cache()
    return sol


# ==================== POMOCNE FUNKCIJE ====================

def _clean_stations(route: Route, instance: Instance):
    """Izbaci nepotrebne posete punionicama iz rute."""
    # Izbaci punionice koje VISE NE TREBAJU (ruta ostaje izvodljiva po gorivu i bez njih).
    changed = True
    while changed:
        changed = False
        for i in range(len(route.nodes) - 2, 0, -1):
            if instance.is_station(route.nodes[i]):
                # Probaj bez te punionice
                test_nodes = route.nodes[:i] + route.nodes[i+1:]
                test_route = Route(test_nodes)
                if _check_fuel_feasibility(test_route, instance):  # ako gorivo i dalje OK...
                    route.nodes = test_nodes              # ...izbaci je
                    changed = True
                    break


def fix_stations(solution: Solution) -> Solution:
    """Sredi posete punionicama u svim rutama: izbaci nepotrebne, dodaj gde treba."""
    # KLJUCNI POMOCNIK: izbaci SVE punionice pa ih vrati SAMO gde gorivo ponestaje.
    # Tako se odvaja "redosled musterija" od "gde se dopuniti".
    sol = solution.copy()
    instance = sol.instance

    for r_idx in range(len(sol.routes)):
        route = sol.routes[r_idx]

        # Prvo izbaci sve punionice (ostaju samo depo + musterije)
        customer_nodes = [0] + [n for n in route.nodes[1:-1] if instance.is_customer(n)] + [0]
        route = Route(customer_nodes)

        # Pa vrati punionice tamo gde trebaju
        route = insert_station_if_needed(route, instance)
        sol.routes[r_idx] = route

    sol.invalidate_cache()
    return sol


# ==================== KOLEKCIJE OKOLINA ====================

# Shaking okoline (za diversifikaciju u GVNS-u)
# 6 NASUMICNIH poteza -- koriste se u SHAKING fazi GVNS-a (bekstvo iz lokalnog minimuma).
SHAKE_NEIGHBORHOODS = [
    relocate_shake,
    swap_shake,
    or_opt_shake,
    two_opt_intra_shake,
    two_opt_star_shake,
    station_change_shake,
]

# Okoline lokalne pretrage (za intenzifikaciju u GVNS-u)
# 3 "NAJBOLJI" poteza -- koriste se u VND fazi (lokalno doterivanje resenja).
LOCAL_SEARCH_NEIGHBORHOODS = [
    relocate_best,
    swap_best,
    two_opt_intra_best,
]
