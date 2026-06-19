"""
Heuristike konstrukcije pocetnog resenja za GrVRP-PCAFS.

Implementira:
1. Pohlepnu konstrukciju (najblizi sused)
2. Konstrukciju po ustedama (Clarke & Wright)
"""

# =============================================================================
# construction.py  -- PRAVLJENJE POCETNOG RESENJA
# -----------------------------------------------------------------------------
# GVNS-u treba neka pocetna tacka. Ovde su dve klasicne heuristike:
#   1) greedy_construction  -- pohlepno: uvek idi na najblizu dostiznu musteriju
#   2) savings_construction -- Clarke & Wright: spajaj rute po "ustedi"
# Plus pomocne funkcije za ubacivanje punionica gde gorivo ponestaje.
# =============================================================================

import math                                          # (uvezeno)
import numpy as np                                    # (uvezeno)
from typing import List, Tuple, Optional, Set         # oznake tipova
from instance import Instance                         # podaci instance
from solution import Route, Solution, evaluate_route  # tipovi resenja + ocena rute


def find_nearest_feasible_station(instance: Instance, from_node: int,
                                   fuel_remaining: float) -> Optional[int]:
    """Nadji najblizu punionicu dostiznu sa trenutnim gorivom iz datog cvora."""
    best_station = None
    best_dist = float('inf')

    for s in instance.station_indices:                # za svaku punionicu...
        d = instance.dist(from_node, s)               # rastojanje do nje
        fuel_needed = instance.consumption_rate * d   # koliko goriva treba da stignemo
        if fuel_needed <= fuel_remaining + 1e-6 and d < best_dist:  # dostizna i bliza od dosadasnje?
            best_dist = d
            best_station = s
    return best_station


def can_return_to_depot(instance: Instance, from_node: int,
                        fuel_remaining: float) -> bool:
    """Proveri da li vozilo moze da se vrati u depo iz trenutnog cvora."""
    # Da li vozilo moze direktno do depoa sa preostalim gorivom?
    d = instance.dist(from_node, instance.depot)
    fuel_needed = instance.consumption_rate * d
    return fuel_needed <= fuel_remaining + 1e-6


def greedy_construction(instance: Instance) -> Solution:
    """
    Pohlepna konstrukcija (najblizi sused).

    Gradi rute jednu po jednu:
    1. Pocni iz depoa sa punim rezervoarom
    2. Idi na najblizu neposecenu musteriju koja je dostizna
    3. Ako je gorivo nisko, ubaci posetu punionici
    4. Ako se nijedna musterija ne moze dodati, zatvori rutu i otvori novu
    """
    # ===== POHLEPNA KONSTRUKCIJA: gradi rutu po rutu, uvek bira najblizu musteriju =====
    solution = Solution(instance)
    unserved = set(instance.customer_indices)          # skup jos neopsluzenih musterija

    while unserved and len(solution.routes) < instance.n_vehicles:  # dok ima posla i slobodnih vozila
        route_nodes = [0]                              # nova ruta krece iz depoa
        current_node = 0
        fuel = instance.tank_capacity                  # pun rezervoar
        current_time = instance.p_start                # pocetno vreme

        while unserved:
            # Nadji najblizu neopsluzenu musteriju
            best_customer = None
            best_dist = float('inf')

            for c in unserved:                         # razmotri svaku neopsluzenu musteriju c
                d = instance.dist(current_node, c)
                fuel_needed = instance.consumption_rate * d
                time_after = current_time + instance.travel_time(current_node, c) + instance.service_time_customer

                # Proveri mozemo li do c, pa onda nazad u depo (mozda preko punionice)
                d_to_depot = instance.dist(c, 0)
                fuel_after = fuel - fuel_needed         # gorivo nakon dolaska do c
                fuel_to_depot = instance.consumption_rate * d_to_depot
                time_to_depot = time_after + instance.travel_time(c, 0)

                can_reach = fuel_needed <= fuel + 1e-6  # mozemo li uopste do c?

                # Mozemo li posle da zatvorimo rutu (vratimo se)? Direktno ili preko punionice
                can_close = False
                if fuel_after >= fuel_to_depot - 1e-6 and time_to_depot <= instance.t_max + 1e-6:
                    can_close = True                    # DIREKTAN povratak moguc (gorivo + vreme OK)
                else:
                    # Ako ne direktno, probaj preko najblize punionice
                    for s in instance.station_indices:
                        d_c_s = instance.dist(c, s)
                        d_s_0 = instance.dist(s, 0)
                        fuel_c_s = instance.consumption_rate * d_c_s
                        if fuel_after >= fuel_c_s - 1e-6:  # imamo goriva do punionice?
                            time_via_s = (time_after + instance.travel_time(c, s)
                                         + instance.refueling_time
                                         + instance.travel_time(s, 0))
                            if time_via_s <= instance.t_max + 1e-6:  # vreme preko punionice OK?
                                # Posle dopune, ima li dovoljno goriva do depoa?
                                fuel_after_refuel = instance.tank_capacity
                                if fuel_after_refuel >= instance.consumption_rate * d_s_0 - 1e-6:
                                    can_close = True     # da, moze povratak preko punionice
                                    break

                if can_reach and can_close and d < best_dist:  # dostizna + moze povratak + najbliza
                    best_dist = d
                    best_customer = c

            if best_customer is None:
                break                                   # nijedna se ne moze dodati -> zatvori rutu

            # Proveri da li treba dopuna PRE posete musteriji
            d_to_cust = instance.dist(current_node, best_customer)
            fuel_needed = instance.consumption_rate * d_to_cust

            if fuel_needed > fuel + 1e-6:               # nemamo goriva ni do izabrane musterije?
                # Prvo dopuna - probaj sve punionice, izaberi najbolju
                best_station = None
                best_station_cost = float('inf')
                for s in instance.station_indices:      # nadji punionicu sa najmanjim obilaskom
                    d_cs = instance.dist(current_node, s)
                    d_sc = instance.dist(s, best_customer)
                    if instance.consumption_rate * d_cs <= fuel + 1e-6:        # stizemo do punionice?
                        if instance.consumption_rate * d_sc <= instance.tank_capacity + 1e-6:  # pa do musterije?
                            cost = d_cs + d_sc
                            if cost < best_station_cost:
                                best_station_cost = cost
                                best_station = s

                if best_station is None:
                    break                                # nijedna punionica ne pomaze -> zatvori rutu

                route_nodes.append(best_station)         # ubaci punionicu u rutu
                current_time += instance.travel_time(current_node, best_station) + instance.refueling_time
                fuel = instance.tank_capacity            # dopuna -> pun rezervoar
                current_node = best_station

                # Ponovo proveri da li je musterija dostizna
                d_to_cust = instance.dist(current_node, best_customer)
                fuel_needed = instance.consumption_rate * d_to_cust
                if fuel_needed > fuel + 1e-6:
                    break

            # Poseti musteriju
            route_nodes.append(best_customer)            # dodaj musteriju u rutu
            fuel -= instance.consumption_rate * d_to_cust  # potrosi gorivo
            current_time += instance.travel_time(current_node, best_customer) + instance.service_time_customer
            current_node = best_customer
            unserved.remove(best_customer)               # vise nije neopsluzena

            # Proveri da li treba dopuna za povratak u depo
            d_to_depot = instance.dist(current_node, 0)
            fuel_to_depot = instance.consumption_rate * d_to_depot

            if fuel < fuel_to_depot - 1e-6:              # nemamo goriva za povratak?
                # Dopuna pre povratka
                station = find_nearest_feasible_station(instance, current_node, fuel)
                if station is not None:
                    d_to_s = instance.dist(current_node, station)
                    d_s_0 = instance.dist(station, 0)
                    time_check = current_time + instance.travel_time(current_node, station) + instance.refueling_time + instance.travel_time(station, 0)
                    if time_check <= instance.t_max + 1e-6:
                        route_nodes.append(station)      # ubaci punionicu pre povratka
                        fuel = instance.tank_capacity
                        current_time += instance.travel_time(current_node, station) + instance.refueling_time
                        current_node = station

        # Zatvori rutu
        route_nodes.append(0)                            # zatvori rutu povratkom u depo
        route = Route(route_nodes)

        if not route.is_empty and len(route.customers(instance)) > 0:
            solution.add_route(route)                    # dodaj rutu ako ima bar jednu musteriju
        else:
            # Ako ruta nema musterija, prekini
            break

    solution.remove_empty_routes()

    # REZERVA: za musterije koje nisu stale, pravi pojedinacne rute [0,c,0] ili [0,s,c,s,0]
    while unserved and len(solution.routes) < instance.n_vehicles:
        c = unserved.pop()
        d_to = instance.dist(0, c)
        d_from = instance.dist(c, 0)
        if instance.consumption_rate * (d_to + d_from) <= instance.tank_capacity + 1e-6:  # ide bez dopune?
            t = instance.p_start + instance.travel_time(0, c) + instance.service_time_customer + instance.travel_time(c, 0)
            if t <= instance.t_max + 1e-6:
                solution.add_route(Route([0, c, 0]))     # prosta ruta depo->c->depo
                continue
        placed = False
        for s in instance.station_indices:               # inace probaj sa dopunom: depo->s->c->s->depo
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
    Konstrukcija po ustedama (Clarke & Wright).

    1. Pocni sa jednom rutom po musteriji: depo -> musterija -> depo
    2. Izracunaj ustede s(i,j) = d(i,0) + d(0,j) - d(i,j)
    3. Spajaj rute redom po najvecoj ustedi
    """
    # ===== CLARKE & WRIGHT: pocni sa rutom po musteriji, pa spajaj po najvecoj ustedi =====
    # Inicijalizacija: jedna ruta po musteriji
    routes = {}
    for c in instance.customer_indices:
        routes[c] = Route([0, c, 0])                     # svaka musterija dobije svoju rutu [0,c,0]

    # Izracunaj ustede
    savings = []
    for i in instance.customer_indices:
        for j in instance.customer_indices:
            if i != j:
                s = (instance.dist(i, 0) + instance.dist(0, j)   # USTEDA spajanja i i j:
                     - instance.dist(i, j))                       # d(i,0)+d(0,j)-d(i,j)
                savings.append((s, i, j))

    savings.sort(reverse=True, key=lambda x: x[0])       # sortiraj od NAJVECE ustede

    # Pridruzivanje: u kojoj je ruti svaka musterija?
    customer_route = {c: c for c in instance.customer_indices}  # u kojoj je ruti svaka musterija

    for saving, i, j in savings:
        if saving <= 0:
            break                                        # negativna usteda -> dalje nema smisla spajati

        route_i_key = customer_route[i]
        route_j_key = customer_route[j]

        if route_i_key == route_j_key:
            continue                                     # vec su u istoj ruti

        route_i = routes[route_i_key]
        route_j = routes[route_j_key]

        # Spajamo samo ako je i poslednja musterija svoje rute, a j prva u svojoj
        customers_i = route_i.customers(instance)
        customers_j = route_j.customers(instance)

        if not customers_i or not customers_j:
            continue

        if customers_i[-1] != i or customers_j[0] != j:  # spajamo samo kraj jedne sa pocetkom druge
            continue

        # Probaj spajanje: route_i + route_j (izbaci depo izmedju)
        merged_nodes = route_i.nodes[:-1] + route_j.nodes[1:]  # spoji rute (izbaci depo izmedju)
        merged_route = Route(merged_nodes)

        # Proveri izvodljivost spojene rute
        eval_data = evaluate_route(merged_route, instance)  # da li je spojena ruta izvodljiva?

        if (eval_data['fuel_feasible'] and eval_data['duration_feasible']
            and len(merged_route.customers(instance)) <= instance.n_customers):

            # Broj ruta nakon spajanja (informativno)
            total_routes_after = len(routes) - 1

            # Izvrsi spajanje
            routes[route_i_key] = merged_route           # zameni i-tu rutu spojenom...
            del routes[route_j_key]                      # ...i izbaci j-tu

            # Azuriraj pridruzivanja
            for c in customers_j:                        # azuriraj kojoj ruti pripadaju musterije iz j
                customer_route[c] = route_i_key

    solution = Solution(instance, list(routes.values()))
    solution.remove_empty_routes()
    return solution


def insert_station_if_needed(route: Route, instance: Instance) -> Route:
    """
    Naknadna obrada rute: ubaci posete punionici tamo gde gorivo ponestaje.
    Koristi ubacivanje najblize dostizne punionice.
    """
    # ===== Prodji rutu i UBACI PUNIONICU svuda gde bi gorivo ponestalo. =====
    nodes = route.nodes
    new_nodes = [nodes[0]]                                # nova ruta krece od istog pocetka (depo)
    fuel = instance.tank_capacity                         # pun rezervoar

    for i in range(1, len(nodes)):
        prev = new_nodes[-1]                              # poslednji cvor u novoj ruti
        curr = nodes[i]                                   # sledeci cvor koji hocemo da dodamo
        d = instance.dist(prev, curr)
        fuel_needed = instance.consumption_rate * d

        if fuel_needed > fuel + 1e-6:                     # nemamo goriva do sledeceg cvora?
            # Treba ubaciti punionicu PRE ovog cvora
            station = find_nearest_feasible_station(instance, prev, fuel)  # nadji dostiznu punionicu
            if station is not None:
                new_nodes.append(station)                 # ubaci punionicu PRE tog cvora
                fuel = instance.tank_capacity             # dopuna -> pun rezervoar
                # Ponovo izracunaj gorivo potrebno od punionice do curr
                fuel_needed = instance.consumption_rate * instance.dist(station, curr)

        new_nodes.append(curr)                            # dodaj trenutni cvor

        if instance.is_station(curr):
            fuel = instance.tank_capacity                 # ako je cvor punionica -> pun rezervoar
        else:
            fuel -= fuel_needed                           # inace potrosi gorivo

    return Route(new_nodes)
