"""
Opsta pretraga promenljivih okolina (GVNS) za GrVRP-PCAFS.

Reference:
- Hansen, P. and Mladenovic, N. (2001). Variable neighborhood search.
  European Journal of Operational Research, 130(3), 449-467.
- Mladenovic, N. and Hansen, P. (1997). Variable neighborhood search.
  Computers & Operations Research, 24(11), 1097-1100.

Struktura GVNS-a:
1. Napravi pocetno resenje
2. Ponavljaj dok ne istekne kriterijum zaustavljanja:
   a. Postavi k = 1
   b. Dok je k <= k_max:
      - Shaking: napravi nasumicnog suseda x' u N_k(x)
      - Lokalna pretraga: primeni VND na x' -> x''
      - Pomeri ili ne: ako x'' poboljsava x, postavi x = x'' i k = 1
        inace, k = k + 1
3. Vrati najbolje nadjeno resenje
"""

# =============================================================================
# gvns.py  -- SAM ALGORITAM (mozak projekta)
# -----------------------------------------------------------------------------
# GVNS = General Variable Neighborhood Search (Hansen & Mladenovic).
# Drzi JEDNO resenje i stalno ga doteruje smenjivanjem "okolina" (vrsta poteza):
#   shaking (nasumican trzaj) -> fix_stations -> VND (lokalno doterivanje) -> odluka.
# Sadrzi i: solution_cost (cena), vnd, _perturb (restart), run_multiple (vise pokretanja).
# =============================================================================

import time                                            # merenje vremena (vremensko ogranicenje)
import random                                           # nasumicnost (sa fiksnim seed-om)
import copy                                             # (uvezeno)
import numpy as np                                      # statistika + seed
from typing import List, Tuple, Optional, Callable      # oznake tipova
from instance import Instance
from solution import (Solution, Route, compute_total_distance,
                      evaluate_solution, compute_penalty, penalized_cost,
                      set_penalty_weight)
from construction import (greedy_construction, savings_construction,
                          insert_station_if_needed, scts_construction)
from neighborhoods import (
    SHAKE_NEIGHBORHOODS, LOCAL_SEARCH_NEIGHBORHOODS,
    fix_stations, relocate_best, swap_best, two_opt_intra_best,
)

# ---- Granice adaptivnog kaznenog faktora ----
PENALTY_INIT = 1000.0      # pocetni faktor (svaki run krece odavde)
PENALTY_MAX = 1.0e5        # gornja granica (da ne eksplodira)
PENALTY_MIN = 1000.0       # donja granica
ADAPT_PERIOD = 30          # na svakih toliko shaking-iteracija proveravamo udeo izvodljivih


def solution_cost(solution: Solution) -> float:
    """Cena resenja = kilometraza + (adaptivni) PENALTY_WEIGHT * kazna."""
    # Jedinstvena tacka racuna cene -> adaptacija kazne deluje i na GVNS i na lokalnu pretragu.
    return penalized_cost(solution)


def is_feasible(solution: Solution) -> bool:
    """Da li je resenje izvodljivo?"""
    eval_data = evaluate_solution(solution)
    return eval_data['feasible']                        # da li resenje postuje SVA ogranicenja?


def vnd(solution: Solution, neighborhoods: List[Callable],
        time_limit: float = None, start_time: float = None) -> Solution:
    """
    Spustanje kroz promenljive okoline (VND).

    Primenjuje okoline lokalne pretrage redom.
    Kad nadje poboljsanje, vraca se na prvu okolinu.
    """
    # ===== VND = LOKALNA PRETRAGA: redom primenjuj "najbolje" poteze; cim neki
    #       poboljsa resenje, vrati se na PRVI potez. Staje kad nista vise ne pomaze. =====
    current = solution
    current_cost = solution_cost(current)
    k = 0

    while k < len(neighborhoods):                       # prolazi kroz okoline (relocate, swap, 2-opt)
        if time_limit and start_time:
            if time.time() - start_time >= time_limit:
                break                                   # stani ako je isteklo vreme

        # Primeni lokalnu pretragu u okolini k
        neighbor = neighborhoods[k](current)            # nadji NAJBOLJI potez okoline k
        neighbor_cost = solution_cost(neighbor)

        if neighbor_cost < current_cost - 1e-6:         # ako je bolje...
            current = neighbor
            current_cost = neighbor_cost
            k = 0                                       # ...prihvati i vrati se na prvu okolinu
        else:
            k += 1                                      # inace probaj sledecu okolinu

    return current


def gvns(instance: Instance,
         time_limit: float = 600.0,
         max_iterations: int = 10000,
         max_no_improve: int = 500,
         k_max: int = 6,
         seed: int = None,
         verbose: bool = True) -> Tuple[Solution, dict]:
    """
    Opsta pretraga promenljivih okolina (GVNS) za GrVRP-PCAFS.

    Parametri:
    -----------
    instance : Instance
        Instanca problema
    time_limit : float
        Maksimalno vreme u sekundama
    max_iterations : int
        Maksimalan broj iteracija
    max_no_improve : int
        Najvise iteracija bez poboljsanja pre restarta
    k_max : int
        Najveci indeks okoline za shaking
    seed : int
        Seme slucajnosti (za ponovljivost)
    verbose : bool
        Ispisuj napredak

    Vraca:
    --------
    (best_solution, stats): najbolje nadjeno resenje i statistiku algoritma
    """
    # ===== GLAVNA FUNKCIJA: nadji najbolje resenje za datu instancu u datom vremenu. =====
    if seed is not None:
        random.seed(seed)                               # fiksiraj slucajnost -> ponovljivost
        np.random.seed(seed)

    set_penalty_weight(PENALTY_INIT)                    # RESETUJ kazneni faktor na pocetak svakog run-a
    penalty_weight = PENALTY_INIT                       # lokalna kopija (za adaptaciju)

    start_time = time.time()                            # pocni stopericu

    # Pracenje statistike
    stats = {                                           # statistika koja se vraca na kraju
        'iterations': 0,
        'improvements': 0,
        'time_to_best': 0.0,
        'best_costs': [],
        'feasible_found': False,
    }

    # ---- Pracenje najboljeg IZVODLJIVOG resenja (vraca se na kraju) ----
    best_feasible = None                                # najbolje IZVODLJIVO resenje do sada
    best_feasible_dist = float('inf')                   # njegova kilometraza
    feas_window = []                                    # poslednjih nekoliko (ne)izvodljivih -> za adaptaciju kazne

    def _track_feasible(sol):
        """Ako je 'sol' izvodljivo i krace od dosadasnjeg najboljeg izvodljivog -> zapamti ga."""
        nonlocal best_feasible, best_feasible_dist
        ev = evaluate_solution(sol)
        if ev['feasible']:
            if ev['total_distance'] < best_feasible_dist - 1e-9:
                best_feasible = sol.copy()
                best_feasible_dist = ev['total_distance']
            return True
        return False

    # Napravi pocetno resenje
    if verbose:
        print("Generating initial solution...")

    # Tri razlicita pocetna resenja -> uzmi najbolje (i zapamti ako je neko izvodljivo)
    candidates = [
        fix_stations(greedy_construction(instance)),    # 1: pohlepno (najblizi sused)
        fix_stations(savings_construction(instance)),   # 2: Clarke & Wright (ustede)
        scts_construction(instance),                    # 3: SCTS podela giant tour-a po Tmax
    ]
    best_solution = min(candidates, key=solution_cost)  # uzmi najbolje od tri pocetna
    for cand in candidates:
        _track_feasible(cand)                           # zapamti ako je neki vec izvodljiv

    best_cost = solution_cost(best_solution)            # cena najboljeg do sada
    current = best_solution.copy()                      # trenutno resenje (radna kopija)
    current_cost = best_cost

    if verbose:
        eval_data = evaluate_solution(best_solution)
        print(f"Initial solution: distance={eval_data['total_distance']:.2f}, "
              f"routes={eval_data['n_routes']}, feasible={eval_data['feasible']}")

    stats['best_costs'].append(best_cost)               # zapamti pocetnu cenu

    no_improve_count = 0                                # brojac iteracija bez napretka
    iteration = 0

    shake_neighborhoods = SHAKE_NEIGHBORHOODS[:k_max]   # uzmi prvih k_max shaking okolina
    ls_neighborhoods = LOCAL_SEARCH_NEIGHBORHOODS       # okoline za lokalnu pretragu (VND)

    while iteration < max_iterations:                   # GLAVNA PETLJA
        elapsed = time.time() - start_time
        if elapsed >= time_limit:
            break                                       # stani kad istekne vreme

        k = 0
        while k < len(shake_neighborhoods):             # prolazi kroz shaking okoline (N_0..N_k_max)
            elapsed = time.time() - start_time
            if elapsed >= time_limit:
                break

            iteration += 1
            stats['iterations'] = iteration

            # ---- SHAKING ----
            neighbor = shake_neighborhoods[k](current)  # NASUMICAN trzaj u okolini k (diversifikacija)

            # Sredi punionice posle trzaja
            neighbor = fix_stations(neighbor)

            # ---- LOKALNA PRETRAGA (VND) ----
            remaining_time = time_limit - elapsed
            neighbor = vnd(neighbor, ls_neighborhoods,  # lokalno doteraj (intenzifikacija)
                          time_limit=remaining_time, start_time=time.time())

            neighbor_cost = solution_cost(neighbor)

            # Zapamti ako je 'neighbor' izvodljivo (i krace) + evidentiraj za adaptaciju kazne
            feas_now = _track_feasible(neighbor)
            feas_window.append(feas_now)
            if best_feasible is not None:
                stats['feasible_found'] = True

            # ---- POMERI ILI NE ----
            if neighbor_cost < current_cost - 1e-6:     # da li je novo resenje BOLJE?
                current = neighbor
                current_cost = neighbor_cost
                k = 0                                   # da -> prihvati i vrati se na prvu okolinu
                no_improve_count = 0

                # Azuriraj najbolje (po ceni, sme biti i privremeno neizvodljivo)
                if current_cost < best_cost - 1e-6:     # da li je novi REKORD?
                    best_solution = current.copy()
                    best_cost = current_cost
                    stats['improvements'] += 1
                    stats['time_to_best'] = time.time() - start_time  # kada smo nasli najbolje

                    if verbose and iteration % 10 == 0:
                        eval_data = evaluate_solution(best_solution)
                        print(f"  Iter {iteration}: cost={best_cost:.2f}, "
                              f"dist={eval_data['total_distance']:.2f}, "
                              f"feasible={eval_data['feasible']}, "
                              f"time={time.time()-start_time:.1f}s")
            else:
                k += 1                                  # nije bolje -> probaj sledecu (vecu) okolinu
                no_improve_count += 1

            stats['best_costs'].append(best_cost)

            # ---- ADAPTACIJA KAZNE: gurni pretragu ka izvodljivom prostoru ----
            if iteration % ADAPT_PERIOD == 0 and len(feas_window) >= ADAPT_PERIOD:
                ratio = sum(feas_window[-ADAPT_PERIOD:]) / ADAPT_PERIOD  # udeo izvodljivih u prozoru
                old_w = penalty_weight
                if ratio < 0.15:                        # skoro nista izvodljivo -> POJACAJ kaznu
                    penalty_weight = min(PENALTY_MAX, penalty_weight * 1.5)
                elif ratio > 0.6:                       # vecina izvodljiva -> POPUSTI (trazi krace rute)
                    penalty_weight = max(PENALTY_MIN, penalty_weight / 1.3)
                if penalty_weight != old_w:
                    set_penalty_weight(penalty_weight)
                    # cene su sad u drugoj skali -> preracunaj referentne cene
                    current_cost = solution_cost(current)
                    best_cost = solution_cost(best_solution)

        # Restart ako smo zaglavili
        if no_improve_count >= max_no_improve:          # ako dugo nema napretka -> RESTART
            if verbose:
                print(f"  Restarting at iter {iteration} (no improvement for {no_improve_count} iters)")

            # Restart: pola puta SCTS ponovna podela (nova struktura ruta, izvodljiva po Tmax),
            # pola puta jak trzaj (ruin-recreate) -> raznovrsnost
            if random.random() < 0.5:
                current = scts_construction(instance, randomized=True)
            else:
                current = _perturb(best_solution, instance)
            current = fix_stations(current)
            _track_feasible(current)                     # mozda je restart vec izvodljiv
            current_cost = solution_cost(current)
            no_improve_count = 0

    # ---- IZBOR REZULTATA: vrati najbolje IZVODLJIVO ako postoji ----
    # (inace najbolje po ceni, koje moze biti neizvodljivo -> da pozivalac to vidi)
    result_solution = best_feasible if best_feasible is not None else best_solution

    elapsed = time.time() - start_time
    stats['total_time'] = elapsed

    eval_data = evaluate_solution(result_solution)      # finalna ocena vracenog resenja
    stats['final_distance'] = eval_data['total_distance']
    stats['final_feasible'] = eval_data['feasible']
    stats['final_n_routes'] = eval_data['n_routes']
    stats['feasible_found'] = best_feasible is not None
    stats['best_feasible_distance'] = best_feasible_dist if best_feasible is not None else None

    if verbose:
        print(f"\nGVNS completed:")
        print(f"  Total time: {elapsed:.2f}s")
        print(f"  Iterations: {stats['iterations']}")
        print(f"  Improvements: {stats['improvements']}")
        print(f"  Best distance: {stats['final_distance']:.2f}")
        print(f"  Feasible: {stats['final_feasible']}")
        print(f"  Routes: {stats['final_n_routes']}")
        print(f"  Time to best: {stats['time_to_best']:.2f}s")

    return result_solution, stats                       # vrati najbolje IZVODLJIVO (ili najbolje po ceni)


def _perturb(solution: Solution, instance: Instance) -> Solution:
    """Jaka perturbacija za restart: nasumicno izbaci pa ponovo ubaci musterije."""
    # ===== JAK TRZAJ (za restart): izbaci 30-50% musterija pa ih pohlepno vrati. =====
    sol = solution.copy()

    # Nasumicno izbaci 30-50% musterija
    all_customers = list(sol.get_all_customers())
    n_remove = max(1, int(len(all_customers) * random.uniform(0.3, 0.5)))  # koliko izbaciti (30-50%)
    to_remove = random.sample(all_customers, min(n_remove, len(all_customers)))

    for c in to_remove:                                 # izbaci izabrane musterije iz ruta
        for route in sol.routes:
            if c in route.nodes:
                route.nodes.remove(c)
                break

    sol.remove_empty_routes()

    # Pohlepno vrati izbacene musterije
    for c in to_remove:                                 # vrati ih na NAJJEFTINIJE mesto
        best_cost_increase = float('inf')
        best_route_idx = -1
        best_pos = -1

        for r_idx, route in enumerate(sol.routes):
            for pos in range(1, len(route.nodes)):
                prev = route.nodes[pos - 1]
                next_n = route.nodes[pos]
                cost_increase = (instance.dist(prev, c) + instance.dist(c, next_n)  # koliko poskupi ruta
                                - instance.dist(prev, next_n))

                if cost_increase < best_cost_increase:
                    best_cost_increase = cost_increase
                    best_route_idx = r_idx
                    best_pos = pos

        # Ili otvori novu rutu ako je jeftinije
        if len(sol.routes) < instance.n_vehicles:
            new_route_cost = instance.dist(0, c) + instance.dist(c, 0)
            if new_route_cost < best_cost_increase:
                sol.routes.append(Route([0, c, 0]))
                continue

        if best_route_idx >= 0:
            sol.routes[best_route_idx].insert(best_pos, c)  # ubaci na najbolje nadjeno mesto

    sol.invalidate_cache()
    return sol


def run_multiple(instance: Instance, n_runs: int = 5,
                 time_limit: float = 600.0, **kwargs) -> Tuple[Solution, dict]:
    """
    Pokreni GVNS vise puta i vrati najbolje resenje.

    Vraca: (best_solution, aggregated_stats)
    """
    # ===== Pokreni GVNS VISE puta (razliciti seed-ovi) i vrati najbolje + statistiku. =====
    best_solution = None
    best_cost = float('inf')
    all_stats = []

    per_run_time = time_limit / n_runs                  # vreme po jednom pokretanju

    for run in range(n_runs):
        print(f"\n{'='*60}")
        print(f"Run {run + 1}/{n_runs} (time limit: {per_run_time:.0f}s)")
        print(f"{'='*60}")

        sol, stats = gvns(instance, time_limit=per_run_time,
                         seed=run * 42 + 7, **kwargs)   # razlicit seed za svako pokretanje
        all_stats.append(stats)

        cost = solution_cost(sol)
        if cost < best_cost:                            # zadrzi NAJBOLJE od svih pokretanja
            best_cost = cost
            best_solution = sol

    # Objedini statistiku
    agg_stats = {                                       # objedinjena statistika (prosek, min, max, std...)
        'n_runs': n_runs,
        'best_distance': min(s['final_distance'] for s in all_stats),
        'avg_distance': np.mean([s['final_distance'] for s in all_stats]),
        'std_distance': np.std([s['final_distance'] for s in all_stats]),
        'worst_distance': max(s['final_distance'] for s in all_stats),
        'avg_time_to_best': np.mean([s['time_to_best'] for s in all_stats]),
        'avg_iterations': np.mean([s['iterations'] for s in all_stats]),
        'feasible_runs': sum(1 for s in all_stats if s['final_feasible']),
        'all_stats': all_stats,
    }

    print(f"\n{'='*60}")
    print(f"Aggregated Results ({n_runs} runs):")
    print(f"  Best distance: {agg_stats['best_distance']:.2f}")
    print(f"  Avg distance:  {agg_stats['avg_distance']:.2f} ± {agg_stats['std_distance']:.2f}")
    print(f"  Worst distance: {agg_stats['worst_distance']:.2f}")
    print(f"  Feasible runs: {agg_stats['feasible_runs']}/{n_runs}")
    print(f"{'='*60}")

    return best_solution, agg_stats
