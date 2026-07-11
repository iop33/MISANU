// solution.hpp -- Ruta, Resenje, ocena + Reschedule (red na pumpi) + kazne.
// Verni port solution.py: ista funkcija cilja i isti uslovi izvodljivosti kao METS.
#pragma once
#include "instance.hpp"
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

using Route = std::vector<int>;                 // niz cvorova [0, ..., 0]
using Solution = std::vector<Route>;            // skup ruta

// ---- Globalni (adaptivni) kazneni faktor -- deli ga cela pretraga ----
inline double PENALTY_WEIGHT = 1000.0;

struct RouteEval {
    double distance = 0, duration = 0;
    bool fuel_feasible = true;
    std::vector<double> arrival, fuel;
};

inline RouteEval evaluate_route(const Route& r, const Instance& ins) {
    RouteEval e;
    int n = (int)r.size();
    if (n <= 2) { e.arrival = {0,0}; e.fuel = {ins.tank_capacity, ins.tank_capacity}; return e; }
    e.arrival.assign(n, 0.0); e.fuel.assign(n, 0.0);
    std::vector<double> dep(n, 0.0);
    e.fuel[0] = ins.tank_capacity;
    dep[0] = ins.p_start;
    for (int i = 1; i < n; ++i) {
        int prev = r[i-1], curr = r[i];
        double d = ins.dist[prev][curr], t = ins.tt[prev][curr];
        e.distance += d;
        e.arrival[i] = dep[i-1] + t;
        e.fuel[i] = e.fuel[i-1] - ins.consumption_rate * d;
        if (e.fuel[i] < -1e-6) e.fuel_feasible = false;
        if (ins.is_customer(curr))      dep[i] = e.arrival[i] + ins.service_time;
        else if (ins.is_station(curr)) { dep[i] = e.arrival[i] + ins.refuel_time; e.fuel[i] = ins.tank_capacity; }
        else                            dep[i] = e.arrival[i];
    }
    e.duration = e.arrival[n-1];
    return e;
}

struct Reschedule {
    bool feasible = true;
    double total_wait = 0.0;
    std::vector<double> route_wait;             // ukupno cekanje po ruti
};

// Red na pumpi (Xu Alg. A1): odlozi vozila koja se preklapaju iznad kapaciteta.
inline Reschedule reschedule(const Instance& ins, const Solution& sol) {
    Reschedule rz;
    rz.route_wait.assign(sol.size(), 0.0);
    std::vector<RouteEval> evals; evals.reserve(sol.size());
    for (auto& r : sol) evals.push_back(evaluate_route(r, ins));

    // grupisi posete po punionici: (route_idx, pos, arrival)
    struct Visit { int r, pos; double arr; };
    std::vector<std::vector<Visit>> per_station(ins.n_total());
    for (int ri = 0; ri < (int)sol.size(); ++ri)
        for (int pos = 0; pos < (int)sol[ri].size(); ++pos) {
            int node = sol[ri][pos];
            if (ins.is_station(node)) per_station[node].push_back({ri, pos, evals[ri].arrival[pos]});
        }

    double refuel = ins.refuel_time;
    for (int st = 0; st < ins.n_total(); ++st) {
        auto& visits = per_station[st];
        if (visits.empty()) continue;
        int cap = ins.station_cap(st);
        if ((int)visits.size() <= cap) continue;
        std::sort(visits.begin(), visits.end(), [](const Visit&a, const Visit&b){ return a.arr < b.arr; });
        int m = (int)visits.size();
        std::vector<double> adj(m);
        for (int k = 0; k < m; ++k) adj[k] = visits[k].arr;

        bool changed = true; int it = 0;
        while (changed && it < 1000) {
            changed = false; ++it;
            for (int i = 0; i < m-1 && !changed; ++i) {
                int overlap = 0;
                double si = adj[i], ei = si + refuel;
                for (int j = i+1; j < m; ++j) {
                    double sj = adj[j], ej = sj + refuel;
                    if (si < ej && sj < ei) {
                        ++overlap;
                        if (overlap > cap - 1) {
                            if (ej <= ei) adj[i] = ej; else adj[j] = ei;
                            changed = true; break;
                        }
                    }
                }
            }
        }
        for (int k = 0; k < m; ++k) {
            double wt = adj[k] - visits[k].arr;
            if (wt > 1e-6) { rz.route_wait[visits[k].r] += wt; rz.total_wait += wt; }
        }
    }
    // izvodljivost: trajanje + cekanje te rute <= Tmax
    for (int ri = 0; ri < (int)sol.size(); ++ri)
        if (evals[ri].duration + rz.route_wait[ri] > ins.t_max + 1e-6) rz.feasible = false;
    return rz;
}

struct SolEval {
    double total_distance = 0;
    bool all_served = false, vehicle_ok = false, fuel_ok = true, duration_ok = true, capacity_ok = true;
    bool feasible = false;
    int n_routes = 0;
};

inline double total_distance(const Instance& ins, const Solution& sol) {
    double td = 0;
    for (auto& r : sol) for (size_t i = 0; i+1 < r.size(); ++i) td += ins.dist[r[i]][r[i+1]];
    return td;
}

inline SolEval evaluate_solution(const Instance& ins, const Solution& sol) {
    SolEval s; s.n_routes = (int)sol.size();
    std::vector<RouteEval> evals; evals.reserve(sol.size());
    for (auto& r : sol) {
        RouteEval e = evaluate_route(r, ins);
        evals.push_back(e);
        s.total_distance += e.distance;
        if (!e.fuel_feasible) s.fuel_ok = false;
        if (e.duration > ins.t_max + 1e-6) s.duration_ok = false;
    }
    // svaka musterija tacno jednom
    std::vector<int> cnt(ins.n_total(), 0);
    for (auto& r : sol) for (int node : r) if (ins.is_customer(node)) cnt[node]++;
    int served = 0; bool dup = false;
    for (int c = 1; c <= ins.n_customers; ++c) { if (cnt[c] > 0) served++; if (cnt[c] > 1) dup = true; }
    s.all_served = (served == ins.n_customers) && !dup;
    s.vehicle_ok = (int)sol.size() <= ins.n_vehicles;

    Reschedule rz = reschedule(ins, sol);
    s.capacity_ok = rz.feasible;
    for (int ri = 0; ri < (int)sol.size(); ++ri)
        if (evals[ri].duration + rz.route_wait[ri] > ins.t_max + 1e-6) s.duration_ok = false;

    s.feasible = s.all_served && s.vehicle_ok && s.fuel_ok && s.duration_ok && s.capacity_ok;
    return s;
}

// Kazna za prekrsaje (mirror compute_penalty).
inline double compute_penalty(const Instance& ins, const Solution& sol) {
    double penalty = 0;
    for (auto& r : sol) {
        RouteEval e = evaluate_route(r, ins);
        if (e.duration > ins.t_max) penalty += (e.duration - ins.t_max);
        for (double fl : e.fuel) if (fl < -1e-6) penalty += std::fabs(fl);
    }
    Reschedule rz = reschedule(ins, sol);
    if (!rz.feasible) penalty += rz.total_wait + 1.0;
    return penalty;
}

inline double penalized_cost(const Instance& ins, const Solution& sol) {
    return total_distance(ins, sol) + PENALTY_WEIGHT * compute_penalty(ins, sol);
}
