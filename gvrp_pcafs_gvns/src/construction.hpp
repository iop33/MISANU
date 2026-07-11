// construction.hpp -- Pocetna resenja: greedy, savings (Clarke&Wright), SCTS.
#pragma once
#include "solution.hpp"
#include <random>
#include <set>

inline int find_nearest_station(const Instance& ins, int from, double fuel) {
    int best = -1; double bestd = 1e18;
    for (int s = ins.first_station(); s < ins.n_total(); ++s) {
        double d = ins.dist[from][s];
        if (ins.consumption_rate * d <= fuel + 1e-6 && d < bestd) { bestd = d; best = s; }
    }
    return best;
}

// Ubaci punionicu gde god bi gorivo ponestalo.
inline Route insert_station_if_needed(const Route& nodes, const Instance& ins) {
    Route out; out.push_back(nodes[0]);
    double fuel = ins.tank_capacity;
    for (size_t i = 1; i < nodes.size(); ++i) {
        int prev = out.back(), curr = nodes[i];
        double need = ins.consumption_rate * ins.dist[prev][curr];
        if (need > fuel + 1e-6) {
            int s = find_nearest_station(ins, prev, fuel);
            if (s != -1) { out.push_back(s); fuel = ins.tank_capacity; need = ins.consumption_rate * ins.dist[s][curr]; }
        }
        out.push_back(curr);
        if (ins.is_station(curr)) fuel = ins.tank_capacity; else fuel -= need;
    }
    return out;
}

// Izbaci sve punionice pa ih vrati samo gde treba (u svakoj ruti).
inline Solution fix_stations(const Solution& sol, const Instance& ins) {
    Solution out;
    for (auto& r : sol) {
        Route cust; cust.push_back(0);
        for (size_t i = 1; i + 1 < r.size(); ++i) if (ins.is_customer(r[i])) cust.push_back(r[i]);
        cust.push_back(0);
        out.push_back(insert_station_if_needed(cust, ins));
    }
    return out;
}

inline bool route_fuel_ok(const Route& r, const Instance& ins) {
    double fuel = ins.tank_capacity;
    for (size_t i = 1; i < r.size(); ++i) {
        fuel -= ins.consumption_rate * ins.dist[r[i-1]][r[i]];
        if (fuel < -1e-6) return false;
        if (ins.is_station(r[i])) fuel = ins.tank_capacity;
    }
    return true;
}

// --- Pohlepna konstrukcija (najblizi sused, uz proveru povratka) ---
inline Solution greedy_construction(const Instance& ins) {
    Solution sol;
    std::set<int> unserved;
    for (int c = 1; c <= ins.n_customers; ++c) unserved.insert(c);

    while (!unserved.empty() && (int)sol.size() < ins.n_vehicles) {
        Route route = {0};
        int cur = 0; double fuel = ins.tank_capacity, t = ins.p_start;
        while (!unserved.empty()) {
            int best = -1; double bestd = 1e18;
            for (int c : unserved) {
                double d = ins.dist[cur][c], need = ins.consumption_rate * d;
                if (need > fuel + 1e-6) continue;                  // ne mozemo do c
                double fuel_after = fuel - need;
                double t_after = t + ins.tt[cur][c] + ins.service_time;
                bool can_close = false;
                if (fuel_after >= ins.consumption_rate * ins.dist[c][0] - 1e-6 &&
                    t_after + ins.tt[c][0] <= ins.t_max + 1e-6) can_close = true;
                else for (int s = ins.first_station(); s < ins.n_total(); ++s) {
                    if (fuel_after >= ins.consumption_rate * ins.dist[c][s] - 1e-6) {
                        double tv = t_after + ins.tt[c][s] + ins.refuel_time + ins.tt[s][0];
                        if (tv <= ins.t_max + 1e-6 &&
                            ins.tank_capacity >= ins.consumption_rate * ins.dist[s][0] - 1e-6) { can_close = true; break; }
                    }
                }
                if (can_close && d < bestd) { bestd = d; best = c; }
            }
            if (best == -1) break;
            double need = ins.consumption_rate * ins.dist[cur][best];
            if (need > fuel + 1e-6) {                              // dopuna pre musterije
                int bs = -1; double bc = 1e18;
                for (int s = ins.first_station(); s < ins.n_total(); ++s)
                    if (ins.consumption_rate*ins.dist[cur][s] <= fuel+1e-6 &&
                        ins.consumption_rate*ins.dist[s][best] <= ins.tank_capacity+1e-6) {
                        double c2 = ins.dist[cur][s] + ins.dist[s][best];
                        if (c2 < bc) { bc = c2; bs = s; }
                    }
                if (bs == -1) break;
                route.push_back(bs); t += ins.tt[cur][bs] + ins.refuel_time; fuel = ins.tank_capacity; cur = bs;
                need = ins.consumption_rate * ins.dist[cur][best];
                if (need > fuel + 1e-6) break;
            }
            route.push_back(best); fuel -= need;
            t += ins.tt[cur][best] + ins.service_time; cur = best; unserved.erase(best);
            if (fuel < ins.consumption_rate * ins.dist[cur][0] - 1e-6) {   // dopuna za povratak
                int s = find_nearest_station(ins, cur, fuel);
                if (s != -1) {
                    double tc = t + ins.tt[cur][s] + ins.refuel_time + ins.tt[s][0];
                    if (tc <= ins.t_max + 1e-6) { route.push_back(s); fuel = ins.tank_capacity; t += ins.tt[cur][s] + ins.refuel_time; cur = s; }
                }
            }
        }
        route.push_back(0);
        int nc = 0; for (int x : route) if (ins.is_customer(x)) nc++;
        if (nc > 0) sol.push_back(route); else break;
    }
    // rezerva: preostale musterije kao pojedinacne rute
    for (int c : unserved) {
        if ((int)sol.size() >= ins.n_vehicles) break;
        sol.push_back(insert_station_if_needed({0, c, 0}, ins));
    }
    return sol;
}

// --- Savings (Clarke & Wright) ---
inline Solution savings_construction(const Instance& ins) {
    std::vector<Route> routes(ins.n_customers + 1);
    std::vector<int> where(ins.n_customers + 1);
    for (int c = 1; c <= ins.n_customers; ++c) { routes[c] = {0, c, 0}; where[c] = c; }

    struct Sv { double s; int i, j; };
    std::vector<Sv> sav;
    for (int i = 1; i <= ins.n_customers; ++i)
        for (int j = 1; j <= ins.n_customers; ++j) if (i != j)
            sav.push_back({ins.dist[i][0] + ins.dist[0][j] - ins.dist[i][j], i, j});
    std::sort(sav.begin(), sav.end(), [](const Sv&a, const Sv&b){ return a.s > b.s; });

    for (auto& e : sav) {
        if (e.s <= 0) break;
        int ri = where[e.i], rj = where[e.j];
        if (ri == rj || routes[ri].empty() || routes[rj].empty()) continue;
        // i mora biti poslednja musterija svoje rute, j prva u svojoj
        auto lastc = [&](const Route& r){ for (int k=(int)r.size()-1;k>=0;--k) if (ins.is_customer(r[k])) return r[k]; return -1; };
        auto firstc = [&](const Route& r){ for (int x : r) if (ins.is_customer(x)) return x; return -1; };
        if (lastc(routes[ri]) != e.i || firstc(routes[rj]) != e.j) continue;
        Route merged(routes[ri].begin(), routes[ri].end()-1);
        merged.insert(merged.end(), routes[rj].begin()+1, routes[rj].end());
        RouteEval ev = evaluate_route(merged, ins);
        if (ev.fuel_feasible && ev.duration <= ins.t_max + 1e-6) {
            for (int x : routes[rj]) if (ins.is_customer(x)) where[x] = ri;
            routes[ri] = merged; routes[rj].clear();
        }
    }
    Solution sol;
    for (int c = 1; c <= ins.n_customers; ++c) if (!routes[c].empty()) sol.push_back(routes[c]);
    return sol;
}

// --- SCTS: giant tour (najblizi sused) + podela po Tmax ---
inline std::vector<int> nn_order(const Instance& ins, std::mt19937& rng, bool randomized) {
    std::set<int> left; for (int c = 1; c <= ins.n_customers; ++c) left.insert(c);
    std::vector<int> order;
    int cur;
    if (randomized) { auto it = left.begin(); std::advance(it, rng() % left.size()); cur = *it; }
    else { cur = 1; double bd = 1e18; for (int c : left) if (ins.dist[0][c] < bd) { bd = ins.dist[0][c]; cur = c; } }
    order.push_back(cur); left.erase(cur);
    std::uniform_real_distribution<double> U(0,1);
    while (!left.empty()) {
        int nxt = *left.begin();
        if (randomized && U(rng) < 0.3) { auto it = left.begin(); std::advance(it, rng() % left.size()); nxt = *it; }
        else { double bd = 1e18; for (int c : left) if (ins.dist[cur][c] < bd) { bd = ins.dist[cur][c]; nxt = c; } }
        order.push_back(nxt); left.erase(nxt); cur = nxt;
    }
    return order;
}

inline Solution scts_construction(const Instance& ins, std::mt19937& rng, bool randomized) {
    std::vector<int> order = nn_order(ins, rng, randomized);
    Solution sol;
    std::vector<int> cur;
    auto make = [&](const std::vector<int>& cs){ Route r; r.push_back(0); for (int c : cs) r.push_back(c); r.push_back(0); return insert_station_if_needed(r, ins); };
    for (int c : order) {
        std::vector<int> trial = cur; trial.push_back(c);
        Route r = make(trial);
        RouteEval ev = evaluate_route(r, ins);
        if (ev.fuel_feasible && ev.duration <= ins.t_max + 1e-6) cur = trial;
        else { if (!cur.empty()) sol.push_back(make(cur)); cur = {c}; }
    }
    if (!cur.empty()) sol.push_back(make(cur));
    return sol;
}
