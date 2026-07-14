// gvns.hpp -- General Variable Neighborhood Search (C++ port).
// Verni port gvns.py: 3 pocetna resenja, VND, adaptivna kazna, restart,
// vraca NAJBOLJE IZVODLJIVO resenje.
#pragma once
#include "solution.hpp"
#include "construction.hpp"
#include "neighborhoods.hpp"
#include <random>
#include <chrono>
#include <functional>

struct GvnsResult {
    Solution solution;
    double distance = 1e18;
    bool feasible = false;
    int n_routes = 0;
    double time_to_best = 0;
};

// VND lokalna pretraga. 'deadline': TVRDI vremenski prekid -- bez njega bi na
// velikim instancama (50-100 musterija) jedan VND mogao da traje minutima i
// visestruko probije --time-limit (uoceno kao visesatno izvrsavanje).
inline Solution vnd(const Instance& ins, Solution sol,
                    std::chrono::steady_clock::time_point deadline
                        = std::chrono::steady_clock::time_point::max()) {
    double cost = penalized_cost(ins, sol);
    int k = 0;
    while (k < 3) {
        if (std::chrono::steady_clock::now() >= deadline) break;   // postuj rok
        Solution nb;
        if (k == 0) nb = relocate_best(ins, sol);
        else if (k == 1) nb = swap_best(ins, sol);
        else nb = two_opt_intra_best(ins, sol);
        double nc = penalized_cost(ins, nb);
        if (nc < cost - 1e-6) { sol = nb; cost = nc; k = 0; }
        else ++k;
    }
    return sol;
}

// Jaka perturbacija: izbaci ~40% musterija i vrati ih pohlepno.
inline Solution perturb(const Instance& ins, const Solution& best, std::mt19937& rng) {
    Solution s = best;
    std::vector<int> custs;
    for (auto& r : s) for (int x : r) if (ins.is_customer(x)) custs.push_back(x);
    std::shuffle(custs.begin(), custs.end(), rng);
    int nrem = std::max(1, (int)(custs.size() * (0.3 + 0.2 * (rng() % 100) / 100.0)));
    std::vector<int> rem(custs.begin(), custs.begin() + std::min(nrem, (int)custs.size()));
    for (int c : rem) for (auto& r : s) { auto it = std::find(r.begin(), r.end(), c); if (it != r.end()) { r.erase(it); break; } }
    remove_empty(s);
    for (int c : rem) {
        double best_inc = 1e18; int br = -1, bp = -1;
        for (int r = 0; r < (int)s.size(); ++r)
            for (int pos = 1; pos < (int)s[r].size(); ++pos) {
                double inc = ins.dist[s[r][pos-1]][c] + ins.dist[c][s[r][pos]] - ins.dist[s[r][pos-1]][s[r][pos]];
                if (inc < best_inc) { best_inc = inc; br = r; bp = pos; }
            }
        if ((int)s.size() < ins.n_vehicles && ins.dist[0][c] + ins.dist[c][0] < best_inc) s.push_back({0, c, 0});
        else if (br >= 0) s[br].insert(s[br].begin() + bp, c);
    }
    return s;
}

inline GvnsResult gvns(const Instance& ins, double time_limit, int max_no_improve, int k_max, unsigned seed) {
    using clock = std::chrono::steady_clock;
    auto t0 = clock::now();
    auto elapsed = [&]{ return std::chrono::duration<double>(clock::now() - t0).count(); };
    // TVRDI rok za ceo run -- postuje ga i VND (prekida se usred pretrage).
    auto deadline = t0 + std::chrono::duration_cast<clock::duration>(
                             std::chrono::duration<double>(time_limit));
    std::mt19937 rng(seed);

    PENALTY_WEIGHT = 1000.0;                         // reset po run-u
    const double PMIN = 1000.0, PMAX = 1e5;
    const int ADAPT = 30;

    Solution best_feasible; double best_feasible_dist = 1e18; bool have_feasible = false;
    double time_to_best = 0;
    auto track = [&](const Solution& s){
        SolEval e = evaluate_solution(ins, s);
        if (e.feasible && e.total_distance < best_feasible_dist - 1e-9) {
            best_feasible = s; best_feasible_dist = e.total_distance; have_feasible = true;
            time_to_best = elapsed();
        }
        return e.feasible;
    };

    // Tri pocetna resenja
    std::vector<Solution> cand = {
        fix_stations(greedy_construction(ins), ins),
        fix_stations(savings_construction(ins), ins),
        scts_construction(ins, rng, false),
    };
    Solution best = cand[0]; double best_cost = penalized_cost(ins, best);
    for (auto& c : cand) { double cc = penalized_cost(ins, c); if (cc < best_cost) { best_cost = cc; best = c; } track(c); }
    Solution current = best; double current_cost = best_cost;

    // shaking okoline
    using Shake = std::function<Solution(const Instance&, const Solution&, std::mt19937&)>;
    std::vector<Shake> shakes = { relocate_shake, swap_shake, or_opt_shake, two_opt_intra_shake, two_opt_star_shake };
    int K = std::min(k_max, (int)shakes.size());

    std::vector<int> feas_window;
    int no_improve = 0; long iteration = 0;
    double penalty_weight = PENALTY_WEIGHT;

    while (elapsed() < time_limit) {
        int k = 0;
        while (k < K) {
            if (elapsed() >= time_limit) break;
            ++iteration;
            Solution neighbor = shakes[k](ins, current, rng);
            neighbor = fix_stations(neighbor, ins);
            neighbor = vnd(ins, neighbor, deadline);
            double nc = penalized_cost(ins, neighbor);

            feas_window.push_back(track(neighbor) ? 1 : 0);

            if (nc < current_cost - 1e-6) {
                current = neighbor; current_cost = nc; k = 0; no_improve = 0;
                if (nc < best_cost - 1e-6) { best = neighbor; best_cost = nc; }
            } else { ++k; ++no_improve; }

            if (iteration % ADAPT == 0 && (int)feas_window.size() >= ADAPT) {
                int cntf = 0; for (int t = (int)feas_window.size() - ADAPT; t < (int)feas_window.size(); ++t) cntf += feas_window[t];
                double ratio = (double)cntf / ADAPT, old = penalty_weight;
                if (ratio < 0.15) penalty_weight = std::min(PMAX, penalty_weight * 1.5);
                else if (ratio > 0.6) penalty_weight = std::max(PMIN, penalty_weight / 1.3);
                if (penalty_weight != old) {
                    PENALTY_WEIGHT = penalty_weight;
                    current_cost = penalized_cost(ins, current);
                    best_cost = penalized_cost(ins, best);
                }
            }
        }
        if (no_improve >= max_no_improve) {
            if ((rng() % 100) < 50) current = scts_construction(ins, rng, true);
            else current = perturb(ins, best, rng);
            current = fix_stations(current, ins);
            track(current);
            current_cost = penalized_cost(ins, current);
            no_improve = 0;
        }
    }

    GvnsResult res;
    Solution ret = have_feasible ? best_feasible : best;
    SolEval e = evaluate_solution(ins, ret);
    res.solution = ret; res.distance = e.total_distance; res.feasible = e.feasible;
    res.n_routes = e.n_routes; res.time_to_best = time_to_best;
    return res;
}
