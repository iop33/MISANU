// neighborhoods.hpp -- Potezi: lokalna pretraga (best) + shaking (random).
#pragma once
#include "solution.hpp"
#include "construction.hpp"
#include <random>

inline void remove_empty(Solution& sol) {
    sol.erase(std::remove_if(sol.begin(), sol.end(),
              [](const Route& r){ return r.size() <= 2; }), sol.end());
}

// pozicije svih musterija: (ruta, indeks)
inline std::vector<std::pair<int,int>> customer_positions(const Instance& ins, const Solution& sol) {
    std::vector<std::pair<int,int>> p;
    for (int r = 0; r < (int)sol.size(); ++r)
        for (int i = 1; i + 1 < (int)sol[r].size(); ++i)
            if (ins.is_customer(sol[r][i])) p.push_back({r, i});
    return p;
}

// ---------- LOKALNA PRETRAGA (best-improvement) ----------
inline Solution relocate_best(const Instance& ins, const Solution& sol) {
    Solution best = sol; double best_cost = penalized_cost(ins, sol);
    for (int r = 0; r < (int)sol.size(); ++r)
        for (int pos = 1; pos + 1 < (int)sol[r].size(); ++pos) {
            if (!ins.is_customer(sol[r][pos])) continue;
            int cust = sol[r][pos];
            for (int tr = 0; tr < (int)sol.size(); ++tr)
                for (int tp = 1; tp < (int)sol[tr].size(); ++tp) {
                    if (tr == r && (tp == pos || tp == pos + 1)) continue;
                    Solution cand = sol;
                    cand[r].erase(cand[r].begin() + pos);
                    int atp = tp; if (tr == r && tp > pos) atp -= 1;
                    cand[tr].insert(cand[tr].begin() + atp, cust);
                    remove_empty(cand);
                    double c = penalized_cost(ins, cand);
                    if (c < best_cost - 1e-6) { best_cost = c; best = cand; }
                }
        }
    return best;
}

inline Solution swap_best(const Instance& ins, const Solution& sol) {
    Solution best = sol; double best_cost = penalized_cost(ins, sol);
    auto pos = customer_positions(ins, sol);
    for (size_t i = 0; i < pos.size(); ++i)
        for (size_t j = i + 1; j < pos.size(); ++j) {
            Solution cand = sol;
            std::swap(cand[pos[i].first][pos[i].second], cand[pos[j].first][pos[j].second]);
            double c = penalized_cost(ins, cand);
            if (c < best_cost - 1e-6) { best_cost = c; best = cand; }
        }
    return best;
}

inline Solution two_opt_intra_best(const Instance& ins, const Solution& sol) {
    Solution best = sol; double best_cost = penalized_cost(ins, sol);
    for (int r = 0; r < (int)sol.size(); ++r) {
        int n = (int)sol[r].size();
        for (int i = 1; i < n - 2; ++i)
            for (int j = i + 1; j < n - 1; ++j) {
                Solution cand = sol;
                std::reverse(cand[r].begin() + i, cand[r].begin() + j + 1);
                double c = penalized_cost(ins, cand);
                if (c < best_cost - 1e-6) { best_cost = c; best = cand; }
            }
    }
    return best;
}

// ---------- SHAKING (random) ----------
inline Solution relocate_shake(const Instance& ins, const Solution& sol, std::mt19937& rng) {
    Solution s = sol;
    auto pos = customer_positions(ins, s);
    if (pos.empty()) return s;
    auto p = pos[rng() % pos.size()];
    int cust = s[p.first][p.second];
    s[p.first].erase(s[p.first].begin() + p.second);
    int nr = (int)s.size();
    int tr = (int)(rng() % (nr + ((int)s.size() < ins.n_vehicles ? 1 : 0)));
    if (tr >= nr) { s.push_back({0, cust, 0}); }
    else { int ip = 1 + rng() % std::max(1, (int)s[tr].size() - 1); s[tr].insert(s[tr].begin() + ip, cust); }
    remove_empty(s);
    return s;
}

inline Solution swap_shake(const Instance& ins, const Solution& sol, std::mt19937& rng) {
    Solution s = sol;
    auto pos = customer_positions(ins, s);
    if (pos.size() < 2) return s;
    int a = rng() % pos.size(), b = rng() % pos.size();
    if (a == b) return s;
    std::swap(s[pos[a].first][pos[a].second], s[pos[b].first][pos[b].second]);
    return s;
}

inline Solution or_opt_shake(const Instance& ins, const Solution& sol, std::mt19937& rng) {
    Solution s = sol;
    std::vector<int> with; for (int r = 0; r < (int)s.size(); ++r) { for (int x : s[r]) if (ins.is_customer(x)) { with.push_back(r); break; } }
    if (with.empty()) return s;
    int r = with[rng() % with.size()];
    std::vector<int> cp; for (int i = 1; i + 1 < (int)s[r].size(); ++i) if (ins.is_customer(s[r][i])) cp.push_back(i);
    if (cp.empty()) return s;
    int seg = 1 + rng() % 3; seg = std::min(seg, (int)cp.size());
    int start = rng() % (cp.size() - seg + 1);
    std::vector<int> seg_nodes; for (int k = 0; k < seg; ++k) seg_nodes.push_back(s[r][cp[start + k]]);
    for (int k = seg - 1; k >= 0; --k) s[r].erase(s[r].begin() + cp[start + k]);
    int tr = rng() % s.size();
    int ip = 1 + rng() % std::max(1, (int)s[tr].size() - 1);
    for (int k = 0; k < seg; ++k) s[tr].insert(s[tr].begin() + ip + k, seg_nodes[k]);
    remove_empty(s);
    return s;
}

inline Solution two_opt_intra_shake(const Instance& ins, const Solution& sol, std::mt19937& rng) {
    Solution s = sol;
    std::vector<int> ok; for (int r = 0; r < (int)s.size(); ++r) if (s[r].size() > 3) ok.push_back(r);
    if (ok.empty()) return s;
    int r = ok[rng() % ok.size()]; int n = (int)s[r].size();
    int i = 1 + rng() % (n - 3); int j = i + 1 + rng() % (n - 2 - i);
    std::reverse(s[r].begin() + i, s[r].begin() + j + 1);
    return s;
}

inline Solution two_opt_star_shake(const Instance& ins, const Solution& sol, std::mt19937& rng) {
    Solution s = sol;
    if (s.size() < 2) return s;
    int a = rng() % s.size(), b = rng() % s.size();
    if (a == b || s[a].size() <= 2 || s[b].size() <= 2) return s;
    int ca = 1 + rng() % (s[a].size() - 2), cb = 1 + rng() % (s[b].size() - 2);
    Route ta(s[a].begin() + ca, s[a].end()), tb(s[b].begin() + cb, s[b].end());
    s[a].resize(ca); s[a].insert(s[a].end(), tb.begin(), tb.end());
    s[b].resize(cb); s[b].insert(s[b].end(), ta.begin(), ta.end());
    remove_empty(s);
    return s;
}
