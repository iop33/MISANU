// localsearch.hpp -- METS "Efficient Local Search" (ELS): granularna prva-poboljsanja
// pretraga nad METS funkcijom cilja. Radimo na rutama SAMO-musterija; AFS dodaje
// canonicalize (u cost_of), pa se kapacitet/domet vide kroz kaznu u ceni.
//
// Skup poteza (odgovara m1..m9 + Depot_* + NewRoute_*):
//   relocate (U posle/pre V), or-opt (segment 1-3), swap(U,V), 2-opt (unutar rute),
//   2-opt* (izmedju ruta), nova ruta za U. Prvo-poboljsanje, red kao u METS-u.
#pragma once
#include "instance.hpp"
#include "eval.hpp"
#include <vector>
#include <random>
#include <algorithm>

// SolC = rute samo-musterija (bez AFS-a). AFS se dodaje u cost_of preko canonicalize.
using SolC = std::vector<std::vector<int>>;

inline Solution to_solution(const SolC& s) { Solution o; for (auto& r : s) if (!r.empty()) o.push_back(r); return o; }

// Cena resenja: biraj AFS pozicije zajedno (vremenski svesno) pa oceni METS ciljem.
inline EvalResult eval_customers(const Instance& I, const SolC& s, const Weights& W) {
    SolC cust; for (auto& r : s) if (!r.empty()) cust.push_back(r);
    if (cust.empty()) return EvalResult{};
    return eval_afs_opt(I, cust, W).second;
}
// vrati i resenje sa najboljim AFS pozicijama
inline std::pair<Solution,EvalResult> eval_customers_full(const Instance& I, const SolC& s, const Weights& W) {
    SolC cust; for (auto& r : s) if (!r.empty()) cust.push_back(r);
    if (cust.empty()) return { {}, EvalResult{} };
    return eval_afs_opt(I, cust, W);
}
inline double cost_of(const Instance& I, const SolC& s, const Weights& W) { return eval_customers(I, s, W).cost; }

// 5 najblizih musterija za svaku musteriju (correlatedVertices).
inline std::vector<std::vector<int>> build_granular(const Instance& I, int K = 5) {
    std::vector<std::vector<int>> g(I.n_customers + 1);
    for (int c = 1; c <= I.n_customers; ++c) {
        std::vector<std::pair<double,int>> nb;
        for (int o = 1; o <= I.n_customers; ++o) if (o != c) nb.push_back({I.dist[c][o], o});
        std::sort(nb.begin(), nb.end());
        for (int k = 0; k < K && k < (int)nb.size(); ++k) g[c].push_back(nb[k].second);
    }
    return g;
}

// nadji (ruta, pozicija) musterije c
inline std::pair<int,int> locate(const SolC& s, int c) {
    for (int r = 0; r < (int)s.size(); ++r)
        for (int p = 0; p < (int)s[r].size(); ++p) if (s[r][p] == c) return {r, p};
    return {-1, -1};
}

inline void drop_empty(SolC& s) {
    s.erase(std::remove_if(s.begin(), s.end(), [](const std::vector<int>& r){ return r.empty(); }), s.end());
}

// ---- Generatori kandidata: vrate praznu SolC ako je potez nevalidan ----

// relocate: izbaci segment [pU..pU+len) iz rute rU i ubaci ga (mozda obrnut) posle V
inline SolC move_segment_after(const SolC& s, int rU, int pU, int len, int rV, int pV, bool reversed) {
    SolC c = s;
    if (pU + len > (int)c[rU].size()) return {};
    std::vector<int> seg(c[rU].begin() + pU, c[rU].begin() + pU + len);
    if (reversed) std::reverse(seg.begin(), seg.end());
    // proveri da V nije u segmentu
    // izbaci segment
    c[rU].erase(c[rU].begin() + pU, c[rU].begin() + pU + len);
    // preracunaj poziciju V (ista ruta, pomeraj)
    int insR = rV, insP = pV + 1;
    if (rV == rU && pV >= pU) insP = pV + 1 - len;
    if (insP < 0) insP = 0;
    if (insP > (int)c[insR].size()) insP = (int)c[insR].size();
    c[insR].insert(c[insR].begin() + insP, seg.begin(), seg.end());
    drop_empty(c);
    return c;
}

// swap dve musterije U i V
inline SolC swap_nodes(const SolC& s, int rU, int pU, int rV, int pV) {
    SolC c = s; std::swap(c[rU][pU], c[rV][pV]); return c;
}

// 2-opt unutar rute: obrni segment [i..j]
inline SolC two_opt_intra(const SolC& s, int r, int i, int j) {
    SolC c = s; std::reverse(c[r].begin() + i, c[r].begin() + j + 1); return c;
}

// 2-opt*: spoji [rU do pU] + [rV posle pV..] i [rV do pV] + [rU posle pU..]
inline SolC two_opt_star(const SolC& s, int rU, int pU, int rV, int pV) {
    if (rU == rV) return {};
    SolC c = s;
    std::vector<int> tailU(c[rU].begin() + pU + 1, c[rU].end());
    std::vector<int> tailV(c[rV].begin() + pV + 1, c[rV].end());
    c[rU].resize(pU + 1); c[rU].insert(c[rU].end(), tailV.begin(), tailV.end());
    c[rV].resize(pV + 1); c[rV].insert(c[rV].end(), tailU.begin(), tailU.end());
    drop_empty(c);
    return c;
}

// nova ruta za musteriju U
inline SolC new_route(const SolC& s, int rU, int pU) {
    SolC c = s; int u = c[rU][pU];
    c[rU].erase(c[rU].begin() + pU);
    c.push_back({u});
    drop_empty(c);
    return c;
}

// Pokusaj sve poteze za par (U,V); ako neki poboljsa -> primeni i vrati true (prvo-poboljsanje).
inline bool try_pair(const Instance& I, SolC& s, const Weights& W, double& cur, int U, int V, int max_routes) {
    auto lu = locate(s, U); auto lv = locate(s, V);
    if (lu.first < 0 || lv.first < 0) return false;
    int rU = lu.first, pU = lu.second, rV = lv.first, pV = lv.second;

    auto consider = [&](const SolC& cand) -> bool {
        if (cand.empty()) return false;
        double cc = cost_of(I, cand, W);
        if (cc < cur - 1e-6) { s = cand; cur = cc; return true; }
        return false;
    };

    // m1: relocate U posle V
    if (consider(move_segment_after(s, rU, pU, 1, rV, pV, false))) return true;
    // Depot_m1 analog: relocate U pre V (na pocetak, tj. posle pV-1)
    if (pV > 0 && consider(move_segment_after(s, rU, pU, 1, rV, pV - 1, false))) return true;
    // m2: relocate par (U,X) posle V
    if (pU + 1 < (int)s[rU].size() && consider(move_segment_after(s, rU, pU, 2, rV, pV, false))) return true;
    // m3: relocate par obrnut (X,U) posle V
    if (pU + 1 < (int)s[rU].size() && consider(move_segment_after(s, rU, pU, 2, rV, pV, true))) return true;
    // or-opt segment 3
    if (pU + 2 < (int)s[rU].size() && consider(move_segment_after(s, rU, pU, 3, rV, pV, false))) return true;
    // m4: swap U i V
    if (!(rU == rV && pU == pV) && consider(swap_nodes(s, rU, pU, rV, pV))) return true;
    // m7: 2-opt unutar rute (ako su U,V u istoj ruti)
    if (rU == rV) { int i = std::min(pU, pV), j = std::max(pU, pV); if (j > i && consider(two_opt_intra(s, rU, i, j))) return true; }
    // m8/m9: 2-opt* izmedju ruta
    if (rU != rV) { if (consider(two_opt_star(s, rU, pU, rV, pV))) return true; }
    // NewRoute za U
    if ((int)s.size() < max_routes && s[rU].size() > 1 && consider(new_route(s, rU, pU))) return true;

    return false;
}

// Glavna ELS petlja: granularno, prvo-poboljsanje, prolazi dok ima poboljsanja.
inline void local_search(const Instance& I, SolC& s, const Weights& W,
                         const std::vector<std::vector<int>>& gran, std::mt19937& rng) {
    double cur = cost_of(I, s, W);
    int max_routes = std::max((int)s.size() + 2, I.n_vehicles);
    bool improved = true; int guard = 0;
    while (improved && guard++ < 100000) {
        improved = false;
        for (int U = 1; U <= I.n_customers; ++U) {
            for (int V : gran[U]) {
                if (U == V) continue;
                if (try_pair(I, s, W, cur, U, V, max_routes)) { improved = true; break; }
            }
        }
    }
    drop_empty(s);
}
