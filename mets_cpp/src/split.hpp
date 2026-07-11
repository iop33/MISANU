// split.hpp -- SCTS: podela "giant tour"-a (permutacije musterija) na rute.
// split_Tmax: rute ogranicene trajanjem (bez AFS-a).
// split_Dmax: rute ogranicene dometom (AFS se ubacuje kad domet zahteva).
// (Verno duhu METS-a; AFS pozicija bira se preko canonicalize_route iz eval.hpp.)
#pragma once
#include "instance.hpp"
#include "eval.hpp"
#include <vector>

// Trajanje rute [depo, custs, depo] bez AFS-a (kao T_able).
inline double route_time_noafs(const Instance& I, const std::vector<int>& custs) {
    if (custs.empty()) return 0.0;
    double d = I.dist[0][custs[0]];
    for (int i = 0; i + 1 < (int)custs.size(); ++i) d += I.dist[custs[i]][custs[i+1]];
    d += I.dist[custs.back()][0];
    return d / I.speed + (int)custs.size() * I.service_time;
}

// Da li skup musterija (u datom redu) moze u JEDNU rutu <= Dmax uz najvise jedan AFS?
inline bool fits_dmax_with_afs(const Instance& I, const std::vector<int>& custs) {
    if (custs.empty()) return true;
    // bez AFS
    double full = I.dist[0][custs[0]];
    for (int i = 0; i + 1 < (int)custs.size(); ++i) full += I.dist[custs[i]][custs[i+1]];
    full += I.dist[custs.back()][0];
    if (full <= I.d_max + 1e-9) return true;
    // sa jednim AFS-om: probaj svaku poziciju, treba obe deonice <= Dmax
    int afs = I.first_station();
    for (int p = 0; p < (int)custs.size(); ++p) {
        Route t(custs.begin(), custs.begin() + p + 1);
        t.push_back(afs);
        t.insert(t.end(), custs.begin() + p + 1, custs.end());
        LegInfo L = get_pd_pt(I, t);
        if (L.d1 <= I.d_max + 1e-9 && L.d2 <= I.d_max + 1e-9) return true;
    }
    return false;
}

// splitTmax: pohlepno najduzi prefiks cije trajanje <= Tmax; rute bez AFS-a.
inline Solution split_Tmax(const Instance& I, const std::vector<int>& tsp) {
    Solution sol;
    std::vector<int> cur;
    for (int c : tsp) {
        std::vector<int> trial = cur; trial.push_back(c);
        if (route_time_noafs(I, trial) <= I.t_max + 1e-9) cur = trial;
        else {
            if (!cur.empty()) sol.push_back(cur);
            cur = {c};                              // nova ruta krece od c
        }
    }
    if (!cur.empty()) sol.push_back(cur);
    return sol;                                     // (AFS dodaje lokalna pretraga)
}

// splitDmax: pohlepno rasti rutu dok staje u Dmax (uz najvise jedan AFS); AFS ubaci na kraju.
inline Solution split_Dmax(const Instance& I, const std::vector<int>& tsp) {
    Solution sol;
    std::vector<int> cur;
    for (int c : tsp) {
        std::vector<int> trial = cur; trial.push_back(c);
        if (fits_dmax_with_afs(I, trial)) cur = trial;
        else { if (!cur.empty()) sol.push_back(cur); cur = {c}; }
    }
    if (!cur.empty()) sol.push_back(cur);
    // ubaci AFS gde domet zahteva (canonicalize_route bira poziciju)
    for (auto& r : sol) canonicalize_route(I, r);
    return sol;
}
