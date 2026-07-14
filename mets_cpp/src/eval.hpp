// eval.hpp -- Jezgro METS evaluacije (verni prevod chromR_detail_all + get_pd_pt +
// get_pc_now + AFSdelay). Puna ponovna procena resenja (bez O(1) delti).
//
// Konvencija cvorova (0-bazno): depo = 0, musterije = 1..nb, JEDNA punionica (AFS) = nb+1.
// Ruta = niz UNUTRASNJIH cvorova (musterije + eventualno AFS nb+1); depo je podrazumevan
// na oba kraja. Svaka ruta ima najvise JEDNU AFS posetu (kao u METS-u).
#pragma once
#include "instance.hpp"
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

using Route = std::vector<int>;        // unutrasnji cvorovi jedne rute (bez depoa)
using Solution = std::vector<Route>;

// Tezine kazni (Penalty_all red 1) + parametar V_nb. wM je uvek 0 (kao u METS-u).
struct Weights { double wT = 527, wC = 195, wD = 430, wM = 0; };

// AFS cvor je bilo koji label > nb (u nasem modelu tacno nb+1).
inline bool is_afs(const Instance& I, int node) { return node > I.n_customers; }

// ---- get_pd_pt: rastojanja/vremena pre- i post-AFS deonice jedne rute ----
struct LegInfo { double d1=0, d2=0, t1=0, t2=0; bool has_afs=false; };

inline LegInfo get_pd_pt(const Instance& I, const Route& r) {
    LegInfo L;
    if (r.empty()) return L;
    double ever = I.service_time, speed = I.speed;
    // pozicija AFS-a (max label); AFS ima najveci label
    int n = 0; int maxlab = r[0];
    for (int i = 0; i < (int)r.size(); ++i) if (r[i] >= maxlab) { maxlab = r[i]; n = i; }  // n = 0-baziran indeks AFS
    if (is_afs(I, maxlab)) {                        // ruta IMA punionicu na poziciji n
        L.has_afs = true;
        double D = 0, d1 = 0;
        for (int i = 0; i + 1 < (int)r.size(); ++i) {
            D += I.dist[r[i]][r[i+1]];
            if (i == n - 1) { d1 = D; D = 0; }      // zatvori pre-deonicu na hop ...->AFS
        }
        d1 += I.dist[0][r[0]];                      // depo -> prvi
        double d2 = D + I.dist[r.back()][0];        // poslednji -> depo
        L.d1 = d1; L.d2 = d2;
        L.t1 = d1 / speed + n * ever;               // (n) musterija pre AFS (0-baz: n cvorova pre)
        L.t2 = d2 / speed + ((int)r.size() - n) * ever; // post: (len - n) cvorova (ukljucuje AFS)
    } else {                                        // ruta BEZ punionice -> jedna deonica
        double D = 0;
        for (int i = 0; i + 1 < (int)r.size(); ++i) D += I.dist[r[i]][r[i+1]];
        L.d1 = D + I.dist[0][r[0]] + I.dist[r.back()][0];
        L.t1 = L.d1 / speed + (int)r.size() * ever;
        L.d2 = 0; L.t2 = 0;
    }
    return L;
}

// ---- get_pc_now: kazna za prekoracenje kapaciteta AFS-a (dato vreme pocetka dopune) ----
// starts = pocetna vremena dopune vozila (jedno po ruti sa AFS-om).
inline double get_pc_now(std::vector<double> starts, double T_Afs, int C_Afs) {
    starts.erase(std::remove(starts.begin(), starts.end(), 0.0), starts.end());
    int m = (int)starts.size();
    if (m == 0) return 0.0;
    std::vector<double> ev;                         // 2m dogadjaja (pocetak i kraj)
    for (double s : starts) { ev.push_back(s); ev.push_back(s + T_Afs); }
    std::sort(ev.begin(), ev.end());
    double pc = 0.0;
    const double eps = 1e-12;
    for (int j = 0; j + 1 < (int)ev.size(); ++j) {
        double left = ev[j], width = ev[j+1] - ev[j];
        if (width <= 0) continue;
        int cnt = 0;                                // koliko vozila se puni u intervalu [left, ...]
        for (double s : starts) { double d = s - left; if (eps < d + T_Afs && d <= 0) ++cnt; }
        int over = std::max(cnt - C_Afs, 0);
        pc += over * width;
    }
    return std::max(pc, 0.0);
}

// ---- AFSdelay: pokusaj da odlozis dolaske u AFS (unutar rezerve) da smanjis pc ----
// afs_arr = vreme dolaska u AFS po ruti; slack = raspoloziva rezerva (Tmax - time_V).
// Vraca rezidualni pc; menja 'afs_arr' na resene (odlozene) pocetke.
inline double afs_delay_resolve(std::vector<double>& afs_arr, const std::vector<double>& slack,
                                double T_Afs, int C_Afs) {
    int m = (int)afs_arr.size();
    if (m == 0) return 0.0;
    std::vector<double> start = afs_arr;            // pocetak dopune = dolazak (mozemo odloziti)
    std::vector<double> rem = slack;                // preostala rezerva po ruti
    double pc = get_pc_now(start, T_Afs, C_Afs);
    if (pc <= 1e-12) { afs_arr = start; return pc; }

    bool changed = true; int guard = 0;
    while (changed && guard++ < 10000) {
        changed = false;
        for (int i = 0; i < m && !changed; ++i)
            for (int j = i + 1; j < m; ++j) {
                double si = start[i], sj = start[j];
                double ei = si + T_Afs, ej = sj + T_Afs;
                if (!(si < ej && sj < ei)) continue;         // ne preklapaju se
                // koliko treba odloziti svako da izadje iza drugog
                double need_i = ej - si;                     // odlozi i iza j
                double need_j = ei - sj;                     // odlozi j iza i
                bool can_i = need_i <= rem[i] + 1e-9;
                bool can_j = need_j <= rem[j] + 1e-9;
                if (!can_i && !can_j) continue;              // par nerazresiv rezervom
                bool move_i;
                if (can_i && can_j) move_i = (rem[i] <= rem[j]); // odlozi onog sa manjom rezervom
                else move_i = can_i;
                if (move_i) { start[i] += need_i; rem[i] -= need_i; }
                else        { start[j] += need_j; rem[j] -= need_j; }
                changed = true;
            }
    }
    afs_arr = start;
    return get_pc_now(start, T_Afs, C_Afs);
}

// ---- Puna ocena resenja ----
struct EvalResult {
    double distance = 0;
    double pT = 0, pC = 0, pD = 0, pM = 0;          // NEtezinske kazne
    double cost = 0;                                // tezinska cena
    bool feasible = false;
    int n_routes = 0;
};

inline EvalResult evaluate(const Instance& I, const Solution& sol, const Weights& W) {
    EvalResult R; R.n_routes = (int)sol.size();
    std::vector<double> afs_arr, slack;
    for (const Route& r : sol) {
        if (r.empty()) continue;
        LegInfo L = get_pd_pt(I, r);
        double distV = L.d1 + L.d2;
        double timeV = L.t1 + L.t2;
        R.distance += distV;
        R.pD += std::max(L.d1 - I.d_max, 0.0) + std::max(L.d2 - I.d_max, 0.0);
        R.pT += std::max(timeV - I.t_max, 0.0);
        if (L.has_afs) { afs_arr.push_back(L.t1); slack.push_back(std::max(I.t_max - timeV, 0.0)); }
    }
    R.pC = afs_delay_resolve(afs_arr, slack, I.refuel_time, I.station_capacity.empty() ? 1 : I.station_capacity[0]);
    R.pM = std::max((double)R.n_routes - I.n_vehicles, 0.0);
    R.cost = R.distance + W.wT * R.pT + W.wC * R.pC + W.wD * R.pD + W.wM * R.pM;
    R.feasible = (R.pT == 0.0 && R.pC == 0.0 && R.pD == 0.0 && R.pM == 0.0);
    return R;
}

// Pomocno: da li ruta ima AFS
inline bool route_has_afs(const Instance& I, const Route& r) {
    for (int n : r) if (is_afs(I, n)) return true;
    return false;
}

// Ukloni suvisne AFS-ove: ako ruta bez AFS staje u Dmax, izbaci AFS (uslovno uklanjanje).
// Ako ruta bez AFS NE staje u Dmax a nema AFS, ubaci AFS na najbolju tacku (dve <=Dmax deonice).
inline void canonicalize_route(const Instance& I, Route& r) {
    // izbaci sve AFS pa odluci ponovo
    Route cust; for (int n : r) if (!is_afs(I, n)) cust.push_back(n);
    if (cust.empty()) { r.clear(); return; }
    // ukupno rastojanje bez AFS
    double full = I.dist[0][cust[0]];
    for (int i = 0; i + 1 < (int)cust.size(); ++i) full += I.dist[cust[i]][cust[i+1]];
    full += I.dist[cust.back()][0];
    if (full <= I.d_max + 1e-9) { r = cust; return; }        // ne treba AFS
    // treba AFS: nadji poziciju p (ubaci AFS posle cust[p]) da obe deonice budu <=Dmax,
    // biraj onu koja minimizuje ukupno rastojanje
    int afs = I.first_station();
    double best = 1e18; int bestp = 0; bool ok = false;
    for (int p = 0; p < (int)cust.size(); ++p) {
        Route t(cust.begin(), cust.begin() + p + 1);
        t.push_back(afs);
        t.insert(t.end(), cust.begin() + p + 1, cust.end());
        LegInfo L = get_pd_pt(I, t);
        double tot = L.d1 + L.d2;
        bool legok = (L.d1 <= I.d_max + 1e-9 && L.d2 <= I.d_max + 1e-9);
        if (legok && tot < best) { best = tot; bestp = p; ok = true; }
    }
    if (!ok) { // nijedna podela ne resava -> stavi AFS na tacku sa min rastojanjem svejedno
        for (int p = 0; p < (int)cust.size(); ++p) {
            Route t(cust.begin(), cust.begin() + p + 1);
            t.push_back(afs);
            t.insert(t.end(), cust.begin() + p + 1, cust.end());
            LegInfo L = get_pd_pt(I, t);
            double tot = L.d1 + L.d2;
            if (tot < best) { best = tot; bestp = p; }
        }
    }
    Route out(cust.begin(), cust.begin() + bestp + 1);
    out.push_back(afs);
    out.insert(out.end(), cust.begin() + bestp + 1, cust.end());
    r = out;
}

inline void canonicalize(const Instance& I, Solution& sol) {
    for (auto& r : sol) canonicalize_route(I, r);
    sol.erase(std::remove_if(sol.begin(), sol.end(), [](const Route& r){ return r.empty(); }), sol.end());
}

// ---- Zajednicka optimizacija pozicija AFS-a (KLJUCNO za kapacitet) ----
// Za date rute samo-musterija, biraj poziciju AFS-a po ruti (ili bez AFS-a) tako da
// se minimizuje ukupna cena UKLJUCUJUCI PC (vremensko rasporedjivanje na punionici).
// Time se izbegava vestacko zagusenje koje pravi "AFS na min-rastojanju".
inline Solution build_with_afs(const Instance& I, const std::vector<std::vector<int>>& cust,
                               const std::vector<int>& pos) {
    Solution s; int A = I.first_station();
    for (int r = 0; r < (int)cust.size(); ++r) {
        Route rt;
        for (int i = 0; i < (int)cust[r].size(); ++i) { if (pos[r] == i) rt.push_back(A); rt.push_back(cust[r][i]); }
        if (pos[r] == (int)cust[r].size()) rt.push_back(A);   // AFS na kraju
        s.push_back(rt);
    }
    return s;
}

// Da li ruta (samo-musterije) sa AFS na poziciji p ima obe deonice <= Dmax? (p=-1 => bez AFS)
inline bool afs_pos_ok(const Instance& I, const std::vector<int>& c, int p) {
    std::vector<int> tmp;
    int A = I.first_station();
    for (int i = 0; i < (int)c.size(); ++i) { if (p == i) tmp.push_back(A); tmp.push_back(c[i]); }
    if (p == (int)c.size()) tmp.push_back(A);
    LegInfo L = get_pd_pt(I, tmp);
    return L.d1 <= I.d_max + 1e-9 && L.d2 <= I.d_max + 1e-9;
}

// Vrati (resenje sa najboljim AFS pozicijama, ocena). Greedy koordinatni spust.
inline std::pair<Solution, EvalResult> eval_afs_opt(const Instance& I,
        const std::vector<std::vector<int>>& cust, const Weights& W) {
    int R = (int)cust.size();
    std::vector<int> pos(R, -1);
    // inicijalno: bez AFS ako ruta staje u Dmax, inace min-rastojanje izvodljiva pozicija
    for (int r = 0; r < R; ++r) {
        Route only = cust[r];
        LegInfo L0 = get_pd_pt(I, only);
        if (L0.d1 <= I.d_max + 1e-9) { pos[r] = -1; continue; }
        double best = 1e18; int bp = 0;
        for (int p = 0; p <= (int)cust[r].size(); ++p) if (afs_pos_ok(I, cust[r], p)) {
            Solution s = build_with_afs(I, cust, std::vector<int>(1, p)); // samo ova ruta
            LegInfo L = get_pd_pt(I, s[0]); double tot = L.d1 + L.d2;
            if (tot < best) { best = tot; bp = p; }
        }
        pos[r] = bp;
    }
    // BRZI IZLAZ: ako pocetni (min-rastojanje) raspored vec nema zagusenja (PC=0),
    // onda je on i optimalan za cenu -> preskoci skupi koordinatni spust.
    {
        Solution s0 = build_with_afs(I, cust, pos);
        EvalResult e0 = evaluate(I, s0, W);
        int afs_routes = 0; for (int p : pos) if (p >= 0) ++afs_routes;
        if (e0.pC <= 1e-9 || afs_routes < 2) return { s0, e0 };
    }
    // Koordinatni spust: menjaj AFS poziciju rute po rute da smanjis ukupnu cenu (sa PC).
    // OGRANICENJA RADI BRZINE (kljucno na 50-100 musterija, gde se ovo zove za svaki
    // kandidat-potez lokalne pretrage): max 3 prolaza, a za duge rute probaj samo
    // pozicije oko trenutne (+-2) plus krajeve -- timing se menja postepeno po
    // pozicijama, pa lokalni pomeraji nose skoro svu korist.
    auto cost_at = [&](const std::vector<int>& p){ return evaluate(I, build_with_afs(I, cust, p), W).cost; };
    double cur = cost_at(pos);
    bool improved = true; int guard = 0;
    while (improved && guard++ < 3) {
        improved = false;
        for (int r = 0; r < R; ++r) {
            int len = (int)cust[r].size();
            // kandidat pozicije: sve ako je ruta kratka; inace okolina trenutne + krajevi
            std::vector<int> cands;
            if (len <= 8) { for (int p = -1; p <= len; ++p) cands.push_back(p); }
            else {
                cands.push_back(-1); cands.push_back(0); cands.push_back(len);
                for (int d = -2; d <= 2; ++d) {
                    int p = pos[r] + d;
                    if (p >= 0 && p <= len) cands.push_back(p);
                }
            }
            int bestp = pos[r]; double bestc = cur;
            for (int p : cands) {
                if (p == pos[r]) continue;
                if (p >= 0 && !afs_pos_ok(I, cust[r], p)) continue;      // preskoci Dmax-nevalidne
                if (p == -1) { LegInfo L0 = get_pd_pt(I, cust[r]); if (L0.d1 > I.d_max + 1e-9) continue; }
                std::vector<int> trial = pos; trial[r] = p;
                double c = cost_at(trial);
                if (c < bestc - 1e-9) { bestc = c; bestp = p; }
            }
            if (bestp != pos[r]) { pos[r] = bestp; cur = bestc; improved = true; }
        }
    }
    Solution s = build_with_afs(I, cust, pos);
    return { s, evaluate(I, s, W) };
}
