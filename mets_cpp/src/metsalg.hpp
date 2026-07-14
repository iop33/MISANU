// metsalg.hpp -- Glavna METS petlja (verni prevod Main_METS): inicijalna populacija
// (giant tour -> split -> lokalna pretraga), zatim genetska petlja (izbor roditelja ->
// OX ukrstanje -> split -> lokalna pretraga -> populacija), adaptacija kazni, repair.
#pragma once
#include "instance.hpp"
#include "eval.hpp"
#include "split.hpp"
#include "localsearch.hpp"
#include "population.hpp"
#include <vector>
#include <random>
#include <chrono>
#include <numeric>

struct MetsParams {
    int mu = 154, lambda = 68;
    double split_prob = 0.5;
    double PT = 527, PC = 195, PD = 430;
    double scaleUp = 1.2, scaleDown = 0.85, targetFeasible = 0.2;
    int nbLast = 20, maxIterNonProd = 300, maxIter = 2000;
    double timeLimit = 60.0;         // sekundi (dodato radi kontrole vremena)
};

struct MetsResult {
    Solution best;                   // najbolje IZVODLJIVO (sa AFS-om)
    double distance = 1e18;
    bool feasible = false;
    double time_to_best = 0;
    long iterations = 0;
};

inline MetsResult run_mets(const Instance& I, const MetsParams& P, unsigned seed) {
    using clock = std::chrono::steady_clock;
    auto t0 = clock::now();
    auto elapsed = [&]{ return std::chrono::duration<double>(clock::now() - t0).count(); };
    // TVRDI rok za CEO run: i lokalna pretraga ga postuje (prekida se usred prolaza).
    auto deadline = t0 + std::chrono::duration_cast<clock::duration>(
                             std::chrono::duration<double>(P.timeLimit));
    std::mt19937 rng(seed);

    Weights W; W.wT = P.PT; W.wC = P.PC; W.wD = P.PD; W.wM = 0;
    auto gran = build_granular(I, 5);

    Population pop; pop.I = &I; pop.W = W; pop.mu = P.mu; pop.lambda = P.lambda;
    pop.eliteNum = (int)(0.5 * P.mu); pop.nClosest = (int)(0.2 * P.mu);

    MetsResult R;
    auto update_best = [&](const Individual& ind) {
        if (ind.feasible && ind.distance < R.distance - 1e-9) {
            R.distance = ind.distance; R.feasible = true; R.time_to_best = elapsed();
            R.best = eval_customers_full(I, ind.routes, W).first;   // resenje sa najboljim AFS pozicijama
            return true;
        }
        return false;
    };

    // prozor poslednjih nbLast jedinki (za adaptaciju kazni) -> pamtimo po-dimenziji izvodljivost
    std::vector<std::array<bool,3>> last_window;   // {pT==0, pC==0, pD==0}
    auto push_window = [&](const EvalResult& e) {
        last_window.push_back({e.pT == 0, e.pC == 0, e.pD == 0});
        if ((int)last_window.size() > P.nbLast) last_window.erase(last_window.begin());
    };

    // repair: privremeno 10x tezine prekrsenih dimenzija, ponovo LS; ako izvodljivo -> u populaciju
    auto try_repair = [&](SolC routes) {
        EvalResult e0 = eval_customers(I, routes, W);
        if (e0.feasible) return;
        Weights Wr = W;
        if (e0.pT > 0) Wr.wT *= 10; if (e0.pC > 0) Wr.wC *= 10; if (e0.pD > 0) Wr.wD *= 10;
        local_search(I, routes, Wr, gran, rng, deadline);
        Individual ind = make_individual(I, routes, W);
        if (ind.feasible) { pop.add(ind); update_best(ind); }
    };

    std::uniform_real_distribution<double> U01(0, 1);

    // ---- Inicijalna populacija: 4*mu nasumicnih permutacija ----
    std::vector<int> base(I.n_customers); std::iota(base.begin(), base.end(), 1);
    int init_n = 4 * P.mu;
    for (int i = 0; i < init_n && (int)R.iterations < P.maxIter; ++i) {
        if (elapsed() > P.timeLimit) break;
        std::vector<int> tsp = base; std::shuffle(tsp.begin(), tsp.end(), rng);
        Solution split = (U01(rng) < P.split_prob) ? split_Dmax(I, tsp) : split_Tmax(I, tsp);
        SolC routes; for (auto& r : split) { std::vector<int> cust; for (int n : r) if (!is_afs(I, n)) cust.push_back(n); if (!cust.empty()) routes.push_back(cust); }
        local_search(I, routes, W, gran, rng, deadline);
        Individual ind = make_individual(I, routes, W);
        pop.add(ind); update_best(ind);
        push_window(eval_customers(I, routes, W));
        if (!ind.feasible && U01(rng) < 0.5) try_repair(routes);
        R.iterations++;
    }

    // ---- Glavna genetska petlja ----
    int nbIterNonProd = 0;
    while ((int)R.iterations < P.maxIter && nbIterNonProd <= P.maxIterNonProd && elapsed() <= P.timeLimit) {
        if (pop.feas.empty() && pop.infeas.empty()) break;
        const Individual& p1 = select_parent(pop, rng);
        const Individual& p2 = select_parent(pop, rng);
        std::vector<int> off = ox_crossover(p1.tour, p2.tour, I.n_customers, rng);
        Solution split = (U01(rng) < 0.5) ? split_Dmax(I, off) : split_Tmax(I, off);
        SolC routes; for (auto& r : split) { std::vector<int> cust; for (int n : r) if (!is_afs(I, n)) cust.push_back(n); if (!cust.empty()) routes.push_back(cust); }
        local_search(I, routes, W, gran, rng, deadline);
        Individual ind = make_individual(I, routes, W);
        pop.add(ind);
        bool newbest = update_best(ind);
        push_window(eval_customers(I, routes, W));
        if (!ind.feasible && U01(rng) < 0.5) try_repair(routes);
        nbIterNonProd = newbest ? 0 : nbIterNonProd + 1;
        R.iterations++;

        // ---- Adaptacija kazni na svakih nbLast iteracija ----
        if (R.iterations % P.nbLast == 0 && (int)last_window.size() >= P.nbLast) {
            double fT = 0, fC = 0, fD = 0; int n = (int)last_window.size();
            for (auto& w : last_window) { fT += w[0]; fC += w[1]; fD += w[2]; }
            fT /= n; fC /= n; fD /= n;
            auto adj = [&](double frac, double& w) {
                if (frac <= P.targetFeasible - 0.05) w = std::min(100000.0, w * P.scaleUp);
                else if (frac >= P.targetFeasible + 0.05) w = std::max(0.1, w * P.scaleDown);
            };
            adj(fT, W.wT); adj(fC, W.wC); adj(fD, W.wD);
            pop.W = W;
            // ponovo oceni cene neizvodljivih (zavise od tezina)
            for (auto& ind2 : pop.infeas) { EvalResult e = eval_customers(I, ind2.routes, W); ind2.cost = e.cost; }
            pop.recompute_fitness(pop.infeas);
        }
    }
    return R;
}
