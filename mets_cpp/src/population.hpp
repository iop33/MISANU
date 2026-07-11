// population.hpp -- Memetska populacija (feasible/infeasible), biased fitness,
// OX ukrstanje, binarni turnir. Verno METS/HGS semantici (§4, §7 spec).
#pragma once
#include "instance.hpp"
#include "eval.hpp"
#include "localsearch.hpp"
#include <vector>
#include <random>
#include <algorithm>

struct Individual {
    SolC routes;                 // rute samo-musterija
    std::vector<int> tour;       // "giant tour" (redosled musterija) za ukrstanje
    double distance = 0, cost = 0;
    bool feasible = false;
    double fitness = 0;          // biased fitness (manje = bolje)
};

// redosled musterija (spljosti rute)
inline std::vector<int> flatten_tour(const SolC& s) {
    std::vector<int> t; for (auto& r : s) for (int c : r) t.push_back(c); return t;
}

inline Individual make_individual(const Instance& I, const SolC& routes, const Weights& W) {
    Individual ind; ind.routes = routes; ind.tour = flatten_tour(routes);
    EvalResult e = eval_customers(I, routes, W);
    ind.distance = e.distance; ind.cost = e.cost; ind.feasible = e.feasible;
    return ind;
}

// successor(c) = sledeca musterija u ruti (ili 0 ako je poslednja)
inline std::vector<int> succ_of(const Instance& I, const SolC& s) {
    std::vector<int> su(I.n_customers + 1, 0);
    for (auto& r : s) for (int p = 0; p + 1 < (int)r.size(); ++p) su[r[p]] = r[p+1];
    return su;
}
inline std::vector<int> pred_of(const Instance& I, const SolC& s) {
    std::vector<int> pr(I.n_customers + 1, 0);
    for (auto& r : s) for (int p = 1; p < (int)r.size(); ++p) pr[r[p]] = r[p-1];
    return pr;
}

// broken-pair rastojanje dve jedinke (udeo razlicitih grana)
inline double broken_dist(const Instance& I, const Individual& a, const Individual& b,
                          std::vector<std::vector<int>>& succ_cache, int ia, int ib) {
    // koristi kes successora
    const std::vector<int>& sa = succ_cache[ia];
    const std::vector<int>& sb = succ_cache[ib];
    // pred b
    // (jednostavna varijanta: razlika successora)
    int diff = 0;
    for (int c = 1; c <= I.n_customers; ++c) if (sa[c] != sb[c]) ++diff;
    return (double)diff / I.n_customers;
}

struct Population {
    std::vector<Individual> feas, infeas;
    const Instance* I; Weights W;
    int mu = 154, lambda = 68, eliteNum = 77, nClosest = 30;

    void recompute_fitness(std::vector<Individual>& pop) {
        int N = (int)pop.size();
        if (N == 0) return;
        // sort po ceni (rastuce)
        std::sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b){ return a.cost < b.cost; });
        if (N == 1) { pop[0].fitness = 0; return; }
        // kes successora
        std::vector<std::vector<int>> succ(N);
        for (int i = 0; i < N; ++i) succ[i] = succ_of(*I, pop[i].routes);
        // fitRank
        std::vector<double> fitRank(N);
        for (int i = 0; i < N; ++i) fitRank[i] = (double)i / (N - 1);
        // avgBrokenDist = -mean(nClosest najmanjih rastojanja)
        std::vector<double> avg(N);
        for (int i = 0; i < N; ++i) {
            std::vector<double> ds;
            for (int j = 0; j < N; ++j) if (j != i) {
                int diff = 0; for (int c = 1; c <= I->n_customers; ++c) if (succ[i][c] != succ[j][c]) ++diff;
                ds.push_back((double)diff / I->n_customers);
            }
            std::sort(ds.begin(), ds.end());
            int k = std::min((int)ds.size(), nClosest);
            double m = 0; for (int t = 0; t < k; ++t) m += ds[t];
            avg[i] = -(k > 0 ? m / k : 0);
        }
        // divRank: sort po avg (rastuce = najizolovaniji prvi), rank normiran
        std::vector<int> order(N); for (int i = 0; i < N; ++i) order[i] = i;
        std::sort(order.begin(), order.end(), [&](int a, int b){ return avg[a] < avg[b]; });
        std::vector<double> divRank(N);
        for (int r = 0; r < N; ++r) divRank[order[r]] = (double)r / (N - 1);
        double coef = 1.0 - (double)eliteNum / N;
        for (int i = 0; i < N; ++i) pop[i].fitness = fitRank[i] + coef * divRank[i];
    }

    void trim(std::vector<Individual>& pop) {
        while ((int)pop.size() > mu) {
            recompute_fitness(pop);
            // ukloni najgoru (najveci fitness), cuvaj najbolju (indeks 0 po ceni)
            int worst = -1; double wf = -1;
            for (int i = 1; i < (int)pop.size(); ++i) if (pop[i].fitness > wf) { wf = pop[i].fitness; worst = i; }
            if (worst < 0) break;
            pop.erase(pop.begin() + worst);
        }
        recompute_fitness(pop);
    }

    void add(const Individual& ind) {
        std::vector<Individual>& pop = ind.feasible ? feas : infeas;
        pop.push_back(ind);
        recompute_fitness(pop);
        if ((int)pop.size() > mu + lambda) trim(pop);
    }
};

// ---- OX ukrstanje (order crossover) na permutacijama musterija ----
inline std::vector<int> ox_crossover(const std::vector<int>& p1, const std::vector<int>& p2,
                                     int nb, std::mt19937& rng) {
    // p1,p2 permutacije 1..nb (mogu imati manje ako neka musterija fali; radimo na zajednickom skupu)
    std::vector<int> a = p1, b = p2;
    // osiguraj da su prave permutacije 1..nb
    std::vector<char> in(nb + 1, 0);
    std::vector<int> pa; for (int x : a) if (x >= 1 && x <= nb && !in[x]) { in[x] = 1; pa.push_back(x); }
    for (int x = 1; x <= nb; ++x) if (!in[x]) pa.push_back(x);
    std::fill(in.begin(), in.end(), 0);
    std::vector<int> pb; for (int x : b) if (x >= 1 && x <= nb && !in[x]) { in[x] = 1; pb.push_back(x); }
    for (int x = 1; x <= nb; ++x) if (!in[x]) pb.push_back(x);

    int n = nb;
    std::uniform_int_distribution<int> D(0, n - 1);
    int p = D(rng), q = D(rng); if (p > q) std::swap(p, q);
    std::vector<int> child(n, -1);
    std::vector<char> used(nb + 1, 0);
    for (int i = p; i <= q; ++i) { child[i] = pa[i]; used[pa[i]] = 1; }
    int idx = (q + 1) % n;
    for (int k = 0; k < n; ++k) {
        int src = pb[(q + 1 + k) % n];
        if (!used[src]) { child[idx] = src; used[src] = 1; idx = (idx + 1) % n; }
    }
    return child;
}

// binarni turnir nad kombinovanim skupom -> jedinka sa manjim fitness-om
inline const Individual& select_parent(const Population& P, std::mt19937& rng) {
    int a = (int)P.feas.size(), b = (int)P.infeas.size();
    auto pick = [&]() -> const Individual& {
        int r = rng() % (a + b);
        return (r < a) ? P.feas[r] : P.infeas[r - a];
    };
    const Individual& x = pick();
    const Individual& y = pick();
    return (x.fitness <= y.fitness) ? x : y;
}
