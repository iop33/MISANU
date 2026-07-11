// instance.hpp -- Model podataka jedne GrVRP-PCAFS instance (C++ port).
// Cita nas jezik-neutralan .txt format (vidi data/instances/*.txt).
// Konvencija indeksa: 0 = depo, 1..n = musterije, n+1..n+s = punionice.
#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>

struct Instance {
    std::string name;
    int n_customers = 0;
    int n_stations  = 0;
    int n_vehicles  = 0;
    double speed = 0, tank_capacity = 0, consumption_rate = 0, d_max = 0;
    double t_max = 0, service_time = 0, refuel_time = 0, p_start = 0;
    std::vector<int> station_capacity;              // eta_s po punionici
    std::vector<double> xs, ys;                     // koordinate
    std::vector<std::vector<double>> dist;          // matrica rastojanja
    std::vector<std::vector<double>> tt;            // matrica vremena putovanja

    int n_total() const { return 1 + n_customers + n_stations; }
    bool is_customer(int i) const { return i >= 1 && i <= n_customers; }
    bool is_station(int i)  const { return i >= n_customers + 1 && i < n_total(); }
    bool is_depot(int i)    const { return i == 0; }
    int  first_station() const { return n_customers + 1; }
    int  station_cap(int node) const { return station_capacity[node - n_customers - 1]; }
    double svc(int node) const { return is_customer(node) ? service_time : 0.0; }
};

inline Instance load_instance(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Ne mogu da otvorim: " + path);
    // Tokenizuj ceo fajl (robustno na raspored po linijama).
    std::vector<std::string> tok;
    std::string w;
    while (in >> w) tok.push_back(w);

    Instance ins;
    std::vector<std::string> types;
    for (size_t i = 0; i < tok.size(); ) {
        const std::string& k = tok[i];
        auto num = [&](size_t j){ return std::stod(tok[j]); };
        if      (k == "NAME")        { ins.name = tok[i+1]; i += 2; }
        else if (k == "N_CUSTOMERS") { ins.n_customers = (int)num(i+1); i += 2; }
        else if (k == "N_STATIONS")  { ins.n_stations  = (int)num(i+1); i += 2; }
        else if (k == "N_VEHICLES")  { ins.n_vehicles  = (int)num(i+1); i += 2; }
        else if (k == "Q")           { ins.tank_capacity = num(i+1); i += 2; }
        else if (k == "R")           { ins.consumption_rate = num(i+1); i += 2; }
        else if (k == "DMAX")        { ins.d_max = num(i+1); i += 2; }
        else if (k == "SPEED")       { ins.speed = num(i+1); i += 2; }
        else if (k == "TMAX")        { ins.t_max = num(i+1); i += 2; }
        else if (k == "PSTART")      { ins.p_start = num(i+1); i += 2; }
        else if (k == "SERVICE")     { ins.service_time = num(i+1); i += 2; }
        else if (k == "REFUEL")      { ins.refuel_time = num(i+1); i += 2; }
        else if (k == "CAP")         { // n_stations celih brojeva
            i += 1;
            for (int s = 0; s < ins.n_stations; ++s) ins.station_capacity.push_back((int)num(i++));
        }
        else if (k == "NODES")       {
            i += 1;
            while (i < tok.size() && tok[i] != "EOF") {
                // idx type x y
                types.push_back(tok[i+1]);
                ins.xs.push_back(num(i+2));
                ins.ys.push_back(num(i+3));
                i += 4;
            }
        }
        else if (k == "EOF")         { i += 1; }
        else                          { i += 1; } // preskoci nepoznato
    }

    // Matrice rastojanja (euklidsko) i vremena.
    int n = ins.n_total();
    ins.dist.assign(n, std::vector<double>(n, 0.0));
    ins.tt.assign(n, std::vector<double>(n, 0.0));
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b) if (a != b) {
            double dx = ins.xs[a] - ins.xs[b], dy = ins.ys[a] - ins.ys[b];
            // METS zaokruzuje rastojanja na 2 decimale (floor(100*d)/100) -> radi tacnog poklapanja
            ins.dist[a][b] = std::floor(100.0 * std::sqrt(dx*dx + dy*dy)) / 100.0;
            ins.tt[a][b]   = ins.dist[a][b] / ins.speed;
        }
    if (ins.d_max <= 0) ins.d_max = ins.tank_capacity / ins.consumption_rate;
    return ins;
}
