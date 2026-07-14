// main.cpp -- Pokretac: nas GVNS (C++) na PRAVIM instancama vs objavljeni METS/GRASP/BKS.
//
// Koriscenje:
//   ./gvrp_pcafs_gvns --set S-Central --n-runs 8 --time-limit 25
//   --set: S-Central | M-Central25 | M-Central50 | M-Central100 | all
#include "instance.hpp"
#include "solution.hpp"
#include "gvns.hpp"
#include "reference.hpp"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <ctime>
#include <sys/stat.h>

static bool file_exists(const std::string& p){ struct stat st; return stat(p.c_str(), &st) == 0; }

// Nadji folder sa instancama (razni radni direktorijumi: CLion build dir, koren MISANU...).
static std::string find_data_dir(const std::string& given) {
    std::vector<std::string> cands = { given, "data/instances", "../data/instances",
                                       "../../data/instances", "../../../data/instances",
                                       // MISANU raspored (projekat je podfolder repo-a):
                                       "gvrp_pcafs_gvns/data/instances",
                                       "../gvrp_pcafs_gvns/data/instances",
                                       "../../gvrp_pcafs_gvns/data/instances",
                                       "../../../gvrp_pcafs_gvns/data/instances",
                                       "/Users/matejacivkaroski/Documents/CLionProjects/MISANU/gvrp_pcafs_gvns/data/instances" };
    for (auto& c : cands) if (!c.empty() && file_exists(c + "/S-Central_1.txt")) return c;
    return given;
}

// ---- AUTOMATSKO CUVANJE ANALIZE ----
// TeeBuf: sve sto ide na ekran (std::cout) upisuje se ISTOVREMENO i u fajl.
struct TeeBuf : std::streambuf {
    std::streambuf *a, *b;
    TeeBuf(std::streambuf* x, std::streambuf* y) : a(x), b(y) {}
    int overflow(int c) override { if (c != EOF) { a->sputc((char)c); b->sputc((char)c); } return c; }
    int sync() override { a->pubsync(); b->pubsync(); return 0; }
};
// Vraca cout na originalni bafer pre gasenja programa (bezbedno unistavanje).
struct CoutRestore { std::streambuf* old; ~CoutRestore(){ std::cout.rdbuf(old); } };

struct Row {
    std::string name; double bks=0, grasp=0, mets=0;
    bool has_feas=false; double best=0, avg=0, sd=0, worst=0; int feas=0, runs=0, nroutes=0;
};

int main(int argc, char** argv) {
    std::string set = "S-Central", data_dir = "data/instances";
    // ISTI PROTOKOL kao mets_cpp: 5 ponavljanja, kontrolisan seed, isti vremenski budzeti.
    int n_runs = 5; double time_limit = -1;
    unsigned base_seed = 1;                       // run r koristi seed = base_seed + r
    // BEZ ARGUMENATA (klik na Run u CLion-u) -> PUN eksperiment nad CELIM dataset-om.
    if (argc == 1) {
        set = "all";
        std::cout << "[Run bez argumenata -> PUN eksperiment: svi setovi, 40 instanci, 5 ponavljanja.\n"
                  << " Isti protokol kao mets_cpp -> rezultati direktno uporedivi. Traje ~3h.]\n"
                  << "[Seed: deterministicki, run r = seed " << base_seed << "+r; promena: --seed N]\n\n";
    }
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--set" && i+1 < argc) set = argv[++i];
        else if (a == "--n-runs" && i+1 < argc) n_runs = std::stoi(argv[++i]);
        else if (a == "--time-limit" && i+1 < argc) time_limit = std::stod(argv[++i]);
        else if (a == "--instances-dir" && i+1 < argc) data_dir = argv[++i];
        else if (a == "--seed" && i+1 < argc) base_seed = (unsigned)std::stoul(argv[++i]);
    }
    data_dir = find_data_dir(data_dir);

    // ---- Ukljuci automatsko cuvanje analize u results/ (koren MISANU repo-a) ----
    // data_dir = .../gvrp_pcafs_gvns/data/instances  ->  koren = data_dir/../../..
    static std::ofstream results_file;
    static TeeBuf tee(nullptr, nullptr);
    static CoutRestore restore{ std::cout.rdbuf() };
    {
        std::string results_dir = data_dir + "/../../../results";
        std::error_code ec; std::filesystem::create_directories(results_dir, ec);
        std::time_t tnow = std::time(nullptr);
        char stamp[32]; std::strftime(stamp, sizeof(stamp), "%Y-%m-%d_%H-%M-%S", std::localtime(&tnow));
        results_file.open(results_dir + "/gvns_" + std::string(stamp) + ".txt");
        if (results_file) {
            tee.a = std::cout.rdbuf(); tee.b = results_file.rdbuf();
            std::cout.rdbuf(&tee);
            std::cout << "[Analiza se automatski cuva u: results/gvns_" << stamp << ".txt]\n\n";
        }
    }

    std::vector<std::string> sets;
    if (set == "all") sets = {"S-Central","M-Central25","M-Central50","M-Central100"};
    else sets = {set};
    // Vremenski budzeti IDENTICNI mets_cpp (fer poredjenje: isti CPU budzet po instanci).
    auto default_tl = [](const std::string& s){
        if (s=="S-Central") return 12.0; if (s=="M-Central25") return 20.0;
        if (s=="M-Central50") return 60.0; return 120.0; };

    for (auto& sname : sets) {
        double tl = (time_limit > 0) ? time_limit : default_tl(sname);
        std::cout << "\n>>> " << sname << ": 10 instanci, " << n_runs
                  << " pokretanja x " << (int)tl << "s  (dir: " << data_dir << ")\n";
        std::vector<Row> rows;
        std::vector<double> gaps_bks, gaps_mets;

        for (int idx = 1; idx <= 10; ++idx) {
            std::string name = sname + "_" + std::to_string(idx);
            std::string path = data_dir + "/" + name + ".txt";
            if (!file_exists(path)) { std::cerr << "  (nedostaje " << path << ")\n"; continue; }
            Instance ins = load_instance(path);
            std::cout << "  [" << idx << "/10] " << name << " ..." << std::flush;

            std::vector<double> feas;
            int nroutes = 0;
            for (int run = 0; run < n_runs; ++run) {
                GvnsResult r = gvns(ins, tl, 400, 6, base_seed + run);   // kontrolisan seed
                if (r.feasible) { feas.push_back(r.distance); nroutes = r.n_routes; }
            }
            Row row; row.name = name; row.runs = n_runs; row.feas = (int)feas.size(); row.nroutes = nroutes;
            auto ref = reference_table().find(name);
            if (ref != reference_table().end()) { row.bks = ref->second[0]; row.grasp = ref->second[1]; row.mets = ref->second[2]; }
            if (!feas.empty()) {
                row.has_feas = true;
                row.best = *std::min_element(feas.begin(), feas.end());
                row.worst = *std::max_element(feas.begin(), feas.end());
                double m = 0; for (double d : feas) m += d; m /= feas.size();
                double v = 0; for (double d : feas) v += (d-m)*(d-m); v /= feas.size();
                row.avg = m; row.sd = std::sqrt(v);
                if (row.bks) gaps_bks.push_back(100.0*(row.best-row.bks)/row.bks);
                if (row.mets) gaps_mets.push_back(100.0*(row.best-row.mets)/row.mets);
            }
            rows.push_back(row);
            std::cout << (row.has_feas ? "  best=" + std::to_string((int)row.best) : "  INF") << "\n";
        }

        // tabela
        std::cout << "\n" << std::string(100,'=') << "\n";
        std::cout << "  NAS GVNS (C++) vs METS/GRASP/BKS  |  set: " << sname << "  (gap>0 = mi smo losiji)\n";
        std::cout << std::string(100,'=') << "\n";
        std::cout << std::left << std::setw(16) << "Instanca" << std::right
                  << std::setw(9) << "BKS" << std::setw(9) << "GRASP" << std::setw(9) << "METS"
                  << std::setw(10) << "OurBest" << std::setw(12) << "OurAvg" << std::setw(9) << "gapBKS%"
                  << std::setw(10) << "gapMETS%" << std::setw(8) << "feas" << "\n";
        std::cout << std::string(100,'-') << "\n";
        std::cout << std::fixed << std::setprecision(2);
        for (auto& r : rows) {
            std::cout << std::left << std::setw(16) << r.name << std::right
                      << std::setw(9) << r.bks << std::setw(9) << r.grasp << std::setw(9) << r.mets;
            if (r.has_feas) {
                std::cout << std::setw(10) << r.best << std::setw(12) << r.avg
                          << std::setw(9) << (r.bks?100.0*(r.best-r.bks)/r.bks:0)
                          << std::setw(10) << (r.mets?100.0*(r.best-r.mets)/r.mets:0);
            } else std::cout << std::setw(10) << "INF" << std::setw(12) << "-" << std::setw(9) << "-" << std::setw(10) << "-";
            std::cout << std::setw(6) << r.feas << "/" << r.runs << "\n";
        }
        std::cout << std::string(100,'-') << "\n";
        auto mean = [](std::vector<double>& v){ double m=0; for(double x:v)m+=x; return v.empty()?0:m/v.size(); };
        std::cout << std::left << std::setw(16) << "PROSEK gap" << std::right << std::setw(45) << ""
                  << std::setw(9) << mean(gaps_bks) << std::setw(10) << mean(gaps_mets) << "\n";
        int ninf = 0; for (auto& r : rows) if (!r.has_feas) ninf++;
        if (ninf) std::cout << "  ! " << ninf << " instanci bez izvodljivog resenja\n";
        std::cout << std::string(100,'=') << "\n";
    }
    return 0;
}
