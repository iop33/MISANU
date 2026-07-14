// main.cpp -- METS (C++ prevod) nad Xu-ovim instancama; poredjenje sa objavljenim brojevima.
//
// Koriscenje:
//   ./mets_cpp --set S-Central --n-runs 5 --time-limit 30
//   --set: S-Central | M-Central25 | M-Central50 | M-Central100 | all
//   --n-runs: broj nezavisnih pokretanja po instanci
//   --time-limit: sekundi po pokretanju
//   --inst 15_5: (opciono) pokreni SAMO jednu instancu (za brzu proveru vernosti)
#include "instance.hpp"
#include "metsalg.hpp"
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

static std::string find_data_dir(const std::string& given) {
    std::vector<std::string> cands = { given, "data/instances", "../data/instances",
                                       "../../data/instances", "../../../data/instances",
                                       // MISANU raspored (projekat je podfolder repo-a):
                                       "mets_cpp/data/instances",
                                       "../mets_cpp/data/instances",
                                       "../../mets_cpp/data/instances",
                                       "../../../mets_cpp/data/instances",
                                       "/Users/matejacivkaroski/Documents/CLionProjects/MISANU/mets_cpp/data/instances" };
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

int main(int argc, char** argv) {
    std::string set = "S-Central", data_dir = "data/instances", single;
    // 5 ponavljanja po instanci (umesto 30 iz rada -- skraceno po dogovoru sa mentorom).
    int n_runs = 5; double time_limit = -1;
    // KONTROLA SEED-a: run r koristi seed = base_seed + r (deterministicno -> lako debagovanje;
    // isti seed = identican rezultat). Promena preko --seed.
    unsigned base_seed = 1;
    // BEZ ARGUMENATA (klik na Run u CLion-u) -> PUN eksperiment nad CELIM dataset-om.
    if (argc == 1) {
        set = "all";
        std::cout << "[Run bez argumenata -> PUN eksperiment: svi setovi (S-Central + M-Central 25/50/100),\n"
                  << " sve 40 instanci, 5 ponavljanja po instanci (skraceno sa 30 iz rada). Traje ~3h.]\n"
                  << "[Seed: deterministicki, run r = seed " << base_seed << "+r; promena: --seed N]\n"
                  << "[Za brzu proveru jedne instance: Program arguments = --inst S-Central_5]\n\n";
    }
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--set" && i+1 < argc) set = argv[++i];
        else if (a == "--n-runs" && i+1 < argc) n_runs = std::stoi(argv[++i]);
        else if (a == "--time-limit" && i+1 < argc) time_limit = std::stod(argv[++i]);
        else if (a == "--instances-dir" && i+1 < argc) data_dir = argv[++i];
        else if (a == "--inst" && i+1 < argc) single = argv[++i];
        else if (a == "--seed" && i+1 < argc) base_seed = (unsigned)std::stoul(argv[++i]);
    }
    data_dir = find_data_dir(data_dir);

    // ---- Ukljuci automatsko cuvanje analize u results/ (koren MISANU repo-a) ----
    // data_dir = .../mets_cpp/data/instances  ->  koren = data_dir/../../..
    static std::ofstream results_file;
    static TeeBuf tee(nullptr, nullptr);
    static CoutRestore restore{ std::cout.rdbuf() };
    {
        std::string results_dir = data_dir + "/../../../results";
        std::error_code ec; std::filesystem::create_directories(results_dir, ec);
        std::time_t tnow = std::time(nullptr);
        char stamp[32]; std::strftime(stamp, sizeof(stamp), "%Y-%m-%d_%H-%M-%S", std::localtime(&tnow));
        results_file.open(results_dir + "/mets_" + std::string(stamp) + ".txt");
        if (results_file) {
            tee.a = std::cout.rdbuf(); tee.b = results_file.rdbuf();
            std::cout.rdbuf(&tee);
            std::cout << "[Analiza se automatski cuva u: results/mets_" << stamp << ".txt]\n\n";
        }
    }

    auto default_tl = [](const std::string& s){
        if (s=="S-Central") return 12.0; if (s=="M-Central25") return 20.0;
        if (s=="M-Central50") return 60.0; return 120.0; };

    // brza provera jedne instance (npr. --inst S-Central_5 ocekivano ~714.55)
    if (!single.empty()) {
        std::string path = data_dir + "/" + single + ".txt";
        if (!file_exists(path)) { std::cerr << "nema " << path << "\n"; return 1; }
        Instance ins = load_instance(path);
        MetsParams P; P.timeLimit = (time_limit > 0) ? time_limit : 30.0;
        std::cout << "METS na " << single << " (" << ins.n_customers << " must., "
                  << n_runs << "x" << (int)P.timeLimit << "s)\n";
        auto ref = reference_table().find(single);
        for (int r = 0; r < n_runs; ++r) {
            MetsResult res = run_mets(ins, P, base_seed + r);          // kontrolisan seed
            std::cout << "  run " << r << " (seed=" << base_seed + r << "): "
                      << (res.feasible ? std::to_string(res.distance) : std::string("INF"))
                      << "  (iter=" << res.iterations << ", t2best=" << std::fixed << std::setprecision(1) << res.time_to_best << "s)\n";
        }
        if (ref != reference_table().end())
            std::cout << "  [ref] BKS=" << ref->second[0] << " GRASP=" << ref->second[1] << " METS=" << ref->second[2] << "\n";
        return 0;
    }

    std::vector<std::string> sets;
    if (set == "all") sets = {"S-Central","M-Central25","M-Central50","M-Central100"};
    else sets = {set};

    for (auto& sname : sets) {
        double tl = (time_limit > 0) ? time_limit : default_tl(sname);
        std::cout << "\n>>> METS(C++) " << sname << ": 10 instanci, " << n_runs << "x" << (int)tl << "s\n";
        std::cout << std::string(96,'=') << "\n";
        std::cout << std::left << std::setw(16) << "Instanca" << std::right
                  << std::setw(10) << "BKS" << std::setw(10) << "METS(rad)"
                  << std::setw(11) << "OurBest" << std::setw(11) << "OurAvg" << std::setw(10) << "gapMETS%" << std::setw(6) << "feas" << "\n";
        std::cout << std::string(96,'-') << "\n";
        std::cout << std::fixed << std::setprecision(2);
        std::vector<double> gaps;
        for (int idx = 1; idx <= 10; ++idx) {
            std::string name = sname + "_" + std::to_string(idx);
            std::string path = data_dir + "/" + name + ".txt";
            if (!file_exists(path)) continue;
            Instance ins = load_instance(path);
            MetsParams P; P.timeLimit = tl;
            std::cout << std::left << std::setw(16) << name << std::flush; // ime -> ziva info pre racunanja
            std::vector<double> feas;
            for (int r = 0; r < n_runs; ++r) {
                MetsResult res = run_mets(ins, P, base_seed + r);      // kontrolisan seed
                if (res.feasible) feas.push_back(res.distance);
            }
            auto ref = reference_table().find(name);
            double bks = ref != reference_table().end() ? ref->second[0] : 0;
            double mets = ref != reference_table().end() ? ref->second[2] : 0;
            std::cout << std::right << std::setw(10) << bks << std::setw(10) << mets;
            if (!feas.empty()) {
                double best = *std::min_element(feas.begin(), feas.end());
                double avg = 0; for (double d : feas) avg += d; avg /= feas.size();
                double g = mets ? 100.0*(best-mets)/mets : 0;
                if (mets) gaps.push_back(g);
                std::cout << std::setw(11) << best << std::setw(11) << avg << std::setw(10) << g;
            } else std::cout << std::setw(11) << "INF" << std::setw(11) << "-" << std::setw(10) << "-";
            std::cout << std::setw(4) << (int)feas.size() << "/" << n_runs << "\n";
        }
        std::cout << std::string(96,'-') << "\n";
        double mg = 0; for (double g : gaps) mg += g; if (!gaps.empty()) mg /= gaps.size();
        std::cout << std::left << std::setw(16) << "PROSEK gapMETS" << std::right << std::setw(52) << mg << "%\n";
        std::cout << std::string(96,'=') << "\n";
    }
    return 0;
}
