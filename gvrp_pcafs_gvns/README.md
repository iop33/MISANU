# GrVRP-PCAFS — GVNS (C++ port)

C++ port našeg **General Variable Neighborhood Search (GVNS)** za Green Vehicle Routing
Problem with Private Capacitated Alternative Fuel Stations. Rešava **prave** benchmark
instance (Xu et al. 2025) i poredi rezultat sa objavljenim **BKS / GRASP / METS** brojevima.

Ovo je brza (C++) verzija Python prototipa iz projekta `MISANU`. Ista logika, isti
uslovi izvodljivosti (Reschedule / red na pumpi), ali ~50–100× brže.

## Poređenje sa METS-om (`mets_cpp`)

Ovaj projekat i `../mets_cpp` koriste **identičan protokol**, pa su tabele direktno
uporedive:
- iste instance (`data/instances/`, 40 Xu instanci) i **ista matrica rastojanja**
  (zaokružena na 2 decimale kao u METS-u),
- **5 ponavljanja** po instanci (skraćeno sa 30 iz rada),
- **isti vremenski budžeti** po setu (12/20/30/45 s po pokretanju),
- **kontrolisan seed**: run `r` = seed `base_seed + r` (isti seed ⇒ identičan
  rezultat; menja se sa `--seed N`).

Pokreni oba (svaki u svom CLion projektu) i uporedi kolone `OurBest/OurAvg` — obe
tabele prikazuju i objavljene BKS/GRASP/METS brojeve i gap%.

## Otvaranje u CLion-u

1. **File → Open** → izaberi ovaj folder (`gvrp_pcafs_gvns`).
2. CLion sam prepozna `CMakeLists.txt` i učita projekat.
3. **Build** (čekić) → **Run** (zeleni trougao).

**Klik na Run (bez argumenata) pokreće KOMPLETAN eksperiment**: svi setovi, 40
instanci, 5 ponavljanja — isti protokol kao `mets_cpp`. Traje ~1–1.5h.

### Argumenti (Run → Edit Configurations → Program arguments)
```
--set S-Central --n-runs 5 --time-limit 20
--seed 42                    # promeni bazni seed (reprodukcija/debagovanje)
```
- `--set` : `S-Central` | `M-Central25` | `M-Central50` | `M-Central100` | `all`
- `--n-runs` : broj nezavisnih pokretanja po instanci (default 5)
- `--time-limit` : sekundi po pokretanju (podrazumevano po setu ako se izostavi)
- `--instances-dir` : putanja do instanci (default se sam pronalazi)
- `--seed N` : bazni seed (default 1)

## Pokretanje iz terminala (bez CLion-a)
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/gvrp_pcafs_gvns --set S-Central --n-runs 8 --time-limit 25
```
(ili direktno: `clang++ -std=c++17 -O2 src/main.cpp -o gvrp && ./gvrp`)

## Struktura
```
CMakeLists.txt            # CMake konfiguracija (C++17, -O2)
src/
  instance.hpp            # model instance + čitanje .txt formata
  solution.hpp            # ruta/rešenje, ocena, Reschedule, kazne
  construction.hpp        # greedy, savings (Clarke&Wright), SCTS
  neighborhoods.hpp       # lokalna pretraga + shaking
  gvns.hpp                # glavni GVNS (adaptivna kazna, najbolje-izvodljivo)
  reference.hpp           # objavljeni BKS/GRASP/METS brojevi
  main.cpp                # CLI + uporedna tabela
data/instances/*.txt      # 40 pravih instanci (S-Central + M-Central 25/50/100)
```

## Šta ovo NIJE
Ovo je **naš** algoritam u C++. METS (Xu-ov algoritam) je u MATLAB-u i pokreće se odvojeno
(treba MATLAB). Ovde METS figuriše samo kroz **objavljene referentne brojeve** u `reference.hpp`.
Za poređenje „na istoj mašini", METS se pokreće u MATLAB-u, a ovaj binar je naša strana.
