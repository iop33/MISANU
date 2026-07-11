# METS u C++ — prevod Xu et al. (2025) algoritma

C++ prevod **METS** (Memetic Search) algoritma za *Green Vehicle Routing Problem with
Private Capacitated Alternative Fuel Stations (GrVRP-PCAFS)*, originalno objavljenog u
MATLAB-u uz rad:

> R. Xu, X. Fan, S. Liu, W. Chen, K. Tang. *Memetic Search for Green Vehicle Routing
> Problem With Private Capacitated Refueling Stations.* IEEE TEVC, 2025.
> Kod: https://github.com/FXBZ-research/METS-Algorithm

Ovaj projekat **prevodi METS iz MATLAB-a (6.438 linija) u C++** i pokreće ga nad
**Xu-ovim originalnim instancama** (S-Central, M-Central 25/50/100), pa poredi sa
objavljenim brojevima iz rada.

---

## Kako je prevedeno (i šta je verno)

METS-ov MATLAB kod održava skupe **O(1) inkrementalne procene** poteza. Ovaj prevod
umesto toga **potpuno ponovo procenjuje** svako rešenje posle poteza — isti rezultat
pretrage, mnogo manji rizik od grešaka. Prevedena je cela logika:

- **SCTS podela** giant tour-a (`split_Tmax`, `split_Dmax`) — `split.hpp`
- **Funkcija cilja i kazne** (PT/PC/PD, `chromR_detail_all` + `get_pd_pt` + `get_pc_now`
  + AFS raspoređivanje) — `eval.hpp`
- **Lokalna pretraga** (granularne okoline, potezi m1–m9 stila: relocate / or-opt /
  swap / 2-opt / 2-opt* / nova ruta), prvo-poboljšanje — `localsearch.hpp`
- **Memetska populacija** (feasible/infeasible, biased fitness, OX ukrštanje, binarni
  turnir) — `population.hpp`
- **Glavna petlja** (init populacija → GA petlja → adaptacija kazni → repair) — `metsalg.hpp`

**Ključna stavka — zajednički raspored AFS-a (kapacitet):** srce ovog problema je da
punionica prima najviše `C_Afs` vozila istovremeno. Zato se pozicija AFS-a u ruti bira
**vremenski svesno** (`eval_afs_opt` u `eval.hpp`), a ne po najkraćem rastojanju —
inače bi se sve dopune grupisale u istom trenutku i pravile veštačko zagušenje.
Ovo je bio kritičan detalj za verno reprodukovanje METS rezultata.

**Validacija:** pošto MATLAB-ov generator slučajnih brojeva ne može bit-po-bit da se
reprodukuje u C++, vernost se proverava **statistički** — da li prevod dostiže iste
najbolje rezultate kao rad.

### Rezultati na S-Central (15 mušt.) — 2 pokretanja × 15s

| Instanca | METS (rad) | Naš C++ | gap |
|---|---|---|---|
| S-Central_1 | 953.94 | 953.94 | 0.00% |
| S-Central_2 | 959.88 | 960.32 | +0.05% |
| S-Central_3 | 958.94 | 958.94 | 0.00% |
| S-Central_4 | 947.98 | 1098.49 | +15.9% ⚠️ |
| S-Central_5 | 714.55 | 714.55 | 0.00% |
| S-Central_6 | 844.43 | 844.43 | 0.00% |
| S-Central_7 | 862.68 | 866.66 | +0.46% |
| S-Central_8 | 712.83 | 712.83 | 0.00% |
| S-Central_9 | 855.43 | 855.43 | 0.00% |
| S-Central_10 | 905.59 | 906.76 | +0.13% |

**9/10 instanci reprodukovano tačno ili unutar 0.5%.** Prosek (bez `S-Central_4`) je
~0.08%. `S-Central_4` je najteža instanca (METS je tu našao specijalni novi BKS 947.98;
naš prevod daje 1098.49 = tačno GRASP-ov objavljeni broj) — traži više vremena/iteracija.

### Rezultati na M-Central25 (25 mušt.) — 2 pokretanja × 20s

Svih 10 instanci izvodljivo, **prosečan gap ka METS-u = 0.29%** (najveći +0.71%;
instanca `_8` je čak i malo bolja od objavljenog METS-a). Prevod je veran i na srednjim
instancama.

> Napomena: M-Central50/100 (50/100 mušt.) rade, ali su sporije — koristi veći
> `--time-limit`. Uzrok je skuplja vremenski-svesna procena AFS pozicija.

---

## Otvaranje i pokretanje u CLion-u

1. **File → Open** → izaberi folder `mets_cpp`. CLion prepozna `CMakeLists.txt`.
2. **Build** (čekić) → **Run** (zeleni trougao).

**Klik na Run (bez ikakvih argumenata) pokreće KOMPLETAN eksperiment** nad celim
dataset-om: sva 4 seta (S-Central + M-Central 25/50/100), svih 40 instanci,
**5 ponavljanja po instanci** (skraćeno sa 30 iz rada, po dogovoru sa mentorom),
i ispisuje uporedne tabele (naš rezultat vs objavljeni METS/BKS). Traje ~1–1.5h;
rezultati se ispisuju u hodu, instancu po instancu.

### Kontrola seed-a (za debagovanje)
Sva slučajnost je **deterministička**: run `r` koristi seed `base_seed + r`
(podrazumevano `base_seed = 1`). **Isti seed ⇒ bit-identičan rezultat** (provereno).
Bazni seed se menja opcijom `--seed N` — pa se svaki problematičan run može tačno
reprodukovati i debagovati.

### Opciono: argumenti (Run → Edit Configurations → Program arguments)
```
--inst S-Central_5           # brza provera JEDNE instance (~10s, očekivano 714.55)
--set S-Central --n-runs 5 --time-limit 20   # samo jedan set
--seed 42                    # promeni bazni seed (reprodukcija/debagovanje)
```
- `--set` : `S-Central` | `M-Central25` | `M-Central50` | `M-Central100` | `all`
- `--n-runs` : broj nezavisnih pokretanja po instanci (default 5)
- `--time-limit` : sekundi po pokretanju
- `--inst <ime>` : pokreni SAMO jednu instancu (ispisuje i seed po run-u)
- `--seed N` : bazni seed (default 1)

**Sve je u ovom folderu** — kod (`src/`), dataset (`data/instances/`, 40 instanci),
uputstvo (`README.md`), i specifikacije prevoda (`specs/`). Ništa spolja nije potrebno.

## Pokretanje iz terminala
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/mets_cpp --inst S-Central_5 --n-runs 3 --time-limit 20
```
(ili direktno: `clang++ -std=c++17 -O2 src/main.cpp -o mets && ./mets --set S-Central`)

---

## Struktura
```
CMakeLists.txt
src/
  instance.hpp      # ucitavanje .txt instance + matrica rastojanja (floor na 2 decimale kao METS)
  eval.hpp          # funkcija cilja, kazne (PT/PC/PD), AFS raspored (eval_afs_opt)
  split.hpp         # SCTS podela (split_Tmax, split_Dmax)
  localsearch.hpp   # granularna lokalna pretraga (m1..m9 stila)
  population.hpp    # memetska populacija, biased fitness, OX ukrstanje, turnir
  metsalg.hpp       # glavna METS petlja (init -> GA -> adaptacija kazni -> repair)
  reference.hpp     # objavljeni BKS/GRASP/METS brojevi (Tab. II-III)
  main.cpp          # CLI + uporedna tabela (default = pun eksperiment)
data/instances/     # 40 pravih Xu instanci (.txt), S-Central + M-Central 25/50/100
specs/              # detaljni specovi iz kojih je radjen prevod (referenca)
```

## Napomene o performansama
- Vremenski-svesan izbor AFS pozicija (`eval_afs_opt`) čini procenu skupljom. Na malim
  instancama (15–25 mušt.) je brzo; na 50/100 mušt. koristi `--time-limit` da ograničiš.
- Za bit-tačno poklapanje po seed-u trebao bi MATLAB-ov RNG; ovde se poklapa statistički
  (isti najbolji rezultati), što je i cilj poređenja „na istoj mašini".

## Instance
40 javnih instanci iz METS repozitorija, izvezene u jednostavan `.txt` format
(depo, mušterije, punionica, parametri). Rastojanja su euklidska, zaokružena na 2
decimale (kao u METS-u).
