# MISANU — GrVRP-PCAFS: naš GVNS vs Xu-ov METS (oba u C++)

Rešavanje **Green Vehicle Routing Problem with Private Capacitated Alternative Fuel
Stations (GrVRP-PCAFS)** — problem rutiranja vozila na alternativna goriva gde
punionica prima ograničen broj vozila istovremeno.

Repozitorijum sadrži **dva algoritma u C++**, na **istim instancama** i **istom
protokolu**, pa su rezultati direktno uporedivi:

| Folder | Algoritam | Poreklo |
|---|---|---|
| `gvrp_pcafs_gvns/` | **GVNS** (General Variable Neighborhood Search) — naš algoritam | prevod našeg Python prototipa (istorija u git-u) |
| `mets_cpp/` | **METS** (Memetic Search) — trenutni state-of-the-art | prevod MATLAB koda iz Xu et al., IEEE TEVC 2025 ([GitHub](https://github.com/FXBZ-research/METS-Algorithm)) |

Instance: **40 javnih benchmark instanci** iz Xu-ovog repozitorijuma
(S-Central 15 mušt. + M-Central 25/50/100), u `*/data/instances/`.

---

## Pokretanje

### Opcija A — jedan skript, oba algoritma (preporučeno)
```bash
./run_all.sh                     # pun eksperiment: svi setovi, 5 ponavljanja (~6h za oba)
./run_all.sh --set S-Central     # brza proba na malom setu (~20 min)
./run_all.sh --set S-Central --n-runs 1 --time-limit 5   # ekspresna provera (~2 min)
```
Izvrši **oba** algoritma sa istim parametrima. **Analize se automatski čuvaju** u
`results/gvns_<vreme>.txt` i `results/mets_<vreme>.txt`.

### Opcija B — CLion
1. **File → Open** → izaberi ovaj folder (`MISANU`). CLion učita koreni `CMakeLists.txt`.
2. U Run meniju biraš metu: `gvrp_pcafs_gvns` (naš) ili `mets_cpp` (Xu).
3. **Run** bez argumenata = pun eksperiment tog algoritma; analiza se **automatski
   čuva** u `results/` (program sam upisuje sve što ispiše).

### Argumenti (oba programa)
```
--set S-Central|M-Central25|M-Central50|M-Central100|all
--n-runs N          # ponavljanja po instanci (default 5, skraćeno sa 30 iz rada)
--time-limit S      # sekundi po pokretanju (default po setu: 12/20/60/120)
--seed N            # bazni seed (default 1); run r koristi seed N+r
```
`mets_cpp` dodatno: `--inst S-Central_5` (jedna instanca, brza provera vernosti —
očekivano ~714.55).

---

## Protokol eksperimenta (dogovoren sa mentorom)
- **Poređenje sa Xu-om**: obe tabele prikazuju objavljene **BKS / GRASP / METS**
  brojeve (Tabele II–III rada) i gap% naspram njih.
- **5 ponavljanja** po instanci (umesto 30 iz rada — kraće, a dovoljno za poređenje).
- **Kontrolisan seed**: sva slučajnost je deterministička; run `r` koristi seed
  `base+r`. Isti seed ⇒ bit-identičan rezultat (provereno) → lako debagovanje.
- **Ista matrica rastojanja** (euklidska, zaokružena na 2 decimale kao u METS-u) i
  **isti CPU budžet** po instanci za oba algoritma.

## Validacija METS prevoda (vernost)
C++ prevod METS-a reprodukuje objavljene brojeve: na S-Central **9/10 instanci
tačno ili unutar 0.5%** (npr. 714.55, 712.83, 953.94 — tačno na 2 decimale);
na M-Central25 prosečan gap **0.29%**. Detalji u `mets_cpp/README.md`;
specifikacije prevoda u `mets_cpp/specs/`.

## Struktura
```
CMakeLists.txt          # koreni build (obe mete)
run_all.sh              # pokreni oba algoritma + sacuvaj analize
gvrp_pcafs_gvns/        # nas GVNS (src/, data/instances/, README.md)
mets_cpp/               # Xu METS prevod (src/, data/instances/, specs/, README.md)
results/                # automatski sacuvane analize (verzionisu se na repou)
```

## Istorija
Python prototip (GVNS + učitavači instanci + eksperimenti) je u git istoriji —
u ranijim commit-ima na `main` grani.
