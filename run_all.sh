#!/bin/zsh
# =============================================================================
# run_all.sh -- Pokreni OBA algoritma (GVNS + METS) nad celim dataset-om.
# -----------------------------------------------------------------------------
# Svaki algoritam SAM automatski cuva svoju analizu u results/<algo>_<vreme>.txt
# (ugradjeno u oba programa), pa posle ovog skripta odmah imas oba izvestaja.
#
# Koriscenje:
#   ./run_all.sh                       # pun eksperiment (svi setovi, 5 ponavljanja, ~2-3h ukupno)
#   ./run_all.sh --set S-Central       # brza proba samo na malom setu (~20 min)
#   ./run_all.sh --set S-Central --n-runs 1 --time-limit 5   # ekspresna provera (~2 min)
# Svi argumenti se prosledjuju OBEMA programima (isti protokol -> fer poredjenje).
# =============================================================================
set -e
ROOT="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$ROOT/build"

echo "== Kompajliranje (clang++, -O2) =="
clang++ -std=c++17 -O2 "$ROOT/gvrp_pcafs_gvns/src/main.cpp" -o "$ROOT/build/gvns_run"
clang++ -std=c++17 -O2 "$ROOT/mets_cpp/src/main.cpp"        -o "$ROOT/build/mets_run"

echo ""
echo "==================== 1/2: GVNS (nas algoritam) ===================="
(cd "$ROOT/gvrp_pcafs_gvns" && "$ROOT/build/gvns_run" "$@")

echo ""
echo "==================== 2/2: METS (Xu algoritam) ====================="
(cd "$ROOT/mets_cpp" && "$ROOT/build/mets_run" "$@")

echo ""
echo "== GOTOVO. Analize su u: $ROOT/results/ =="
ls -1t "$ROOT/results" | head -4
