# GrVRP-PCAFS: Green Vehicle Routing Problem with Private Capacitated Alternative Fuel Stations

## Overview

This project implements a **General Variable Neighborhood Search (GVNS)** metaheuristic for solving the Green Vehicle Routing Problem with Private Capacitated Alternative Fuel Stations (GrVRP-PCAFS), also known as GVRP-CAFS.

### Problem Description

Route a fleet of Alternative Fuel Vehicles (AFVs) to serve a set of customers while:
- Minimizing total travel distance
- Respecting limited driving range (vehicles must refuel at Alternative Fuel Stations)
- Respecting maximum route duration
- Managing **limited capacity at fuel stations** (only η_s vehicles can refuel simultaneously)

### Key References

1. **Bruglieri, M., Mancini, S. & Pisacane, O.** (2019). The green vehicle routing problem with capacitated alternative fuel stations. *Computers & Operations Research*, 112, 104759.
   - Introduces the GVRP-CAFS problem with MILP formulations

2. **Xu, R., Fan, X., Liu, S., Chen, W. & Tang, K.** (2025). Memetic Search for Green Vehicle Routing Problem with Private Capacitated Refueling Stations. *arXiv:2504.04527*.
   - State-of-the-art METS algorithm; benchmark instances

3. **Bruglieri, M., Ferone, D., Festa, P. & Pisacane, O.** (2025). An effective and efficient matheuristic for the electric vehicle routing problem with capacitated recharging stations. *International Transactions in Operational Research*.
   - MILP-based ALNS matheuristic for the EVRPTW-CRS variant

## Project Structure

```
gvrp_cafs/
├── instance.py              # Instance data model and I/O
├── solution.py              # Solution representation, evaluation, Reschedule procedure
├── construction.py          # Greedy and savings construction heuristics
├── neighborhoods.py         # Neighborhood structures (relocate, swap, 2-opt, etc.)
├── gvns.py                  # GVNS algorithm implementation
├── generate_instances.py    # Benchmark instance generator
├── main.py                  # Main experiment runner
├── results_analysis.py      # Results formatting, CSV/LaTeX export
└── README.md                # This file
```

## Algorithm: GVNS

The GVNS metaheuristic proceeds as follows:

1. **Initialization**: Greedy nearest-neighbor + savings-based construction
2. **Main loop** (until time limit):
   - **Shaking**: Random perturbation in neighborhood N_k
   - **Local search (VND)**: Sequential best-improvement with relocate, swap, 2-opt
   - **Move decision**: Accept if improvement, otherwise try next neighborhood
3. **Restart**: Perturbation when stuck

### Neighborhood Structures

**Shaking (diversification):**
- N1: Customer relocate (random)
- N2: Customer swap (random)
- N3: Or-opt (move segment of 1-3 customers)
- N4: 2-opt intra-route (reverse segment)
- N5: 2-opt* inter-route (exchange tails)
- N6: Station insert/remove/change

**Local Search (intensification):**
- Best relocate
- Best swap
- Best 2-opt intra-route

### AFS Capacity Management

The Reschedule procedure (from Xu et al., 2025) handles the limited station capacity:
- For each AFS, all vehicle visits are collected with arrival times
- If overlaps exceed station capacity, vehicles are delayed
- Feasibility is checked against T_max after rescheduling

## Benchmark Instance Sets

| Set | Customers | AFSs | Vehicles | η_s | T_max | D_max |
|-----|-----------|------|----------|-----|-------|-------|
| S-Central | 15 | 1 | 15 | 1 | 7h | 250 mi |
| M-Central25 | 25 | 1 | 7 | 2 | 7.5h | 250 mi |
| M-Central50 | 50 | 1 | 13 | 3 | 7.5h | 250 mi |
| M-Central100 | 100 | 1 | 25 | 8 | 7.5h | 250 mi |
| Triangle | 15 | 3 | 10 | 1 | 11h | 250 mi |
| EMH | 20 | 6 | 8 | 1 | 11h | 300 mi |

## Usage

### Quick Test
```bash
cd gvrp_cafs
python main.py --quick
```

### Full Experiment
```bash
python main.py --n-runs 5 --time-limit 120 --n-instances 10
```

### Specific Instance Set
```bash
python main.py --set S-Central --n-runs 10 --time-limit 300
```

### Command Line Options
```
--quick          Quick test mode (reduced parameters)
--set            Instance set to run (all, S-Central, M-Central25, etc.)
--n-runs         Number of independent runs per instance (default: 5)
--time-limit     Time limit per run in seconds (default: 120)
--n-instances    Number of instances per set (default: 10)
--output-dir     Output directory (default: results)
--quiet          Minimal output
```

## Output

Results are saved to the `results/` directory:
- `results/<set>/results.csv` - Results in CSV format
- `results/<set>/results_table.txt` - Formatted text table
- `results/<set>/results_latex.tex` - LaTeX table for paper
- `results/<set>/<instance>_solution.json` - Best solution details
- `results/summary.txt` - Summary across all sets
- `results/all_results.json` - All results in JSON format

## Requirements

- Python 3.8+
- NumPy
