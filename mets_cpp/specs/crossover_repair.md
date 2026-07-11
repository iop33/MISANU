I now have a thorough understanding of the entire system: the three target files plus every dependency needed to make the spec self-contained. I have enough to write the exhaustive specification. I'll write it as structured plain text.

Since the task explicitly says "Do NOT Write report/summary/findings/analysis .md files" and to "Return findings directly as your final assistant message," I'll return the full spec as my final response text.

---

# METS → C++ PORT SPECIFICATION

Scope: this spec is the sole implementation reference for `Crossover.m`, `selectparents.m`, and `Repair_sol.m`, plus every data structure, penalty, objective, and control-flow element those three files touch. All supporting behavior needed to implement them faithfully is included. MATLAB is 1-indexed; all conversions to 0-indexed C++ are called out explicitly.

Source files read in full: `Crossover.m`, `selectparents.m`, `Repair_sol.m`, `get_chromR.m`, `phrase_chromR.m`, `chromR_detail_all.m`, `get_pd_pt.m`, `get_pc_now.m`, `AFSdelay_new.m`, `AFSdelay_recursion.m`, `deleteAFS_node.m`, `m1.m`, `ELS_mian.m`, `Main_METS.m`, `PopManagement.m`, `add2Pop.m`, `infeasiblePop_updateBiasedFitnesses.m`, `uti_updateBestSol.m`, `uti_addSol2Last100.m`, `get_vrp.m`.

---

## 1. NODE NUMBERING AND INDEXING CONVENTIONS

There are **two distinct numbering schemes** in this codebase. Getting this right is the single most important thing.

### 1.1 The `vrp.distance_table` / "detail" scheme (1-indexed, depot = 1)
Used everywhere `vrp.distance_table` is indexed, and in the `chromR` (route cell) structure produced by the segmentation `split_*` and consumed by `chromR_detail_all`.
- Index `1` = the depot.
- Index `c+1` = customer `c` (customers are `1..nb_customer`). So customer 1 is table-index 2, etc.
- Indices `> nb_customer + 1` = AFS (alternative fuel station) nodes. `nb_customer+2` is the first AFS table index. In practice the model uses a single physical AFS location whose distance row/column is table index `nb_customer+2` (see `dAll(1,nbClients+2)`, `dAll(vrp.nb_customer+2,nodeY+1)` in m1). Each *visit* to an AFS is a distinct node, but they all share the same distance-table row/column `nb_customer+2`.

So a `chromR{i}` route is a **column vector** like `[1; c1+1; c2+1; ...; 1]` — it **starts and ends at depot index 1**, customer values are `customer+1`, AFS visit values are `> nb_customer+1`.

### 1.2 The "move" scheme (0-indexed relative, depot = 0), used inside ELS / linked-list
`chromR_move` = `chromR` with **1 subtracted from every element and all resulting zeros (the depot markers) removed**. So a move-route holds customers as their raw id `1..nb_customer`, and AFS visits as `> nb_customer`. Depot is implicit (value 0, removed). This is the representation the crossover, the linked-list, and the local search operate on.

Conversion, both directions, appears repeatedly:
- detail→move: `a = chromR{i} - 1; a(a==0) = []` (drop depot).
- move→table index when indexing `distance_table`: `distance_table(node+1, ...)` — i.e. add 1 back. Depot is `distance_table(1,...)`.

### 1.3 The linked-list arrays (the core solution state used by ELS)
Produced by `phrase_chromR` from `chromR_move`, and carried through ELS as `Node_related` and `Route_related`. This is the canonical mutable solution representation. Four parallel arrays over all nodes:

- `predecessors[k]` — predecessor node id of node `k` (0 = depot / route start). Length = total number of nodes = `nb_customer + (number of AFS visits)`.
- `successor[k]` — successor node id of node `k` (0 = depot / route end).
- `routeID[k]` — 1-based route number that node `k` belongs to.
- `node_location[k]` — position marker relative to the route's AFS: **-1** = node is before the AFS (or route has no AFS at all), **1** = node is after the AFS, **100** = this node *is* the AFS. (See §6.5.)

**Node index layout.** Indices `1..nb_customer` are the customers (customer id == its index). Indices `nb_customer+1 .. end` are the AFS-visit nodes, one per AFS visit in the solution; the AFS node for route `r` (when built by ELS bookkeeping) is at index `r + nb_customer` in the canonical ordering. When operators add an AFS they append a new node at index `numel(node_location)+1`.

`Node_related` is the matrix `[predecessors, successor, routeID, node_location']` (columns 1..4). `Route_related` is a per-route matrix, one row per route (see §3.2).

**C++ recommendation:** represent a solution as (a) `vector<vector<int>> chromR_move` (routes of raw ids, customers `1..N`, AFS `>N`), plus derived linked-list arrays. Because the task says to fully re-evaluate every candidate rather than do incremental deltas, the linked-list arrays are needed mainly to know route membership and node_location; you can rebuild `chromR_move` from them via `get_chromR` (§6.1) and recompute all costs via `chromR_detail_all` logic (§3).

---

## 2. THE `vrp` STRUCT (problem instance) — fields actually used

All are loaded from a `.mat` file (`get_vrp` just `switch`es on an integer 1..60 and loads `vrp`). Fields referenced by the code we must port:

- `vrp.nb_customer` (int) — number of customers `N`.
- `vrp.last_customer` (int) — used only as a normalizing denominator in broken-pair distance (§4.4). Treat as a scalar from the instance (equal to `nb_customer` in effect).
- `vrp.distance_table` (matrix, size `(N+1+#AFS) × (N+1+#AFS)`) — symmetric distance matrix in the 1-indexed depot-first scheme (§1.1). `distance_table(a,b)` = distance between table-nodes `a` and `b`.
- `vrp.V_speed` (double) — vehicle speed. Travel time = distance / speed.
- `vrp.T_Customer` (double) — per-customer service time (`everTime`).
- `vrp.T_Afs` (double) — refueling duration at an AFS.
- `vrp.C_Afs` (int) — AFS capacity = number of vehicles that may refuel simultaneously.
- `vrp.V_Dmax` (double) — max distance a vehicle may travel **on a single leg between depot/AFS boundaries** (per-segment distance cap; drives PD).
- `vrp.T_max_V` (double) — max route duration per vehicle (drives PT).
- `vrp.V_nb` (int) — target/allowed number of vehicles (drives PM).
- `vrp.correlatedVertices` (matrix `N × nbGranular`) — granular neighbor lists: for each customer, its `nbGranular` nearest correlated customers (values are customer ids `1..N`). Used to enumerate candidate `nodeV` per `nodeU` in ELS.
- `vrp.ALL_brokenDIS` (matrix) — pairwise broken-pair distances among **infeasible** population members (maintained incrementally). Preallocated width `popSizeLambda+popSizeMu+1`.
- `vrp.ALL_brokenDIS_feasible` (matrix) — same for the **feasible** population.

Penalty base weights are **not** in `vrp`; they live in `Penalty_all` and `par_hgs` (§3.1, §5).

---

## 3. OBJECTIVE, PENALTIES, AND THE `Penalty_all` / `Route_related` / `sol_table` STRUCTURES

### 3.1 `Penalty_all` — 2×4 matrix
- Row 1 (weights, the `w*`): `[wT, wC, wD, wM]` = `[PT, PC, PD, 0]` initially. `wM` starts at 0 and stays 0 (vehicle-count penalty is effectively disabled by weight = 0, but the code still computes it).
- Row 2 (weighted penalty values): `[penalty_T, penalty_C, penalty_D, penalty_m]` where each is `w* × raw_penalty`.
- Initial values (from `Main_METS`): `PT = 527, PC = 195, PD = 430, wM = 0`.

### 3.2 `Route_related` — one row per route, 7 columns
Column meaning (fixed throughout):
- Col 1,2 = `penalty_D_v` = `[max(d1 - V_Dmax, 0), max(d2 - V_Dmax, 0)]` where `d1` = distance of the pre-AFS segment, `d2` = distance of the post-AFS segment (0 if no AFS).
- Col 3 = `penalty_T_v` = `max(time_V - T_max_V, 0)` for that route.
- Col 4 = `afs_time` = arrival time at the AFS on that route (0 if no AFS), **without** waiting/delay.
- Col 5 = `time_V` = total route duration.
- Col 6,7 = `distance_pre_su` = `[d1, d2]` (pre-AFS distance, post-AFS distance).

### 3.3 Per-route distance & time (`get_pd_pt` and `chromR_detail_all`)
Given a `route_now` in **move scheme** (customers `1..N`, AFS `>N`, depot implicit at both ends):

Let the max element be `a` at position `n` (`[a,n] = max(route_now)`).

**Case A — route contains an AFS (`a > nb_customer`):** the AFS split is at position `n`.
- `d1` (pre-AFS) = `distance_table(1, route_now(1)+1)` + sum of consecutive `distance_table(route_now(i)+1, route_now(i+1)+1)` for `i = 1..n-1`. (Depot → ... → AFS.)
- `d2` (post-AFS) = sum of consecutive leg distances for `i = n..end-1` + `distance_table(route_now(end)+1, 1)`. (AFS → ... → depot.)
- `t1` = `d1/speed + (n-1)*everTime` (service times for the `n-1` nodes before the AFS; the AFS itself and the pre-AFS depot are not customers).
- `t2` = `d2/speed + (numel(route_now) - n + 1)*everTime`.
- Route total time = `t1 + t2`. **Note**: `t1`+`t2` counts `everTime` for `numel(route_now)` positions including the AFS position; in `chromR_detail_all` the total route time is computed separately as `time_V = distance_V/speed + (numel(full_route_with_depots) - 2)*T_Customer` — see below; there is a subtle difference in service-time accounting between the two code paths, so reproduce each path exactly where it is used.

**Case B — no AFS (`a <= nb_customer`):**
- `d1` = full loop distance: `sum of legs` + `distance_table(1, route_now(1)+1)` + `distance_table(route_now(end)+1, 1)`.
- `t1` = `d1/speed + numel(route_now)*everTime`.
- `d2 = 0`, `t2 = 0`.

**In `chromR_detail_all`, the authoritative per-route time** (`time_V`) is computed on the **detail-scheme** route `a = chromR{i}` (with depot 1s at both ends): `distance_V = sum of distance_table(a(j),a(j+1))`, `time_V = distance_V/V_speed + (numel(a)-2)*T_Customer` (the `-2` excludes the two depot endpoints; note the AFS node **is** counted as a service stop here, i.e. it adds a `T_Customer`). `overtime_V = max(time_V - T_max_V, 0)`.

`afs_time` (arrival time at AFS, no wait) in `chromR_detail_all`: build `time_window(i,j)` = `distance_table(a(j-1),a(j))/V_speed + T_Customer` for each edge into node `j` (but **no** `T_Customer` when the node is the depot `==1`). Then `afs_time(i) = sum(time_window(i, 1:afs_position)) - T_Customer` where `afs_position` is the index (within the route) of the AFS node.

### 3.4 The four penalties (raw, then weighted)

**PT — time / T_max penalty.**
`penalty_T = max(overtime_Total, 0)` where `overtime_Total = sum over routes of max(time_V - T_max_V, 0)`.
Weighted: `Penalty_all(2,1) = wT * penalty_T`. Also stored per-route as `penalty_T_v = max(time_V - T_max_V, 0)` (Route_related col 3).

**PD — per-segment distance / V_Dmax penalty.**
For each route, split into segments at AFS boundaries: segment boundaries `z = [1, afs_positions..., route_end]`. For each segment, `D` = sum of its leg distances; `penalty_D_v(route, seg) = max(D - V_Dmax, 0)`. (A route with no AFS is one segment = the whole route.)
`penalty_D = sum(sum(penalty_D_v))`. Weighted: `Penalty_all(2,3) = wD * penalty_D`.

**PM — vehicle-count penalty.**
`penalty_m = max(nb_V - vrp.V_nb, 0)` where `nb_V = number of routes`. Weighted: `Penalty_all(2,4) = wM * penalty_m`. Since `wM = 0`, this is always 0 but must still be computed and carried.

**PC — AFS-capacity (simultaneous-refueling) penalty.** This is the most involved. See §3.5–§3.6.

### 3.5 PC base computation (`get_pc_now`) — given the finalized arrival-delay times
Inputs: `afs_time_delay` (per-route arrival time at AFS, **after** any delay scheduling; 0 means route has no AFS), `T_Afs`, `C_Afs`.
1. Drop zero entries from `afs_time_delay` (routes with no AFS).
2. `afs_time_end = afs_time_delay + T_Afs` (each refueling finishes `T_Afs` after arrival).
3. `c = sort([afs_time_delay, afs_time_end])` — merge all start/end timestamps, ascending.
4. `time_during = diff(c)` — width of each elementary interval between consecutive events.
5. For each interval, count how many vehicles are refueling: a vehicle with start `s` is refueling in interval starting at `c(k)` iff `eps < (s - c(k)) + T_Afs` AND `(s - c(k)) <= 0`, i.e. `s <= c(k) < s + T_Afs`. `nb_fueling(k) = ` that count.
6. `nb_fueling = max(nb_fueling - C_Afs, 0)` — vehicles exceeding capacity.
7. `pc_now = max( nb_fueling · time_duringᵀ , 0 )` — sum over intervals of (excess vehicles × interval width). This is "vehicle-minutes of over-capacity".

Weighted: `Penalty_all(2,2) = wC * pc_now`.

`chromR_detail_all` computes the identical quantity inline (lines 244–279) using `afs_time_delay` after the conflict-resolution pass; it counts with `-T_Afs < d & d <= 0` (note the strict/loose boundaries differ slightly from `get_pc_now`'s `eps < d+T_Afs & d <= 0` — reproduce each exactly in its own path). When there are no AFS at all, `penalty_C = 0`.

### 3.6 AFS delay scheduling (`AFSdelay_new` → `AFSdelay_recursion`)
The vehicles that visit the (shared) AFS may be scheduled to **wait** (delay their arrival) to reduce simultaneous over-capacity, as long as the delay fits within each route's slack. This computes both the resulting `pc_now` **and** the adjusted `afs_time_delay` array.

`AFSdelay_new(time_afs, time_v, T_max, T_Afs, C_Afs, afs_time_delay)`:
- If `time_afs` is empty → `pc_now = 0`, return.
- `time_V_shifting = max(T_max - time_v, 0)` — per-route slack (how much a route can be delayed before it violates T_max).
- `afs_time_delay = time_afs` (start from un-delayed arrivals).
- Call `AFSdelay_recursion`.

`AFSdelay_recursion` (pairwise conflict resolution with backtracking):
- Compute initial `pc_now = get_pc_now(...)`.
- Repeat (`delay_over` loop) scanning all ordered pairs `(i, ii)`, `i<ii`:
  - Order the pair so `t1 = smaller arrival`, `t2 = larger arrival`; `q1,q2` = corresponding slacks; `t3` records which was which (1 = route i is the later one; 2 = route i is the earlier one).
  - If `t1 + T_Afs <= t2` → no overlap, `continue`.
  - If **both** `t1+T_Afs-t2 > q2` **and** `t2+T_Afs-t1 > q1` → unresolvable conflict; bump `conflict_table`, `continue` (this pair's overlap is unavoidable and will be charged in PC).
  - Else if only `t1+T_Afs-t2 > q2` → must delay the earlier one: set the earlier vehicle's `afs_time_delay` to `t2 + T_Afs`, reduce its slack, restart the scan (`i=1; ii=1`).
  - Else if only `t2+T_Afs-t1 > q1` → must delay the later one similarly.
  - Else (both delays feasible) → the fourth `if` (`t1+T_Afs-t2 <= q2 && t2+T_Afs-t1 <= q1`) plus the duplicated-condition branch triggers `backtrack(...)`, which tries delaying vehicle 1 then vehicle 2, recursively re-solving, and keeps whichever yields lower `pc_now`. Set `delay_over=1`, break, and rescan.
- After stabilization, recompute `pc_now = get_pc_now(...)` on the final delays and return.

**Porting note:** this is a heuristic delay-assignment search. Reproduce the branch structure and ordering exactly; there are known quirks (e.g. `i=1; ii=1` reassignments inside nested for-loops don't reset the MATLAB loop counters — MATLAB `for` ignores in-loop index reassignment — so those lines have **no effect on iteration**; only the `afs_time_delay`/`time_V_shifting` updates and the `continue` matter). `conflict_table` is computed but its value is not returned/used by callers; you may omit it.

`chromR_detail_all` (lines 136–221) performs an equivalent one-pass (non-recursive) version of the same pairwise delay logic to initialize `afs_time_delay` for a freshly built solution. Reproduce that path for the initial-solution route detailing.

### 3.7 Total cost & feasibility
- `distance_Total = sum of all route distances` = `sum(sum(Route_related(:,[6,7])))` (pre + post segment distances over all routes).
- `cost_Total = distance_Total + sum(Penalty_all(2,:))` = `distance_Total + wT·pT + wC·pC + wD·pD + wM·pM`.
- `IsFeasible = 1` iff `sum(Penalty_all(2,:)) == 0`, i.e. all weighted penalties are zero (equivalently all raw penalties zero). Otherwise `IsFeasible = 0`.
- (In `chromR_detail_all`: `IsFeasible = 1` iff `penalty_T + penalty_C + penalty_D + penalty_m > 0` is **false**.)

### 3.8 `sol_table` — the per-individual record (fields the ported code touches)
A row per individual (`ID = tspid`). Fields set by `chromR_detail_all` and `ELS_mian` that downstream code (the three target files + pop management) reads:
- `ID` (int) = tspid.
- `chromR` (cell of route column-vectors, detail scheme).
- `chromR_move` (cell of route vectors, move scheme). ← used by `Crossover` indirectly and by `Repair_sol` (passes `sol_table.chromR_move{tspid}`) and by final result extraction in Main.
- `predecessors`, `successor`, `routeID` (column vectors, linked-list). ← used by `Crossover` (`get_chromR`), by broken-pair distance, by biased fitness.
- `node_location` (row vector).
- `distance_Total`, `cost_Total` (doubles).
- `penalty_T, penalty_C, penalty_D, penalty_m` (weighted doubles). ← `penalty_T/C/D` used by the penalty-adaptation loop and feasibility fractions.
- `IsFeasible` (0/1).
- `Fitness` (biased fitness, double), `fitRank`, `divRank`, `avgBrokenDist`, `brokenPairDistance` (set by pop management).
- `time` (double, wall-clock at insertion).
- Plus many diagnostic fields (`distance_window`, `time_window`, `afs_time`, `afs_time_delay`, `overtime_V`, etc.) not needed for the three target files but produced by `chromR_detail_all`.

---

## 4. BIASED FITNESS, POPULATION STRUCTURE, BROKEN-PAIR DISTANCE

### 4.1 Two subpopulations
`feasiblePop` and `infeasiblePop`, each an array (MATLAB table) of `sol_table` rows, **kept sorted by ascending `cost_Total`** (with ties broken so that a newer/higher-`ID` clone sorts earlier — see §5.6). Sizes are capped: after insertion, if a subpop exceeds `popSizeMu + popSizeLambda` it is trimmed back down to `popSizeMu` by removing worst-biased-fitness / clone individuals (§4.5).

### 4.2 Ranks
For a subpop of size `n` sorted by cost ascending:
- `fitRank(i)` (cost rank, normalized): assign `1..n` by cost order, then `fitRank = (rank-1)/(n-1)` ∈ [0,1]. Lower is better (cheaper).
- `divRank(i)` (diversity rank, normalized): sort individuals by `avgBrokenDist` ascending, assign `1..n`, normalize `(rank-1)/(n-1)` ∈ [0,1]. Lower rank = **more** diverse (because `avgBrokenDist` is stored negated — see §4.3).

### 4.3 `avgBrokenDist`
For each individual, `avgBrokenDist = -mean( the maxSize smallest broken-pair distances to the other members )`, where `maxSize = min(nClosest, n-1)`, `nClosest = floor(0.2 · popSizeMu)`. Negated so that "smaller value ⇒ closer to others ⇒ less diverse."

### 4.4 Biased fitness formula
```
Fitness(i) = fitRank(i) + (1 - eliteNum/n) * divRank(i)
```
where `eliteNum = floor(0.5 · popSizeMu)`, `n = current subpop size`. **Lower Fitness = better.** This is the value `selectparents` compares and the value pop-management uses to pick the worst individual to delete.

### 4.5 Survivor selection / clone removal (`PopManagement`, when size > μ+λ)
While `n > popSizeMu`:
- Scan `i = 2..n` (index 1, the best, is never removed).
- `isClone(i)` = `sum(mink(brokenPairDistance(i,:), 2)) < 1e-8` — the two smallest pairwise distances (which include distance-to-self = 0 plus nearest neighbor) sum below threshold ⇒ i is essentially a duplicate of some other member.
- Track the worst: prefer to delete a clone; among equals (both clone or both non-clone) delete the one with the **largest** (worst) `Fitness`. Rule:
  `if (isClone && !isWorstIndividualClone) || (isClone==isWorstIndividualClone && Fitness(i) > worstFitness)` → update worst = i.
- Remove that individual, remove its row/col from the broken-pair matrix, repeat.
- After trimming, recompute `fitRank`, `avgBrokenDist`, `divRank`, `Fitness` for the survivors.

### 4.6 Broken-pair distance between two solutions
`calculateBrokenPairDistance(solA, solB)` (also inlined in `add2Pop`):
1. Take `successor`/`predecessor` arrays of both; truncate both to the common length `a = min(len_A, len_B)` (AFS counts may differ).
2. Map every AFS reference to a single canonical id: any value `> nb_customer` becomes `nb_customer+1`.
3. For each node `jj` in `1..a`:
   - `cc(jj)=1` if `succA(jj) != succB(jj)` **and** `succA(jj) != predB(jj)` (the pair (jj, its successor) is "broken" — not present in the other solution in either orientation).
   - `dd(jj)=1` if `predA(jj)==0` (jj starts a route in A) **and** `predB(jj)!=0` **and** `succB(jj)!=0` (but jj is mid-route in B).
4. `distance = (sum(cc) + sum(dd)) / vrp.last_customer`.

Range ~[0, ~2]. Distance 0 ⇒ identical structure (clone). Maintained incrementally in `vrp.ALL_brokenDIS[_feasible]`; for a full C++ re-eval you may recompute pairwise on demand.

---

## 5. CONTROL FLOW — `Main_METS` driver (the context the three files run in)

### 5.1 Constants (exact)
```
split_prob            = 0.5
PT = 527, PC = 195, PD = 430           (PM weight wM = 0)
penaltyScaleFactor    = 1.2
penaltyDecreaseFactor = 0.85
popSizeMu             = 154
popSizeLambda         = 68
targetFeasible        = 0.2
nbLast                = 20
maxIterNonProd        = 300
maxIter               = 2000
timeLimit             = 100000   (seconds; effectively unbounded)
el = 0.5   → eliteNum   = floor(0.5 * 154)   = 77
nc = 0.2   → nClosest   = floor(0.2 * 154)   = 30
nbGranular = 20
```
`Penalty_all = [PT PC PD 0; 0 0 0 0]`.

### 5.2 RNG
`rng(SEED + tspid)` is re-seeded at specific points (init loop, per-ELS-call, main loop). To reproduce results exactly you must mirror MATLAB's Mersenne-Twister stream, which is generally infeasible in C++. For a faithful-but-not-bit-identical port, use a seeded MT19937 and re-seed with `SEED+tspid` at the same call sites. The algorithm's correctness does not depend on exact RNG, only its statistical behavior.

### 5.3 Initial population construction (loop `i = 1 .. popSizeMu*4`, i.e. 616)
`tsp_all` = `popSizeMu*4` random permutations of `1..N` (each row a giant-tour customer ordering). For each `i` (with `tspid = i`):
1. Break if `toc > timeLimit` or `tspid > maxIter`.
2. `tsp = tsp_all(i,:)`.
3. **Segmentation ("Split"):** with prob `split_prob` call `split_Dmax(vrp, tsp, par_hgs)`, else `split_Tmax(...)`. Produces `chromR` (routes, detail scheme) — a Split that partitions the giant tour into vehicle routes, inserting AFS where needed (respecting Dmax or Tmax respectively). (Bodies not in the three target files; treat as: partition `tsp` into routes each `[1 ... 1]` with conditional AFS insertion so segment distance ≤ V_Dmax / route time ≤ T_max_V.)
4. `chromR_detail_all` → fills `sol_table(tspid)`, updates `Penalty_all` row 2, returns `node_location`, `Route_related`.
5. `ELS_mian` → local search improves the individual (§6), rewrites `sol_table(tspid)`, `Penalty_all`, `Route_related`.
6. `PopManagement` → insert into feasible/infeasible subpop, recompute biased fitness, trim if oversized.
7. `uti_updateBestSol` → update `bestSolRestart`/`bestSolOverall`.
8. `uti_addSol2Last100` → push into the sliding `Last100Sol` window (capacity `nbLast = 20`; FIFO).
9. **Repair** (§5.5): if `IsFeasible==0` and `rand-0.5 > 1e-6` (≈ 50% probability, biased slightly toward triggering), call `Repair_sol`.

### 5.4 Main loop (`tspid = popSizeMu*4+1 .. maxIter`)
Terminate when `nbIterNonProd > maxIterNonProd (300)` or `toc > timeLimit`. Each iteration:
1. `p1 = selectparents(...)`, `p2 = selectparents(...)` (§7).
2. `offspring_tsp = Crossover(p1, p2, vrp)` (§8) — an offspring **giant-tour permutation** of customers.
3. `rng(SEED+tspid)`.
4. Segmentation on `tsp` (**note a source bug: it splits `tsp`, the stale variable from the init loop, not `offspring_tsp`** — reproduce faithfully: the Split uses `tsp`, but the offspring is threaded into `chromR_detail_all` as its `tsp` argument). Practically, `chromR` comes from splitting `tsp`; `offspring_tsp` is recorded as the individual's `tsp`.
5. `chromR_detail_all` → detail the solution.
6. `ELS_mian` → local search.
7. `PopManagement`, `uti_updateBestSol`, `uti_addSol2Last100`.
8. Repair if infeasible (same condition as §5.3.9).
9. `nbIterNonProd`: reset to 1 if `isNewBest`, else `+1`.
10. **Penalty adaptation** (§5.7): every `nbLast` iterations once `tspid >= 100`.

### 5.5 `Repair_sol` — EXACT (target file)
Signature (MATLAB): `[isrepair, feasiblePop, infeasiblePop, sol_table, bestSolRestart, bestSolOverall, isNewBest, vrp] = Repair_sol(par_hgs, sol_table, vrp, tspid, test, isrepair, feasiblePop, infeasiblePop, bestSolRestart, bestSolOverall, nbClients, Penalty_all, Route_related, seednum)`.

Logic (pseudocode):
```
isrepair = 1
isNewBest = 0
WP = 10                                   # penalty multiplier for repair
for each column i of Penalty_all:         # i = 1..4
    if Penalty_all(2,i) > 0:              # this penalty is currently active (>0)
        Penalty_all(:,i) = Penalty_all(:,i) * WP   # multiply BOTH weight (row1) and value (row2) by 10
# Re-run local search on the SAME individual, now with amplified penalties,
# in "repair mode" (isrepair=1). Passing sol_table.chromR_move{tspid} as the route set,
# and sol_table.node_location{tspid}. In repair mode ELS does NOT re-derive chromR_move
# (skips the "-1 / drop depot" step) because chromR is already in move scheme.
[~, sol_table] = ELS_mian(sol_table, /*sol_individual*/1, vrp, tspid, par_hgs, test,
                          isrepair=1, sol_table.chromR_move{tspid}, nbClients,
                          sol_table.node_location{tspid}, Penalty_all, Route_related, seednum)
if sol_table.IsFeasible(tspid) == 1:      # repair succeeded → now feasible
    [feasiblePop, infeasiblePop, vrp] = PopManagement(sol_table(tspid,:), par_hgs, tspid,
                                                      feasiblePop, infeasiblePop, vrp)
    [bestSolRestart, bestSolOverall, isNewBest] =
        uti_updateBestSol(bestSolRestart, bestSolOverall, sol_table(tspid,:))
isrepair = 0
```
Key semantics for the port:
- Repair = **local search again with each active penalty weight ×10**, forcing the search to drive penalties toward 0. Only currently-active penalties (row-2 value > 0) are amplified; inactive ones untouched.
- The amplified `Penalty_all` is **local** to Repair (MATLAB pass-by-value; not returned to Main), so the main loop's `Penalty_all` is unchanged after repair. Only `sol_table`, populations, and best-sols persist.
- In repair mode, `ELS_mian` is entered with `isrepair=1`, which (per `ELS_mian` lines 5–14) means `chromR_move = chromR` is used **as-is** (no `-1`/depot-drop), because the passed `chromR` is `sol_table.chromR_move{tspid}` already in move scheme.
- If still infeasible after repair, the individual is **not** inserted anywhere (it was already inserted pre-repair in the caller's PopManagement call; note the improved-but-still-infeasible `sol_table(tspid)` row is left updated but not re-inserted).

### 5.6 `uti_updateBestSol` — EXACT
```
if bestSolRestart empty: bestSolRestart.cost_Total = 99999999
if bestSolOverall empty: bestSolOverall.cost_Total = 99999999
if sol.IsFeasible AND sol.cost_Total < bestSolRestart.cost_Total - 1e-9:
    bestSolRestart = sol
    if sol.cost_Total < bestSolOverall.cost_Total - 1e-9:
        bestSolOverall = sol
    isNewBest = true
else:
    isNewBest = false
```
Only **feasible** solutions can become best. `bestSolOverall` initialized with `IsFeasible=0, cost_Total=999999`.

### 5.7 Penalty adaptation (main loop, `mod(tspid, nbLast)==0 && tspid>=100`)
Over the sliding window `Last100Sol` (last ≤20 individuals):
```
fractionFeasible_T = (# window members with penalty_T == 0) / window_size
fractionFeasible_C = ... penalty_C == 0 ...
fractionFeasible_D = ... penalty_D == 0 ...
```
For each of T, C, D independently (with `target = targetFeasible = 0.2`):
```
if fraction <= target - 0.05:      # too few feasible on this constraint → tighten
    w = min(100000, w * penaltyScaleFactor)        # ×1.2, capped at 1e5
    rescale every infeasiblePop member's penalty_* by (old_penalty/origin_w)*new_w
elif fraction >= target + 0.05:    # too many feasible → loosen
    w = max(0.1, w * penaltyDecreaseFactor)        # ×0.85, floored at 0.1
    rescale infeasiblePop penalty_* similarly
```
`origin_PT/PC/PD` capture the pre-update weights so stored per-individual weighted penalties can be rescaled to the new weight. After updating all three:
- Recompute every infeasible member's `cost_Total = penalty_D + penalty_C + penalty_T + distance_Total`.
- If `>1` infeasible member: recompute biased fitness (`infeasiblePop_updateBiasedFitnesses`), re-sort by `cost_Total` (ties: higher ID first — done via negating ID), and rebuild `vrp.ALL_brokenDIS`.

**Note:** weights only ratchet for the **infeasible** population's stored penalties and the global `Penalty_all(1,:)`. The feasible population has all penalties 0 so is unaffected.

---

## 6. `ELS_mian` — LOCAL SEARCH DRIVER (context for Repair)

You do **not** need the O(1) delta bookkeeping; you will fully re-evaluate. But you must reproduce the **neighborhood enumeration order, acceptance rule, and AFS insertion/deletion semantics**.

### 6.1 `get_chromR` — rebuild routes from linked list
Input `predecessors, successor, routeID` (columns of `Node_related`), and `nb = nb_customer`.
```
chromR = cell(number of distinct routeIDs)
candidates = 1 .. max(max(pred), max(succ))
depots = candidate nodes k with predecessors(k)==0   # route starts
for each depot d:
    walk: x = [d]; repeatedly x_next = successor(x_last); stop when 0
    route = the walked node list (in move/linked-list ids, i.e. customers 1..N, AFS >N)
    chromR{ routeID(d) } = route
remove empty cells
```
Result: `chromR` in **move scheme** (no depot markers; depot is implicit). This is what `Crossover` receives.

### 6.2 `phrase_chromR` — linked list from routes
Input `chromR_move` (routes of ids). For each route, set `predecessors(route(k)) = route(k-1)`, `successor(route(k-1)) = route(k)`, `routeID(route) = j`. First node's predecessor and last node's successor stay 0 (depot). Output arrays sized to total node count.

### 6.3 Entry / mode
- If **not** repair (`isrepair==0`): derive `chromR_move` from `chromR` by subtracting 1 and dropping zeros (detail→move).
- If repair: use `chromR` as-is (already move scheme).
- Build `predecessors, successor, routeID` via `phrase_chromR`; assemble `Node_related = [pred, succ, routeID, node_location']`.

### 6.4 Main search loop
- `whenLastTestedRI` (per customer) and `whenLastModified` (per route) track staleness for the "don't-revisit" pruning (RI = route-improvement memory). `nbMoves` counter, `loopID` counter.
- Outer `while ~searchCompleted`: after the first full pass (`loopID>0`) it sets `searchCompleted=true` so it does **at least two** passes (empty-route moves may only be reachable on the 2nd pass).
- At the top of each pass: for each customer `i`, with prob `1/nbGranular` (`mod(randi, nbGranular)==0`) shuffle its `correlatedVertices(i,:)` row (randomize neighbor order).
- **Time cutoff (non-first individuals):** `if tspid ~= 1 && toc - sol_table.time(end-1) > sol_table.time(end-1)*20/numel(sol_table.ID): break`.
- For each `nodeU = 1..nb_customer`:
  - `correlatedU = correlatedVertices(nodeU,:)`; record `lastTestRINodeU`.
  - There is a large **"patch" bookkeeping block** (lines 60–144) that runs whenever routes were modified since last check (`ispatch != sum(whenLastModified)`): it removes empty routes, re-sorts AFS nodes to the canonical `route+nb_customer` positions, and recomputes PT/PC when a route's AFS was dropped because merging its two segments now fits under `V_Dmax`. For a full-re-eval C++ port, you can **replace this entire block** by: after each accepted move, rebuild `chromR_move` via `get_chromR`, drop degenerate routes, drop an AFS whose two segments re-merge under `V_Dmax`, and recompute `Route_related`/`Penalty_all` from scratch via the §3 formulas. The *observable semantics* are: (a) empty routes are deleted and routeIDs compacted; (b) a route whose only reason for its AFS split was distance may have its AFS removed if the merged segment fits `V_Dmax`, saving `2*dAll(1,nb_customer+2)` of AFS-access distance.
  - For each `nodeV` in `correlatedU` (in order): if the RI-pruning gate passes (`loopID==0 || max(whenLastModified[routeU], whenLastModified[routeV]) > lastTestRINodeU`), try operators **in this fixed order**, taking the **first improving** one (first-improvement, accept if `delta < -1e-6`), then `continue`:
    `m1, m2, m3, m4, m5, m6, m7, m8, m9`, and — **only if `nodeV` is a route start (`predecessors(nodeV)==0`)** — `Depot_m1, Depot_m2, Depot_m3, Depot_m8, Depot_m9`.
  - After exhausting `correlatedU` (when `nodeV == correlatedU(end)` and `loopID != 1`), try **new-route** operators `NewRoute_m1, NewRoute_m2, NewRoute_m3` (open a fresh route for `nodeU`).
- `loopID++` each outer pass.

### 6.5 Acceptance rule (all operators)
Each operator computes `delta = Δdistance + Δ(all weighted penalties)` for the specific move (see m1 example §6.7). **Accept iff `delta < -1e-6`** (strictly improving beyond a tiny tolerance). On accept: apply the linked-list edits, update `Route_related`, `Penalty_all`, `whenLastModified`, `nbMoves`, set `searchCompleted=false`, `isSuccess=true`. On reject: `isSuccess=false`, return unchanged.

**For the C++ full-re-eval port:** for each candidate move, construct the resulting `chromR_move`, run the §3 evaluation to get `cost_Total_new` and feasibility, and accept iff `cost_Total_new < cost_Total_current - 1e-6`. This is behaviorally equivalent to the delta test since `delta = cost_new - cost_old`.

### 6.6 The nine base operators + depot + new-route variants — moves
These are standard HGS/route-improvement moves; each has a plain-move variant plus AFS-aware handling. Node `X = successor(U)`, `Y = successor(V)`.
- **m1**: relocate `U` to just after `V` (insert single node U after V). 
- **m2**: relocate the pair `(U, X)` after `V`.
- **m3**: relocate the pair `(U, X)` after `V` but **reversed** → insert as `(X, U)`.
- **m4**: swap `U` and `V`.
- **m5**: swap pair `(U,X)` with single `V`.
- **m6**: swap pair `(U,X)` with pair `(V,Y)`.
- **m7**: 2-opt within a route — reverse the segment between U and V (`if routeU==routeV`).
- **m8**: 2-opt* between two routes (reconnect `U→Y` and `V→X`, swapping tails).
- **m9**: 2-opt* other reconnection variant (reverse-and-splice across routes).
- **Depot_m1/2/3/8/9**: the same moves but where `V` is a depot/route-start position (inserting at the head of a route / relative to the depot). Triggered only when `predecessors(nodeV)==0`.
- **NewRoute_m1/2/3**: move `U` (or `U,X`) into a **brand-new empty route** (spawn a new vehicle). Triggered once per `nodeU` at the end of its neighbor list.

### 6.7 Conditional AFS Insertion / Deletion (CRI rule) — semantics (from m1)
When a move places customers into a route, the route may or may not need an AFS:
- **Target route already has an AFS** (`max(node_location(route==routeV))==100`): just insert; recompute the affected segment's distance/time; PD/PT/PC recomputed. The pre/post-AFS assignment of the inserted node is inherited from `nodeV`'s `node_location` (`-1` before-AFS, `1` after-AFS).
- **Target route has NO AFS (CRI insertion):** insert `U` after `V` **and append a new AFS node** so the route becomes `[V, U, AFS, Y...]`. Concretely:
  - `costTwo` for this branch = `-dAll(V+1,Y+1) + dAll(V+1,U+1) + dAll(nb_customer+2, Y+1) + dAll(U+1, nb_customer+2)` — i.e. remove edge V→Y, add V→U, add U→AFS, add AFS→Y (AFS table index = `nb_customer+2`).
  - Rebuild the candidate route in move scheme, call `get_pd_pt` to get `d1,d2,t1,t2` (segments split at the new AFS), set `distance_pre/su`, `time_afs/su` accordingly.
  - A **new AFS node** is appended at index `numel(node_location)+1` with `node_location=100`, wired into the linked list between `U` and `Y`, `routeID = routeV`. Node_location of the route is re-marked: nodes up to the AFS → `-1`, from the AFS on → `1`, the AFS itself → `100`.
  - The move is accepted only if the full `delta` (including the extra AFS-access distance) still improves.
- **AFS Deletion (isdelete):** if moving `U` out of `routeU` leaves that route with **only 2 nodes and it contains an AFS** (`sum(routeID==routeU)==2 && node has 100`), the route degenerates to just an AFS visit and must be deleted. `deleteAFS_node` removes the AFS node, compacts `routeID` (all `>delete_idx` decremented) and node indices (`pred/succ > node_deleteAfs` decremented), and the delta credits back `2*dAll(1, nb_customer+2)` (the depot↔AFS round trip that no longer exists). `whenLastModified[routeU]` is removed; `whenLastModified[routeV]` set to `nbMoves`.

**Route bookkeeping after any accepted move:** if a route ends up with a single node (`sum(Node_related(:,3)==routeU)==1`), that route is deleted and higher routeIDs decremented; `Route_related` row removed; `whenLastModified` entry removed.

### 6.8 ELS output (writes back into `sol_table(tspid)`)
```
chromR = get_chromR(Node_related cols)                       # move scheme routes
distance_Total = sum(sum(Route_related(:,[6,7])))            # total distance
cost_Total = distance_Total + sum(Penalty_all(2,:))
IsFeasible = (sum(Penalty_all(2,:)) == 0)
store: chromR_move, distance_Total, cost_Total, predecessors,
       successor, routeID, node_location, per-route Route_related columns,
       penalty_T/C/D/m (= Penalty_all(2,1..4)), IsFeasible, afs_time_delay.
```

---

## 7. `selectparents` — EXACT (target file)

Tournament selection of one parent from the union of both subpopulations.
```
a = (feasiblePop empty)   ? 0 : numel(feasiblePop.ID)
b = (infeasiblePop empty) ? 0 : numel(infeasiblePop.ID)

p1 = randi(a + b)                 # uniform integer in 1..a+b
p1 = (p1 > a) ? infeasiblePop(p1 - a) : feasiblePop(p1)

p2 = randi(a + b)                 # independent uniform
p2 = (p2 > a) ? infeasiblePop(p2 - a) : feasiblePop(p2)

p = (p1.Fitness < p2.Fitness) ? p1 : p2      # lower biased fitness wins; tie → p2
return p
```
Notes: binary tournament over the combined feasible∪infeasible pool. Indices `1..a` map to feasible members (in their sorted order), `a+1..a+b` map to infeasible members. `Fitness` is the biased fitness (§4.4); **strictly-less** comparison means on a tie the second draw (`p2`) is returned. Called twice per main-loop iteration to get `p1` and `p2` for crossover. Requires at least one non-empty subpop (`a+b >= 1`), which holds after the init phase.

---

## 8. `Crossover` — EXACT (target file)

Ordered crossover (OX) on the customer giant-tours of two parents, producing two offspring permutations (only the first, `y1`, is used by the caller).

Signature: `[y1, y2] = Crossover(p1, p2, vrp)`. Inputs are two parent **individuals** (rows with `predecessors{1}`, `successor{1}`, `routeID{1}`). Output: two customer-permutation row vectors of length = number of customers present.

```
# 1. Rebuild each parent's routes (move scheme), flatten to a giant tour,
#    then strip AFS nodes → pure customer sequence.
chromR_p1 = get_chromR(p1.predecessors{1}, p1.successor{1}, p1.routeID{1}, vrp.nb_customer)
chromR_p2 = get_chromR(p2.predecessors{1}, p2.successor{1}, p2.routeID{1}, vrp.nb_customer)
x1 = concatenation of all routes of chromR_p1   (row vector)
x2 = concatenation of all routes of chromR_p2
x1(x1 > nb_customer) = []      # drop AFS nodes, keep customers 1..N only
x2(x2 > nb_customer) = []
nPoint = length(x1)            # == nb_customer (all customers present exactly once)

y1 = zeros(1, nPoint); y2 = zeros(1, nPoint)

# 2. Two cut points (MATLAB: c = randi(nPoint,2,1); point1=min(c); point2=max(c))
c = two independent uniform ints in 1..nPoint
point1 = min(c); point2 = max(c)          # 1 <= point1 <= point2 <= nPoint

# 3. Copy the middle slice directly.
y1(point1:point2) = x1(point1:point2)
y2(point1:point2) = x2(point1:point2)

# 4. Build the donor sequences starting AFTER point2, wrapping around.
x2sorted = [ x2(point2+1 : nPoint), x2(1 : point2) ]
x1sorted = [ x1(point2+1 : nPoint), x1(1 : point2) ]

# 5. Remove from each donor the customers already inherited, then fill the
#    remaining offspring positions (the "outside" of the slice) in donor order.
x2sorted(ismember(x2sorted, y1)) = []     # drop customers already in y1's slice
y1([point2+1:nPoint, 1:point1-1]) = x2sorted   # fill wrap-around positions

x1sorted(ismember(x1sorted, y2)) = []
y2([point2+1:nPoint, 1:point1-1]) = x1sorted
```
Semantics: classic **Order Crossover (OX1)**. `y1` inherits `x1`'s middle segment `[point1..point2]`; the rest of `y1` is filled with `x2`'s customers (starting just after `point2`, wrapping) in their `x2` order, skipping those already inherited. `y2` is symmetric (inherits `x2`'s middle, fills from `x1`). AFS nodes are discarded before crossover — offspring is a pure customer permutation that is subsequently re-split into routes.

**C++ indexing:** convert the MATLAB 1-based ranges carefully. With 0-based arrays: `point1, point2 ∈ [0, nPoint-1]`, `point1 = min`, `point2 = max`; inherited slice is indices `[point1..point2]` inclusive; donor sequence is `x2[point2+1 .. nPoint-1] ++ x2[0 .. point2]`; fill positions are `[point2+1 .. nPoint-1] ++ [0 .. point1-1]`. Filtering "already in slice" uses membership over the inherited customer set. Guard the degenerate case `point1==point2` (single-element slice) — the fill still covers all remaining positions.

**Caller usage (Main line 129):** `[offspring_tsp, ~] = Crossover(p1,p2,vrp)` — only `y1` is kept as the offspring giant-tour.

---

## 9. IMPLEMENTATION CHECKLIST / EDGE CASES

- **Depot is node index 1 (detail) / value 0 (move)**; customer `c` is `c+1` (detail) / `c` (move); AFS visits are `> nb_customer(+1)`; the shared AFS distance-table index is `nb_customer+2`.
- **Feasibility** = all four raw penalties zero. `wM=0` makes PM inert but still computed.
- **Cost** = total distance + weighted (PT+PC+PD+PM).
- **Acceptance** in ELS = strict improvement by `> 1e-6` (delta `< -1e-6`).
- **Best-sol** update tolerance = `1e-9`; only feasible solutions qualify.
- **Repair** = ELS re-run with each active penalty weight ×10 (`WP=10`), local `Penalty_all` copy; re-insert into pop + update best only if it becomes feasible.
- **selectparents** = binary tournament over feasible∪infeasible by biased fitness (lower wins, tie→second draw).
- **Crossover** = OX on customer permutations; AFS stripped; only `y1` used.
- **Penalty adaptation** every `nbLast=20` iters from `tspid>=100`, per-constraint, using `targetFeasible=0.2 ± 0.05`, scale `1.2` (cap `1e5`) / decrease `0.85` (floor `0.1`); rescales stored infeasible penalties and re-sorts infeasible pop.
- **Survivor selection** trims to `popSizeMu=154` when a subpop exceeds `μ+λ=222`, deleting clones first, then worst biased fitness; index-1 (best-cost) is protected.
- **Known source quirks to reproduce**: main loop splits the stale `tsp` rather than `offspring_tsp`; `disp('...PD...', PC)` typo in a printout (cosmetic); `i=1;ii=1` reassignments in the AFS-delay loops are MATLAB no-ops (don't reset loop counters).
- **`Repair_sol` return `isNewBest`** is only set true when repair produces a feasible solution strictly better than `bestSolRestart`.

Relevant source paths (absolute):
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Crossover selectparents/Crossover.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Crossover selectparents/selectparents.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/Repair_sol.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/ELS_mian.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/get_chromR.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/phrase_chromR.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/get_pd_pt.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/get_pc_now.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/AFSdelay_new.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/AFSdelay_recursion.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/deleteAFS_node.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/m1.m (representative operator with CRI AFS insertion/deletion)
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Load/chromR_detail_all.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Load/get_vrp.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Main_METS.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Population management/PopManagement.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Population management/add2Pop.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Population management/infeasiblePop_updateBiasedFitnesses.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Population management/uti_updateBestSol.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Population management/uti_addSol2Last100.m

Remaining operators not individually transcribed (m2–m9, Depot_m*, NewRoute_m*, split_Dmax, split_Tmax, deleteAFS_delta): their move semantics are described in §6.6–§6.7; if the C++ port needs the exact edge-arithmetic of each (e.g. the precise `costOne/costTwo` edge sets for m4–m9 or the exact Split partition rule), those specific files should be read next — they were not required to specify the three target files but are needed for a complete ELS/Split port.