I now have everything needed to write the exhaustive spec for the 8 population-management files, with all cross-referenced data structures (solution-table schema, node encoding, predecessor/successor/routeID representation, penalty formulas, and control flow). The task is scoped to the population-management module, so I'll focus the spec there while fully documenting the data structures and formulas those files depend on. I don't need to read the local-search operator internals (m1-m9) since section 4's operators are not in the target file set — but I'll note the semantics the population files rely on.

Let me write the specification.

# METS — Population Management Module: Exhaustive C++ Implementation Specification

This spec covers the 8 files in `Population management/` plus every data structure, constant, and formula they depend on. It is self-contained: a C++ programmer can implement the module from this alone. All array indices below are described in **MATLAB 1-based convention**; conversion notes to 0-based C++ are given where it matters.

---

## 0. INDEXING & NODE ENCODING (READ FIRST — load-bearing)

There are TWO different node-numbering conventions in this codebase. Getting them wrong breaks the broken-pair distance.

### 0.1 `chromR` convention (route arrays, "with depot")
- Node **1** = depot.
- Node **c+1** = customer `c` (customers are `1..nb_customer`, stored as `2..nb_customer+1`).
- Nodes **> nb_customer + 1** = AFS (refueling stations), appended after customers.
- A route is a MATLAB row/col vector like `[1, 5, 3, 1]` (starts and ends at depot = 1). `chromR` is a cell array of such vectors, one per vehicle/route.
- `chromR_move{i}` = same routes but with depot removed and all indices shifted down by 1: `a = chromR{i} - 1; a(a==0) = [];`. So in `_move` convention: customer `c` = index `c`, AFS = index `> nb_customer`, depot does not appear.

### 0.2 Linked-list convention (`predecessors`, `successor`, `routeID`, `node_location`)
Built by `phrase_chromR` from `chromR_move`. These are **column vectors indexed by node id in `_move` convention** (customer `c` at position `c`; AFS at positions `nb_customer+1 ...`). Length = total number of distinct non-depot nodes across all routes (customers + all AFS visits).

For node `k` (1-based position = its `_move` id):
- `predecessors(k)` = the node that immediately precedes `k` in its route, in `_move` ids; **`0` means the predecessor is the depot** (k is first customer/AFS in its route).
- `successor(k)` = node immediately following `k`; **`0` means the successor is the depot** (k is last in its route).
- `routeID(k)` = 1-based index of the route/vehicle containing `k`.
- `node_location(k)` (also called `I`): position marker relative to the route's AFS:
  - `-1` = node is **before** the AFS in its route, OR the route has no AFS at all.
  - `1` = node is **after** the AFS.
  - `100` = node **is** the AFS (the marker sentinel for "this node is a refueling station"). Exactly one node per route that has an AFS carries `100`.

`phrase_chromR(chromR_move)` pseudocode:
```
a = total count of nodes over all routes in chromR_move
predecessors = successor = routeID = zeros(a,1)
for j = 1..num_routes:
    route = chromR_move{j}          # e.g. [5, 3, 12] (12 is an AFS if >nb_customer)
    for k = 2..len(route):
        predecessors[ route[k]   ] = route[k-1]
        successor[   route[k-1] ] = route[k]
    routeID[ route[all] ] = j
# positions never written keep 0 (first node's pred = 0, last node's succ = 0)
```

**These four column vectors, plus `chromR`/`chromR_move`, are exactly what the population files read.** The population files NEVER recompute routes; they treat `predecessors`/`successor`/`node_location` as opaque per-solution vectors.

---

## 1. DATA STRUCTURES

### 1.1 `vrp` (problem instance; read-only inside population mgmt except the two distance caches)
Fields used by the 8 target files:
- `vrp.nb_customer` (int): number of customers `n`.
- `vrp.last_customer` (int): **the divisor used to normalize broken-pair distance.** Treat as the customer count for normalization (equal to `nb_customer` for these instances). Store and use exactly as given by the instance — do NOT substitute `nb_customer` blindly if the instance sets it differently.
- `vrp.ALL_brokenDIS_feasible` (double matrix): cached pairwise broken-pair distance matrix for the FEASIBLE population. Persisted across calls in `vrp`. Preallocated width `popSizeLambda + popSizeMu + 1`.
- `vrp.ALL_brokenDIS` (double matrix): same, for the INFEASIBLE population.

(Other `vrp` fields — `V_speed`, `T_Customer`, `T_max_V`, `T_Afs`, `C_Afs`, `V_Dmax`, `V_nb`, `distance_table`, `correlatedVertices` — are used by cost evaluation / local search, documented in §3 and §6, not by the population files directly.)

### 1.2 `par_hgs` (algorithm parameters)
Computed in `Main_METS`; constants from the parameter block:
- `popSizeMu = 154` — target population size μ (the size a subpopulation is trimmed back down to).
- `popSizeLambda = 68` — offspring/overflow allowance λ. A subpopulation is allowed to grow to `popSizeMu + popSizeLambda` before survivor selection triggers.
- `el = 0.5`; `eliteNum = floor(el * popSizeMu) = floor(0.5*154) = 77`.
- `nc = 0.2`; `nClosest = floor(nc * popSizeMu) = floor(0.2*154) = 30`. Number of nearest neighbors used in the diversity contribution.
- `nbGranular = 20` (local search granularity; not used by pop mgmt).
- `targetFeasible = 0.2` (penalty adaptation target).
- `nbLast = 20` (size of the sliding "Last100" window — despite the name it is 20).
- `Penalty_all` (2×4 double): row 1 = weights `[wT, wC, wD, wM]`, row 2 = weighted penalty values `[wT*penalty_T, wC*penalty_C, wD*penalty_D, wM*penalty_m]`. Initial weights `[527, 195, 430, 0]`.

### 1.3 The solution "table row" (`sol_table` row) — one individual
In MATLAB this is a `table` where each column may be a scalar or a cell (holding an array). In C++, model as a `struct Solution`. Fields consumed/produced by the population files:

Per-individual identity & cost:
- `ID` (int): unique, monotonically increasing = the `tspid` iteration index when created. Larger ID = newer solution. Used as a tie-breaker.
- `cost_Total` (double): penalized objective (see §3). This is the sort key for cost.
- `distance_Total` (double): pure total travelled distance (no penalties).
- `IsFeasible` (0/1): 1 iff all penalties are 0.
- `penalty_T`, `penalty_C`, `penalty_D`, `penalty_m` (double): **weighted** penalty values (`w* * raw_penalty`). Note `penalty_*==0` is the per-constraint feasibility test used in penalty adaptation.

Route representation (all in the conventions of §0):
- `chromR` (cell of vectors): routes with depot, `chromR` convention.
- `chromR_move` (cell of vectors): routes without depot, `_move` convention.
- `tsp`, `tsp_now` (giant-tour permutations): customer order; used by crossover, not pop mgmt.
- `predecessors` (col vector): §0.2.
- `successor` (col vector): §0.2.
- `routeID` (col vector): §0.2.
- `node_location` (ROW vector): §0.2. **Note the transpose**: stored as a row (`Node_related(:,4)'`). `sum(node_location == 100)` = number of AFS visits in the solution.

Diversity / fitness bookkeeping (managed BY the population files):
- `Fitness` (double): **biased fitness** — LOWER is better. Formula §3.4.
- `fitRank` (double in `[0,1]`): normalized rank by cost. `0` = cheapest.
- `divRank` (double in `[0,1]`): normalized rank by diversity contribution. `0` = most diverse (largest avg distance to neighbors, since `avgBrokenDist` is stored negated).
- `avgBrokenDist` (double): **negated** mean distance to the `nClosest` nearest population members. More negative = more diverse/isolated. (Negation makes "sort ascending → most diverse first".)
- `brokenPairDistance` (cell holding a row vector) [add2Pop/PopManagement path]: this solution's row of the cached distance matrix.
- `brokenDist` (cell holding a row vector) [update_*/infeasiblePop_updateBiasedFitnesses path]: this solution's distances to all OTHERS (self removed). Two different field names for essentially the same per-row data — see §2 notes.
- `time` (double): wall-clock timestamp; carried along, not used in pop logic.

### 1.4 The two subpopulations
- `feasiblePop`: array (MATLAB table; C++ `std::vector<Solution>`) of feasible individuals, **kept sorted by `cost_Total` ascending** (position 1 = best cost). `add2Pop` maintains this sorted-insert invariant.
- `infeasiblePop`: same, for infeasible individuals.
- Both are capped: after `PopManagement`, size ≤ `popSizeMu`. They may transiently reach `popSizeMu + popSizeLambda + 1` before trimming.

### 1.5 The distance caches (critical, easy to get wrong)
`vrp.ALL_brokenDIS_feasible` and `vrp.ALL_brokenDIS` are **square symmetric matrices** (padded with trailing zero rows/cols up to `popSizeLambda+popSizeMu+1`). Entry `(i,j)` = broken-pair distance between population member at position `i` and position `j` (positions match the current order of the corresponding subpopulation). Diagonal = 0. These are incrementally maintained by `add2Pop`/`PopManagement` so the full O(N²) recompute is avoided.

**C++ porting note (per your instructions): you may drop the incremental cache bookkeeping and simply RE-EVALUATE the full pairwise matrix each time from `predecessors`/`successor`/`node_location`.** The semantics of the matrix content are fully specified in §2.1 (the broken-pair distance formula). Only the *content* matters for correctness; the shifting/insertion gymnastics in `add2Pop`/`PopManagement` are pure performance plumbing to keep that matrix consistent under insert/delete.

---

## 2. THE BROKEN-PAIR DISTANCE (the one formula shared by all files)

This is THE core diversity primitive. Every population file computes it identically (up to the two conventions of clipping). Implement once.

### 2.1 `brokenPairDistance(solA, solB)` → double in `[0, ~2]`

Inputs: two solutions' `successor` (`sA`, `sB`) and `predecessors` (`pA`, `pB`) vectors, `vrp.nb_customer`, `vrp.last_customer`.

```
a = min(length(sA), length(sB))          # clip to shorter (AFS counts can differ)
sA = sA[1..a]; sB = sB[1..a]; pA = pA[1..a]; pB = pB[1..a]

# Collapse ALL AFS ids to a single canonical id so different AFS choices
# do not count as structural differences:
for each vector V in {sA, pA, sB, pB}:
    V[V > nb_customer] = nb_customer + 1

differences = 0
for jj = 1..a:
    # (1) "extra points": successor of jj differs, and is not merely the reversed edge
    if sA[jj] != sB[jj]  AND  sA[jj] != pB[jj]:
        differences += 1
    # (2) "missing points": jj starts a route in A (pred=depot=0) but is mid-route in B
    if pA[jj] == 0  AND  pB[jj] != 0  AND  sB[jj] != 0:
        differences += 1

distance = differences / vrp.last_customer
return distance
```

Notes / tie-breaks / edge cases:
- The `sA[jj] != pB[jj]` clause makes an edge count as "unbroken" even if traversed in the opposite direction in B. This is the standard broken-pairs symmetry allowance.
- Clipping to `min` length means the value can be mildly **asymmetric** depending on which solution is longer; the caller (`add2Pop`) deliberately orders the two so the SHORTER-derived value is stored in both `(i,j)` and `(j,i)` — see §5.1. In your full-recompute port, pick a fixed convention (e.g. always clip to min, always put the smaller-index solution as "A") and store it symmetrically. The precise A/B ordering rule from `add2Pop`: when comparing new solution `new_solnum` with member `ii`, **if `ii < new_solnum`, then A = member `ii`, B = new; if `ii > new_solnum`, then A = new, B = member `ii`.** (i.e. A = the lower position.)
- Divisor is `vrp.last_customer` (NOT `a`, NOT `nb_customer` literally — use the field).
- In `update_feasiblePop`/`update_infeasiblePop` the loop upper bound is written differently but is equivalent: `jj = 1 .. min(#AFS_in_p, #AFS_in_q) + nbClients`, where `#AFS = sum(node_location == 100)`. That equals `min(len(sA),len(sB))` = `a`. Those two functions **do NOT collapse AFS ids** and do **not** clip individual vectors — they just bound the loop; but since customers occupy the first `nbClients` positions and AFS positions beyond are compared directly, results match in practice for equal-length solutions. **For the C++ port, use the single canonical formula in §2.1 everywhere** (the collapse-to-`nb_customer+1` version) — it is the most-used and most-correct variant (`add2Pop`, `infeasiblePop_updateBiasedFitnesses`).

### 2.2 `avgBrokenDist` (diversity contribution of one member)
Given the full pairwise matrix `M` (N×N, diagonal 0):
```
maxSize = min(nClosest, N-1)
for each member c:
    row = M[c, :] with entry c removed          # distances to all OTHERS
    kSmallest = the maxSize SMALLEST values of row   # mink(row, maxSize)
    avgBrokenDist[c] = -mean(kSmallest)          # NEGATED
```
- `mink(x, k)` = the `k` smallest elements (MATLAB). Nearest neighbors = smallest broken-pair distances.
- Negation: most diverse (large distances → large mean → most negative) sorts first when sorting ascending.
- Edge case `N == 2`: `mink` on a 1-element column returns a scalar; MATLAB then duplicates it to length 2 (`avgBrokenDist(1:2) = avgBrokenDist`). In C++ just compute per-member normally; the duplication is a MATLAB reshape artifact.
- Edge case `N == 1`: no diversity; avg undefined → treated as 0 (single-member init path sets `avgBrokenDist = 0`, `divRank = 0`).

---

## 3. OBJECTIVE & PENALTIES (from `chromR_detail_all` / `ELS_mian`)

The population files consume `cost_Total`, `distance_Total`, `IsFeasible`, and the four weighted `penalty_*`. Here is exactly how those are formed, so the port is faithful.

### 3.1 Weights and adaptation state
`Penalty_all(1,:) = [wT, wC, wD, wM]`, initially `[527, 195, 430, 0]`. `wM = 0` throughout (vehicle-count penalty is disabled by weight). Weights adapt in the main loop (§5.7).

### 3.2 Raw penalties (per solution)
- **PM (vehicle count):** `penalty_m = max(nb_V - vrp.V_nb, 0)` where `nb_V = numel(chromR)` (number of routes). Weighted: `wM*penalty_m` (= 0 since wM=0).
- **PD (distance / Dmax):** For each route split into segments by its AFS (segment before AFS, segment after AFS; if no AFS the whole route is one segment): compute segment distance `D`; `penalty_D_v(route,seg) = max(D - vrp.V_Dmax, 0)`. `penalty_D = sum over all routes,segments`. Weighted `wD*penalty_D`.
- **PT (time / Tmax):** For each route, `time_V = distance_V / V_speed + (num_nodes_in_route - 2) * T_Customer` (service time per non-depot stop). `overtime_V = max(time_V - vrp.T_max_V, 0)`. `penalty_T = sum(overtime_V)`. Weighted `wT*penalty_T`. (Per-route `penalty_T_v = max(time_V - T_max_V, 0)`.)
- **PC (AFS shared-capacity / refuel congestion):** The complex part. Each vehicle's AFS arrival time `afs_time` is computed; a conflict-resolution pass (`chromR_detail_all` lines 136-221, and `AFSdelay_new` in ELS) shifts arrival times within each vehicle's slack (`time_V_shifting = max(T_max_V - time_V, 0)`) to deblock the shared station of capacity `C_Afs` with per-refuel duration `T_Afs`. After scheduling, for each distinct time boundary the number of vehicles simultaneously refueling is counted; `nb_fueling = max(count - C_Afs, 0)` per interval; `penalty_C = max( nb_fueling · time_during', 0 )` (dot product of overflow counts and interval durations). Weighted `wC*penalty_C`.

**Porting note:** you said full re-evaluation is fine and you do NOT need the O(1) deltas. For PC you DO still need to reproduce the scheduling/conflict semantics (it is not a delta — it is the definition of the penalty). Treat `chromR_detail_all` lines 65-279 + `AFSdelay_new` as the reference algorithm for PC. This is outside the 8 target files but is the definition of a field they read.

### 3.3 `cost_Total` and `IsFeasible`
```
cost_Total   = wT*penalty_T + wC*penalty_C + wD*penalty_D + wM*penalty_m + distance_Total
distance_Total = sum of all route distances (pure, unpenalized)
IsFeasible = 1  iff  (penalty_T + penalty_C + penalty_D + penalty_m) == 0   (raw, i.e. all zero)
             else 0
```
(In `ELS_mian` the equivalent test is `sum(Penalty_all(2,:)) != 0 → infeasible`.)

### 3.4 Biased fitness (the survivor/selection key) — LOWER IS BETTER
For a subpopulation of size `N`:
```
Fitness[i] = fitRank[i] + (1 - eliteNum/N) * divRank[i]
```
where:
- `fitRank[i] ∈ [0,1]`: normalized cost rank. Sort members by `cost_Total` ascending, assign raw rank `1..N`, then `fitRank = (rank - 1)/(N - 1)`. Best cost → 0.
- `divRank[i] ∈ [0,1]`: normalized diversity rank. Sort by `avgBrokenDist` ascending (recall it's negated, so most-diverse first), assign raw rank `1..N`, then `divRank = (rank - 1)/(N - 1)`. Most diverse → 0.
- `eliteNum = 77`. The factor `(1 - eliteNum/N)` shrinks the diversity weight as the population grows; when `N ≤ eliteNum` the factor can go ≤ 0 (diversity actively *rewarded* into negative territory for elites — this is intended HGS behavior).

`Fitness` drives: (a) tournament parent selection (lower wins), (b) survivor selection (highest = worst = removed).

---

## 4. FUNCTION-BY-FUNCTION SPEC

### 4.1 `add2Pop(sol_table_row, par_hgs, tspid, feasiblePop, infeasiblePop, vrp) → [feasiblePop, infeasiblePop, vrp]`
Purpose: insert ONE new solution into its subpopulation (feasible or infeasible per `IsFeasible`), keep it cost-sorted, update the cached distance matrix, recompute all ranks/fitness. The feasible and infeasible branches are byte-for-byte identical logic on different data (`feasiblePop`+`ALL_brokenDIS_feasible` vs `infeasiblePop`+`ALL_brokenDIS`). Implement once, parameterize by `b = IsFeasible`.

Let `POP` = target subpopulation, `M` = its cached distance matrix, `L = popSizeLambda + popSizeMu + 1` (preallocation width = `68+154+1 = 223`).

```
if POP is empty:                              # first member
    POP = [ sol ]
    sol.Fitness   = 1
    sol.fitRank   = 1
    sol.brokenPairDistance = zeros(L,1)
    sol.avgBrokenDist = 0
    sol.divRank   = 0
    M = zeros(1, L)                            # single zero row
    return

else:
    N_old = size(POP)
    M = M[1..N_old, :]                          # trim cached matrix to current size

    # --- sorted insertion by cost_Total (ascending) ---
    append sol at end of POP; sol.fitRank = 99 (sentinel)
    new_solnum = size(POP)                       # provisional position = last
    M add a new zero row at index new_solnum

    if sol.cost_Total < POP[end-1].cost_Total:   # not already the worst → find slot
        for ii = 1..N_new:                        # scan from best
            if sol.cost_Total <= POP[ii].cost_Total:
                # insert sol at position ii, shift [ii..] down by one
                move sol into position ii (remove the temporary last copy first)
                new_solnum = ii
                # mirror the same insert-shift on matrix M: insert a zero
                #   row AND zero column at index new_solnum
                M: delete last row; shift rows [new_solnum..] down; 
                   shift cols [new_solnum..] right; set row/col new_solnum = 0
                break
        # (if sol is <= nobody it stays at end; new_solnum stays = last)
    else:
        # sol is the worst → stays at end; just init its brokenPairDistance cell

    # --- recompute cost ranks for whole pop ---
    fitRank[i] = i        for i=1..N            # raw ascending (pop already cost-sorted)
    fitRank    = (fitRank - 1)/(N - 1)          # normalize to [0,1]

    # --- update distance matrix: only the NEW row/col need real values ---
    for ii = 1..N, ii != new_solnum:
        # order A/B so the SHORTER-derived value is stored (see §2.1 ordering rule):
        if ii < new_solnum:  A = POP[ii],  B = sol
        else (ii>new_solnum):A = sol,      B = POP[ii]
        d = brokenPairDistance(A, B)             # §2.1
        M[new_solnum, ii] = d;  M[ii, new_solnum] = d

    save M back to vrp cache

    # --- recompute avgBrokenDist for ALL members from M ---
    maxSize = min(nClosest, N - 1)
    for c = 1..N:
        row = M[c,1..N] with entry c removed
        avgBrokenDist[c] = -mean( mink(row, maxSize) )
    if N == 2: duplicate the single value (MATLAB artifact; see §2.2)

    store per-member brokenPairDistance{c} = M[c, 1..N]   (full row incl. self=0)

    # --- diversity rank + biased fitness ---
    sortedIndex = argsort(avgBrokenDist ascending)
    divRank[sortedIndex] = 1..N                   # raw
    divRank = (divRank - 1)/(N - 1)               # normalize
    Fitness = fitRank + (1 - eliteNum/N) * divRank
```

Key invariants for the port:
- After `add2Pop`, `POP` is sorted by `cost_Total` ascending and `M` is the consistent N×N (padded to L) symmetric distance matrix in the SAME order.
- `fitRank`, `divRank`, `avgBrokenDist`, `Fitness` are fully populated for every member.
- **Full-recompute simplification:** the entire "shift rows/cols, only fill new row" machinery can be replaced by: insert sol at correct cost-sorted position, then recompute the full N×N `M` from scratch via §2.1, then §2.2 + fitness. Result is identical.

### 4.2 `PopManagement(sol_table_row, par_hgs, tspid, feasiblePop, infeasiblePop, vrp) → [feasiblePop, infeasiblePop, vrp]`
Purpose: call `add2Pop`, then if the subpopulation overflowed (`> popSizeMu + popSizeLambda`), run **survivor selection** trimming it back to `popSizeMu`, removing clones first then worst-fitness. Finally recompute ranks/fitness.

```
[feasiblePop, infeasiblePop, vrp] = add2Pop(sol, ...)

b = sol.IsFeasible
subpop = (b ? feasiblePop : infeasiblePop)
M      = (b ? vrp.ALL_brokenDIS_feasible : vrp.ALL_brokenDIS)

if size(subpop) > popSizeMu + popSizeLambda:      # overflow → trim
    IsUpdateFitness = 1
    del_ids = []                                   # IDs removed (for matrix col cleanup)
    while size(subpop) > popSizeMu:                # remove one per iteration
        isWorstClone = false
        worstFitness = -999999999                  # larger = worse
        worstPos = (none)
        for i = 2 .. size(subpop):                 # NOTE: starts at 2 → never removes best-cost member (position 1)
            # clone test: distance to nearest OTHER member ~ 0
            twoSmallest = mink(M[i, :], 2)         # smallest two incl. self(=0)
            isClone = (sum(twoSmallest) < 1e-8)    # self is 0, so this tests nearest-other ≈ 0
            if (isClone && !isWorstClone) OR
               (isClone == isWorstClone && subpop[i].Fitness > worstFitness):
                worstFitness = subpop[i].Fitness
                isWorstClone = isClone
                worstPos     = i
        del_ids.append( subpop[worstPos].ID )
        remove subpop[worstPos]                     # also remove row worstPos from M
    # remove the corresponding COLUMNS from M (matched by original ID → column index),
    # then pad zero columns back at the end to keep width
    save M back to vrp cache

else:
    IsUpdateFitness = 0

if IsUpdateFitness:                                 # ranks changed → recompute
    fitRank[i] = i ; fitRank = (fitRank-1)/(N-1)     # pop still cost-sorted
    maxSize = min(nClosest, N-1)
    for c: avgBrokenDist[c] = -mean(mink(M[c,:] without self, maxSize))
    sortedIndex = argsort(avgBrokenDist asc)
    divRank[sortedIndex] = 1..N ; divRank = (divRank-1)/(N-1)
    Fitness = fitRank + (1 - eliteNum/N)*divRank

write subpop back to feasiblePop or infeasiblePop
```

Survivor-selection semantics (the 3 cases the comments enumerate — reproduce EXACTLY):
- **Clone** = a member whose nearest *other* member is at distance ≈ 0 (`< 1e-8`). Detected via `sum(mink(M[i,:], 2)) < 1e-8` (member's own diagonal is 0, so the two smallest are `{0, nearest-other}`; their sum ≈ 0 means nearest-other ≈ 0).
- Selection of the "worst" (to delete), scanning `i = 2..N`:
  - Case 1 (no clone yet, current not clone): pick by worst (highest) `Fitness`.
  - Case 2 (current IS clone, no clone recorded yet): clone always beats non-clone → take it regardless of fitness.
  - Case 3 (both current and recorded are clones): tie-break by worst `Fitness`.
  - Encoded by: `(isClone && !isWorstClone) || (isClone==isWorstClone && Fitness[i] > worstFitness)`.
- **Position 1 (lowest cost) is never eligible for removal** (loop starts at `i=2`). The single best-cost solution is protected.
- Remove exactly one member per while-iteration; repeat until size == `popSizeMu`.

**Full-recompute simplification:** you may (a) run the clone/worst selection loop against a freshly recomputed `M`, and (b) after all removals, recompute the whole N×N `M` + ranks + fitness once. The per-removal matrix row/col surgery is only there to keep the cache valid mid-loop.

### 4.3 `update_feasiblePop(feasiblePop, nbClients, eliteNum, nClosest) → feasiblePop`
Purpose: FULL from-scratch recomputation of the feasible population's ranks, distances, and fitness, then re-sort by `Fitness`. Used on restart / bulk refresh (not the incremental hot path).

```
feasiblePop = sort by cost_Total ascending
fitRank = (1..N)'; fitRank = (fitRank-1)/(N-1)
N = size(feasiblePop)

# full O(N^2) pairwise matrix
brokenPairDistance = zeros(N, N)      # (allocated (N-1)x(N-1) in code but filled as NxN symmetric)
for pp = 1..N-1:
  for qq = pp+1..N:
     # loop bound = min(#AFS_pp, #AFS_qq) + nbClients  == min(len_pp, len_qq)
     differences = 0
     for jj = 1 .. min( sum(node_location{pp}==100), sum(node_location{qq}==100) ) + nbClients:
        if successors_pp[jj] != successors_qq[jj] AND successors_pp[jj] != predecessors_qq[jj]:
            differences++
        if predecessors_pp[jj]==0 AND predecessors_qq[jj]!=0 AND successors_qq[jj]!=0:
            differences++
     d = differences / nbClients
     brokenPairDistance[pp,qq] = brokenPairDistance[qq,pp] = d

if matrix nonempty:
  for c = 1..N:
     row = brokenPairDistance[c,:] with entry c removed   # brokenDist{c}
     maxSize = min(nClosest, N-1)
     avgBrokenDist[c] = -mean( mink(row, maxSize) )
  feasiblePop = sort by avgBrokenDist ascending
  divRank[cc] = cc for cc=1..N ; divRank = (divRank-1)/(N-1)
  Fitness[ff] = fitRank[ff] + (1 - eliteNum/N)*divRank[ff]
  feasiblePop = sort by Fitness ascending          # FINAL order = biased fitness
```
Differences vs the incremental path: (a) divisor is `nbClients` (not `vrp.last_customer`) and (b) AFS ids are NOT collapsed to a single sentinel and vectors are NOT clipped to min-length individually — only the loop bound clips. **For the port, keep §2.1 as the single distance function; the discrepancy only manifests when two compared solutions have different AFS multisets, and the intended behavior is the collapse version.** Final sort here is by **Fitness ascending** (unlike `add2Pop`/`PopManagement`, which keep the pop cost-sorted). Callers of `update_*` must not assume cost order afterward.

### 4.4 `update_infeasiblePop(infeasiblePop, nbClients, eliteNum, nClosest) → infeasiblePop`
**Identical** to `update_feasiblePop` with `feasiblePop → infeasiblePop`. Same pseudocode. Implement as one templated/parameterized function.

### 4.5 `infeasiblePop_updateBiasedFitnesses(infeasiblePop, vrp, par_hgs) → infeasiblePop`
Purpose: after penalty weights change (§5.7), the infeasible pop's `cost_Total` values shifted, so re-rank and re-sort. Same structure as `update_infeasiblePop` but with: the canonical §2.1 distance (collapse AFS, clip to min, divide by `vrp.last_customer`), stable tie-break by ID, and `linspace` divRank.

```
# cost sort with newest-first tie-break:
infeasiblePop.ID = -infeasiblePop.ID
infeasiblePop = sortrows by {cost_Total ASC, ID ASC}   # ID negated → larger original ID first on ties
infeasiblePop.ID = -infeasiblePop.ID
fitRank = (1..N)'                                        # raw (normalized later)
N = size

brokenPairDistance = zeros(N,N)
for pp=1..N-1, qq=pp+1..N:
    brokenPairDistance[pp,qq] = brokenPairDistance[qq,pp] = calculateBrokenPairDistance(pp,qq)
        # == §2.1 EXACTLY: clip both to a=min(len), collapse ids >nb_customer to nb_customer+1,
        #    count cc + dd, divide by vrp.last_customer

if any(brokenPairDistance != 0):
    for c=1..N:
        row = brokenPairDistance[c,:] with c removed     # brokenDist{c}
        maxSize = min(nClosest, N-1)
        avgBrokenDist[c] = -mean(mink(row, maxSize))
    infeasiblePop = sort by avgBrokenDist ascending
    divRank = linspace(0, 1, N)'                          # 0, 1/(N-1), ..., 1
    fitRank = (fitRank - 1)/(N - 1)                        # normalize (note: fitRank was NOT re-derived after the avgBrokenDist sort — it carries the cost-sorted raw ranks; see caveat)
    Fitness[ff] = fitRank[ff] + (1 - eliteNum/N)*divRank[ff]
    # final sort by fitness, newest-first tie-break:
    infeasiblePop.ID = -infeasiblePop.ID
    infeasiblePop = sortrows by {Fitness ASC, ID ASC}
    infeasiblePop.ID = -infeasiblePop.ID
```
CAVEAT to reproduce faithfully: `fitRank` is assigned `(1..N)` in the **cost-sorted** order, but is only **normalized** AFTER the population has been re-sorted by `avgBrokenDist`. Because MATLAB `sortrows` reorders the whole table (including the already-assigned `fitRank` column), each row keeps its cost-based raw rank value through the diversity re-sort — so `fitRank` correctly stays tied to the row's cost position, then gets normalized in place. Net effect: `fitRank[row] = (costRank(row) - 1)/(N-1)`, `divRank[row] = (divRank position)/(N-1)`. Implement by: compute both integer ranks per solution, then normalize; do NOT re-derive fitRank from post-diversity order.

Tie-break rule (used 3×): on equal cost or equal Fitness, the solution with the **larger ID (newer)** comes first. Implement as comparator `(key ASC, ID DESC)`.

### 4.6 `uti_updateBestSol(bestSolRestart, bestSolOverall, sol) → [bestSolRestart, bestSolOverall, isNewBest]`
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
- Only **feasible** solutions can become best.
- Strict improvement by margin `1e-9`.
- `bestSolRestart` = best since last restart; `bestSolOverall` = global best. Overall only updated when restart-best is (and it also strictly improves overall).

### 4.7 `uti_addSol2Last100(Last100, individual, par_hgs) → Last100`
Sliding window (FIFO) of the last `nbLast (=20)` solutions.
```
if Last100 empty:
    Last100 = [individual]
else if size(Last100) < nbLast:
    append individual
else:
    drop the oldest (front), append individual   # Last100[1..end-1] = Last100[2..end]; Last100[end]=individual
```
Used for penalty adaptation (fraction of recent solutions satisfying each constraint). Only `penalty_T`, `penalty_C`, `penalty_D`, `ID` fields are actually read downstream, but store full rows.

### 4.8 `bestsolfind(best_sol, feasiblePop, instance_id) → best_sol`
```
b = argmin(feasiblePop.cost_Total)      # index of cheapest feasible member
best_sol[instance_id] = feasiblePop[b]
```
Trivial: pick the minimum-cost feasible individual into a results table row. (Utility, called at end/reporting.)

---

## 5. CONTROL FLOW (driver context for the module)

From `Main_METS.m`. Provided so the port wires the population functions correctly. Constants: `split_prob=0.5`, `PT=527, PC=195, PD=430`, `penaltyScaleFactor=1.2`, `penaltyDecreaseFactor=0.85`, `popSizeMu=154`, `popSizeLambda=68`, `targetFeasible=0.2`, `nbLast=20`, `maxIterNonProd=300`, `maxIter=2000`, `timeLimit=100000` (s).

### 5.1 Initial population construction
```
rng(SEED + 1)
tsp_all = popSizeMu*4 (=616) random customer permutations (randperm(nbClients))
for i = 1 .. popSizeMu*4:
    tspid = i
    if elapsed > timeLimit: break
    if tspid > maxIter:     break
    tsp = tsp_all[i]
    if rand <= split_prob: chromR = split_Dmax(vrp, tsp, par_hgs)     # 50%
    else:                  chromR = split_Tmax(vrp, tsp, par_hgs)     # 50%
    sol = chromR_detail_all(...)          # evaluate cost/penalties (§3)
    sol = ELS_mian(...)                    # local search (§6) → sets predecessors/successor/etc
    sol.time = elapsed
    [feasiblePop, infeasiblePop, vrp] = PopManagement(sol, ...)
    [bestSolRestart, bestSolOverall] = uti_updateBestSol(...)
    Last100Sol = uti_addSol2Last100(...)
    # repair with probability 0.5 if infeasible:
    if sol.IsFeasible == 0 AND (rand - 0.5 > 1e-6):     # ≈ 50% chance
        Repair_sol(...)   # attempts to make feasible, re-runs pop mgmt/best inside
```

### 5.2 Main evolutionary loop
```
nbIterNonProd = 1
for tspid = popSizeMu*4 + 1 .. maxIter (=617..2000):
    if nbIterNonProd > maxIterNonProd (300) OR elapsed > timeLimit: break
    p1 = selectparents(feasiblePop, infeasiblePop, ...)
    p2 = selectparents(feasiblePop, infeasiblePop, ...)
    offspring_tsp = Crossover(p1, p2, vrp)
    rng(SEED + tspid)
    chromR = (rand<=0.5 ? split_Dmax : split_Tmax)(...)
    sol = chromR_detail_all(...)
    sol = ELS_mian(...)
    PopManagement(sol, ...)
    [.., isNewBest] = uti_updateBestSol(...)
    Last100Sol = uti_addSol2Last100(...)
    if sol.IsFeasible==0 AND rand-0.5>1e-6: Repair_sol(...)  # may set isNewBest
    if isNewBest: nbIterNonProd = 1  else: nbIterNonProd++   # non-productive counter
    if (tspid mod nbLast == 0) AND (tspid >= 100): ADAPT PENALTIES  (§5.7)
```

### 5.3 Parent selection (`selectparents`) — binary tournament over the UNION
```
a = size(feasiblePop) (0 if empty); b = size(infeasiblePop)
draw p1 = randi(a+b): if p1>a → infeasiblePop[p1-a] else feasiblePop[p1]
draw p2 likewise
return the one with LOWER Fitness (ties → p2)
```
Two independent draws per parent slot are done in the driver (it calls `selectparents` twice for `p1` and twice-equivalent for `p2`); each call itself already does a 2-candidate tournament. Net: each parent is the better of two uniformly-random members from the combined feasible+infeasible pool.

### 5.4 Termination
- Init phase stops at `popSizeMu*4` iterations (or maxIter/timeLimit).
- Main loop stops when `nbIterNonProd > 300` (no new overall-best for 300 iters), or `tspid` reaches `maxIter=2000`, or `timeLimit`.

### 5.5 Survivor selection
Handled entirely inside `PopManagement` (§4.2): trim to `popSizeMu` whenever size exceeds `popSizeMu + popSizeLambda`, clone-first then worst-biased-fitness, protecting the single best-cost member.

### 5.6 Repair probability
`rand - 0.5 > 1e-6` ⇒ probability ≈ 0.5 that an infeasible offspring is sent to `Repair_sol`. (Not in the 8 files; noted for completeness.)

### 5.7 Penalty adaptation (`targetFeasible` logic) — triggers every `nbLast=20` iters once `tspid>=100`
For each constraint X ∈ {T, C, D} independently:
```
fractionFeasible_X = (# of Last100 solutions with penalty_X == 0) / size(Last100)
origin_wX = current weight
if fractionFeasible_X <= targetFeasible - 0.05 (=0.15):     # too few feasible → raise
    wX = min(100000, wX * penaltyScaleFactor(1.2))
elif fractionFeasible_X >= targetFeasible + 0.05 (=0.25):   # too many feasible → lower
    wX = max(0.1, wX * penaltyDecreaseFactor(0.85))
# rescale every infeasible member's stored weighted penalty to the new weight:
for each s in infeasiblePop:
    s.penalty_X = (s.penalty_X / origin_wX) * wX
```
Then:
```
for each s in infeasiblePop:
    s.cost_Total = s.penalty_D + s.penalty_C + s.penalty_T + s.distance_Total
if size(infeasiblePop) > 1:
    infeasiblePop = infeasiblePop_updateBiasedFitnesses(infeasiblePop, vrp, par_hgs)   # §4.5
    # then re-sort by {cost_Total ASC, ID DESC} and rebuild vrp.ALL_brokenDIS from brokenDist
```
`wM` never adapts (stays 0). Weights are clamped to `[0.1, 100000]`.

### 5.8 Final result extraction
`bestSolOverall`: if infeasible → `Result = 99999`; else recompute exact distance from `chromR_move` over `vrp.distance_table` (with `+1` to go from `_move` ids back to `distance_table` 1-based node ids where depot=1), and `distance_table` is rounded to 2 decimals (`floor(100*d)/100`) at the very end.

---

## 6. LOCAL-SEARCH & AFS INSERTION (semantics the population files depend on)

The 8 target files contain **no local-search operators** — they only consume the `predecessors`/`successor`/`node_location`/penalty fields that ELS produces. For completeness of the module's contract:

- `ELS_mian` runs a granular neighborhood local search using operators `m1..m9`, `Depot_m1/2/3/8/9`, `NewRoute_m1/2/3`, ordered exactly as listed in `ELS_mian` lines 148-205, first-improvement (each returns `isSuccess`; on success it `continue`s the correlated-vertex loop). Candidate `V` for each `U` comes from `vrp.correlatedVertices(U,:)` (granular neighbor list, size `nbGranular=20`), optionally shuffled with probability `1/nbGranular` per `U`.
- **Conditional AFS Insertion** is embedded in the move operators and in `AFSdelay_new`/`deleteAFS_*`/`get_pd_pt`: after a move, an AFS is retained only if removing it (the "patch" pass in `ELS_mian` lines 60-144) would push a route segment distance over `V_Dmax`; i.e. an AFS is deleted (`deleteAFS_node`) when the merged segment `b = D_pre + D_su - d(pre,afs) - d(afs,su) + d(pre,su) <= V_Dmax`, otherwise kept. `node_location==100` marks the surviving AFS node.
- **What the population files require from this:** each stored solution must have consistent `predecessors`, `successor`, `routeID` (via `phrase_chromR`/`get_chromR`), `node_location` (row vector with `-1/1/100` markers), and correct weighted `penalty_*`, `distance_Total`, `cost_Total`, `IsFeasible`. Given those, §2–§4 are fully determined and independent of HOW the local search produced them.

Per your instruction, in the C++ port each candidate move is FULLY RE-EVALUATED (recompute route distances, PT/PC/PD, AFS schedule) rather than using the MATLAB O(1) deltas — the population-management semantics above are unaffected by that choice.

---

## 7. IMPLEMENTATION CHECKLIST / GOTCHAS

1. **`Fitness`: lower is better** everywhere (parent tournament, survivor removal picks the highest). Do not invert.
2. **`avgBrokenDist` is negated**; sort ascending ⇒ most-diverse first.
3. **`fitRank`/`divRank` normalize by `(N-1)`** ⇒ undefined for `N==1` (handled by the empty/first-member init paths that hardcode `0`). Guard division by zero.
4. **Clone threshold `1e-8`**, best-cost improvement margin `1e-9` — keep these exact.
5. **Penalty-adaptation band `targetFeasible ± 0.05`** and factors `1.2 / 0.85`, clamps `[0.1, 100000]`.
6. **Tie-breaks by ID DESC (newer first)** in `infeasiblePop_updateBiasedFitnesses` and the post-adaptation re-sort. Cost-sort tie-break also newer-first there.
7. **Depot = pred/succ value 0** in the linked-list representation; **depot = node 1** in `chromR`. AFS ids collapse to `nb_customer+1` inside the distance function only.
8. **Position 1 (best cost) is protected from survivor removal** (loop starts at index 2).
9. **`nbLast` is 20** despite "Last100" naming; window is FIFO.
10. **Broken-pair divisor**: use `vrp.last_customer` for the incremental/canonical path (§2.1), `nbClients` in the `update_*` path — but standardize on §2.1 (`vrp.last_customer`, collapse AFS) for the whole port; the only observable difference is when compared solutions differ in AFS multiset, and the collapse-version is the intended semantics.
11. You may replace ALL incremental matrix row/col shifting in `add2Pop`/`PopManagement` with a from-scratch N×N recompute; results are identical and much simpler in C++.

---

Reference files (absolute paths):
- Target module: `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Population management/{add2Pop,PopManagement,update_feasiblePop,update_infeasiblePop,infeasiblePop_updateBiasedFitnesses,uti_updateBestSol,uti_addSol2Last100,bestsolfind}.m`
- Solution schema & objective/penalties: `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Load/chromR_detail_all.m`
- Linked-list representation: `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/{phrase_chromR,get_chromR}.m`
- Control flow / constants / penalty adaptation: `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Main_METS.m`
- Parent selection: `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Crossover selectparents/selectparents.m`
- Local-search driver (field contract): `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/ELS_mian.m`