I now have a complete and grounded understanding of all five requested files plus the supporting context needed to make the spec self-contained. Here is the exhaustive specification.

---

# GrVRP-PCAFS — METS Efficient Local Search: C++ Implementation Specification

This spec covers the five core files (`ELS_mian`, `get_chromR`, `phrase_chromR`, `get_pd_pt`, `get_pc_now`) plus the supporting model semantics (`AFSdelay_new`/`AFSdelay_recursion`, `deleteAFS_node`, `deleteAFS_delta`, and the operator convention exemplified by `m1`) that these files depend on. It is written to be implementable **without** the MATLAB.

Absolute source paths:
- `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/ELS_mian.m`
- `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/get_chromR.m`
- `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/phrase_chromR.m`
- `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/get_pd_pt.m`
- `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/get_pc_now.m`
- (supporting) `AFSdelay_new.m`, `AFSdelay_recursion.m`, `deleteAFS_node.m`, `deleteAFS_delta.m`, `m1.m`

---

## 0. INDEXING CONVENTION (critical — read first)

MATLAB is 1-indexed. The problem has three node classes and TWO parallel numbering systems: the **abstract node ID** (used inside routes/predecessor/successor arrays) and the **distance-table index** (always `abstractID + 1`).

Let `nbClients = vrp.nb_customer` = number of customers.

**Abstract node IDs** (values stored in routes, predecessors, successor, chromR):
- Customer `c` has abstract ID `c`, for `c = 1 … nbClients`.
- Depot is abstract ID `0` inside `predecessors`/`successor` (0 means "no node / connects to depot"). In a chromR route array the depot is NOT stored as an element; routes are pure customer/AFS sequences bracketed implicitly by the depot.
- AFS (refueling station) nodes have abstract IDs `> nbClients`. The k-th route's AFS node (when present) has abstract ID `nbClients + k`. There is a **single physical AFS location** shared by all routes (see distance-table note below); each route that uses it gets its own AFS *node record* numbered `nbClients + routeIndex`.

**Distance-table index** = `abstractID + 1`:
- `vrp.distance_table` is `(nbClients + 2) × (nbClients + 2)` (or larger if multiple AFS locations, but this model uses one AFS).
- Row/col `1` = depot.
- Row/col `c+1` = customer `c`.
- Row/col `nbClients + 2` = the (single) AFS location. In code, `dAll(1, nbClients+2)` is depot↔AFS distance; `dAll(vrp.nb_customer+2, nodeY+1)` is AFS→nodeY. So **every AFS node, regardless of its abstract ID `nbClients+k`, maps to distance-table index `nbClients+2`.** When you need the distance table index of an AFS node, use `nbClients+2`, NOT `abstractID+1`.
- IMPORTANT MATLAB→C++: everywhere the MATLAB reads `dAll(A+1, B+1)` (nodes A,B abstract), translate to `distance_table[A][B]` if you 0-index your table with depot=0, customer c = c, AFS = nbClients+1. Keep a single consistent scheme and convert once.

In C++ prefer 0-based arrays but keep an explicit `toDistIdx(abstractID)` helper:
```
toDistIdx(id): if id == 0 return depot_index;
               else if id <= nbClients return customer_index(id);
               else return afs_index;   // the single AFS location
```

---

## 1. DATA STRUCTURES

### 1.1 `vrp` (problem instance, read-only)
Fields referenced by the five files + m1:
- `vrp.nb_customer` (int) — number of customers = `nbClients`.
- `vrp.distance_table` (matrix, `(nbClients+2)²`) — symmetric distance matrix, distance-table indexing (Section 0).
- `vrp.T_Customer` (double) — fixed service/handling time added per customer visit (`everTime`).
- `vrp.V_speed` (double) — vehicle speed; travel time = distance / speed.
- `vrp.V_Dmax` (double) — max distance per **leg** (a leg = depot→AFS or AFS→depot half of a route; if no AFS, the whole route is one leg). Distance penalty PD is per-leg over `Dmax`.
- `vrp.T_max_V` (double) — max route duration per vehicle (`Tmax`). Time penalty PT is over this.
- `vrp.T_Afs` (double) — refueling duration at the AFS (`T_Afs`), a fixed occupancy time.
- `vrp.C_Afs` (int) — AFS capacity: number of vehicles that may refuel simultaneously.
- `vrp.correlatedVertices` (matrix `nbClients × nbGranular`) — for each customer, its granular neighbor list (candidate `nodeV`s), abstract customer IDs. Row `i` lists the neighbors of customer `i`, ordered by proximity/relevance.

### 1.2 `par_hgs`
- `par_hgs.nbGranular` (int) — granular neighborhood size; also used as the modulus in the per-client neighbor-shuffle test (Section 5.2).

### 1.3 Route representation: `chromR`
`chromR` is a cell array (vector of int-vectors), one entry per route. Each `chromR{i}` is a column vector of abstract node IDs giving the visit order **within** route i (customers and, if present, the AFS node), NOT including the depot. Example route with an AFS: `chromR{i} = [3; 7; (nbClients+i); 2]` meaning depot → cust3 → cust7 → AFS → cust2 → depot. On the way INTO the ELS from an external representation, `chromR{i}` may contain a depot marker; `ELS_mian` strips it (see 2.1).

### 1.4 Adjacency representation: `Node_related` (the working state)
`Node_related` is an `N × 4` matrix where `N = nbClients + (number of AFS nodes currently in use)`. Rows `1..nbClients` are the customers (row `c` = customer `c`). Rows `nbClients+1 … N` are the AFS node records, row `nbClients+k` being the AFS of route k. Columns:
- Col 1 `predecessors[id]` — abstract ID of predecessor of node `id` within its route; `0` means "route starts here (predecessor is depot)". For an AFS node record that is *not currently placed*, this is `-1`.
- Col 2 `successor[id]` — abstract ID of successor; `0` means "route ends here (successor is depot)"; `-1` for an unplaced AFS node.
- Col 3 `routeID[id]` — 1-based route index the node belongs to; `-1` for an unplaced AFS node.
- Col 4 `node_location[id]` — position tag:
  - `-1` = node is on the **pre-AFS leg** (between depot and AFS).
  - `1` = node is on the **post-AFS leg** (between AFS and depot).
  - `100` = this node **is** the AFS node.
  - For a route with **no AFS**, all its customer nodes have `node_location = -1` (the whole route is treated as the pre-leg). (This is the convention `get_pd_pt` relies on: a route with `max(route) <= nbClients`, i.e. no AFS, is one leg.)

In C++, represent `Node_related` as 4 parallel `std::vector<int>` (`predecessors`, `successor`, `routeID`, `node_location`), each length `N`, 1-indexed logically (leave index 0 unused, or offset consistently). The vectors grow when an AFS is inserted and shrink when a route/AFS is deleted.

### 1.5 Per-route metrics: `Route_related`
`Route_related` is an `R × 7` matrix, R = number of routes. Row k describes route k. Columns:
- Col 1 `pd_v(k,1)` — distance penalty of the **pre-leg** of route k: `max(distance_pre(k) − V_Dmax, 0)`.
- Col 2 `pd_v(k,2)` — distance penalty of the **post-leg**: `max(distance_su(k) − V_Dmax, 0)`.
- Col 3 `pt_v(k)` — time penalty of route k: `max(time_afs(k) + time_su(k) − T_max_V, 0)`.
- Col 4 `time_afs(k)` — the **arrival time at the AFS** = duration of the pre-leg (travel + service). For a route with no AFS this holds the whole route duration `t1` and column 5 equals it (see below). ALSO doubles as an **AFS-usage flag holder** in patch code: value `2` marks "route has no AFS" in some patch branches (line 140-142: sets `Route_related(k,4)=2` when the AFS record has predecessor 0, meaning route degenerate). Treat col 4 primarily as `time_afs`.
- Col 5 `time_v(k)` — total route duration = `time_afs(k) + time_su(k)`.
- Col 6 `distance_pre(k)` — total distance of the pre-leg (depot → … → AFS), or the whole route if no AFS.
- Col 7 `distance_su(k)` — total distance of the post-leg (AFS → … → depot); `0` if no AFS.

Derived within operators: `time_su(k) = time_v(k) − time_afs(k)`.

`distance_Total = sum(sum(Route_related(:,[6,7])))` — sum of all leg distances = pure travel distance of the whole solution.

### 1.6 Penalty matrix: `Penalty_all` (2 × 4)
Row 1 = **weights**, Row 2 = **current penalty values**. Columns: 1=T, 2=C, 3=D, 4=M.
- `Penalty_all(1,1)=wt`, `Penalty_all(1,2)=wc`, `Penalty_all(1,3)=wd`, `Penalty_all(1,4)=wm` — penalty weights (adapted outside ELS).
- `Penalty_all(2,1)=pt` (weighted time penalty), `(2,2)=pc` (weighted AFS-capacity penalty), `(2,3)=pd` (weighted distance penalty), `(2,4)=pm` (weighted vehicle-count / missing-route penalty).

### 1.7 `afs_time_delay`
Column vector, length = number of AFS nodes = `R` (one per route that has an AFS). `afs_time_delay(k)` = the (possibly delayed) start-of-refueling time for route k's vehicle at the shared AFS, after resolving capacity conflicts. Entries equal `0` for routes with no AFS and are stripped inside `get_pc_now`. Carried in and out of ELS as `Route_related(:,4)`-seeded state. In `ELS_mian`, `afs_time_delay = Route_related(:,4)` initially (line 3).

### 1.8 `sol_table`, `sol_individual`
`sol_table` is a table (struct-of-arrays) of solutions indexed by `tspid`. ELS writes back into row `tspid` (see 2.1 output section). Relevant timing field: `sol_table.time` (vector of timestamps) used for the time-budget cutoff. `sol_individual` is passed through (returned unchanged by ELS body except as a slot). `sol_table.ID` is the vector of solution IDs (its length used in the cutoff formula).

### 1.9 Bookkeeping scalars (local to ELS)
- `whenLastTestedRI` (`nbClients × 1`, init 0) — per client, `nbMoves` value at which it was last tested.
- `whenLastModified` (`R × 1`, init 0) — per route, `nbMoves` value at which it was last modified. Grows/shrinks with routes.
- `nbMoves` (int, init 0) — count of accepted moves.
- `loopID` (int, init 0) — outer pass counter.
- `searchCompleted` (bool, init false).
- `isSuccess` (bool) — set by each operator; if true the inner V-loop `continue`s.
- `ispatch` (double) — snapshot of `sum(whenLastModified)`, used to detect that a modification happened since the last "patch/normalize" pass.

---

## 2. FUNCTION SPECS

### 2.1 `ELS_mian` — the local-search driver

**Signature (MATLAB):**
`[sol_individual, sol_table, Penalty_all, Route_related] = ELS_mian(sol_table, sol_individual, vrp, tspid, par_hgs, test, isrepair, chromR, nbClients, node_location, Penalty_all, Route_related, SEED)`

**Inputs:** current solution (`chromR`, `node_location`, `Route_related`, `Penalty_all`), `vrp`, `par_hgs`, `tspid` (index of this solution / offspring), `isrepair` (bool: are we in repair mode), `SEED` (RNG seed base), `sol_table` (for timing).

**Outputs:** updated `sol_individual`, `sol_table` (row `tspid` filled), `Penalty_all`, `Route_related`.

**Step-by-step:**

**(A) Setup (lines 3–22):**
1. `afs_time_delay = Route_related(:,4)`.
2. `correlatedVertices = vrp.correlatedVertices` (local copy; it will be shuffled in place per pass).
3. Build `chromR_move`:
   - If NOT repair: for each route, take `chromR{i}`, subtract 1 from every element, then delete all resulting zeros. (This removes the depot marker which was stored as `1`; i.e. the incoming chromR used 1-based-with-depot encoding, and `−1` maps depot(=1)→0 which is then dropped and shifts customer IDs. **Note:** after this transform, the customer previously stored as `c+1` becomes `c`. Implement to match: incoming route arrays are `[1, c1+1, c2+1, …, 1]`; output is `[c1, c2, …]`.)
   - If repair: `chromR_move = chromR` unchanged.
4. `[predecessors, successor, routeID] = phrase_chromR(vrp, tspid, nbClients, chromR_move)` (Section 2.3).
5. Build `Node_related`. Determine orientation of `node_location` (row vs column) by matching dimensions:
   - `a` = number of rows of `node_location`; `bb` = number of columns of `predecessors` (=1); `aa` = number of rows of `predecessors`.
   - If `a == bb` (node_location is a row vector): `Node_related = [predecessors, successor, routeID, node_location']`.
   - Elseif `a == aa`: `Node_related = [predecessors, successor, routeID, node_location]`.
   - Net effect: `Node_related` is `N×4` with columns as in 1.4, `node_location` as a column. In C++ just store node_location as length-N and skip the orientation dance.

**(B) Loop bookkeeping init (lines 24–30):**
`whenLastTestedRI = zeros(nbClients,1)`; `whenLastModified = zeros(R,1)`; `ispatch = sum(whenLastModified)` (=0); `nbMoves=0`; `loopID=0`; `searchCompleted=false`; `isSuccess=0`.

**(C) Outer loop `while ~searchCompleted` (lines 31–208):**

1. **Two-pass guard (lines 32–34):** if `loopID > 0` set `searchCompleted = true` at the *top* of the pass. This guarantees at least TWO passes (pass 0 and pass 1); pass 1 runs to completion but a successful move inside sets `searchCompleted=false` again, so more passes happen until a full pass yields no move. (Because operators set `searchCompleted=false` on success — see m1 line 309 — the loop terminates only when a whole pass makes zero accepted moves *and* `loopID>0`.)

2. **Per-pass neighbor shuffle (lines 37–42):** `rng(SEED + tspid)` (re-seed the RNG at the start of every pass so the shuffle sequence is deterministic and reproducible). Then for `i = 1..nbClients`: draw `r = randi(999999999)`; if `mod(r, par_hgs.nbGranular) == 0`, randomly permute row `i` of `correlatedVertices`. (So on average `1/nbGranular` of clients get their neighbor order reshuffled each pass.) **C++: reproduce MATLAB's `rng`/`randi`/`randperm` exactly only if bit-exact reproduction is required; otherwise use a seeded PRNG with the same modulus test and Fisher–Yates permutation.**

3. **Time-budget cutoff (lines 43–47):** if `tspid ~= 1`: if `toc − sol_table.time(end-1) > sol_table.time(end-1) * 20 / numel(sol_table.ID)` then `break` the outer loop. (`toc` = elapsed wall-clock since the matching `tic`. This aborts LS when it has consumed more than a per-solution share of time; the first solution `tspid==1` is exempt.)

4. **U-loop (lines 50–206):** for `ii = 1..nbClients`:
   - `nodeU = ii`; `correlatedU = correlatedVertices(nodeU,:)`; `lastTestRINodeU = whenLastTestedRI(nodeU)`; then set `whenLastTestedRI(nodeU) = nbMoves`.
   - **V-loop (lines 56–194):** for `jj = 1..length(correlatedU)`: `nodeV = correlatedU(jj)`; `afs_time_delay1 = afs_time_delay` (snapshot).

     - **(4a) NORMALIZE / PATCH pass (lines 60–144).** Runs only when `ispatch ~= sum(whenLastModified)` (i.e. a modification has occurred since the last normalize). On entry set `ispatch = sum(whenLastModified)`. This block re-canonicalizes `Route_related`/`Node_related` after moves that left them in a non-canonical state, and recomputes affected penalties. It does these things (with `tolerance = 1e-10`):
       - **Empty-route removal (lines 63–69):** find the first route `p` with `Route_related(p,6)≈0 && Route_related(p,7)≈0` (both leg distances zero → empty route). Delete row `p` from `Route_related`; decrement by 1 all `Node_related` routeIDs `> p`. (`break` after the first — only one per normalize entry.)
       - **Route reordering / AFS-only-route collapse loop (lines 70–120):** iterate `patch2 = 1 …` over routes. Maintain `sort_ind` (an `R×3` accumulator, later compacted). For each `patch2`:
         - **Case A — the route's max node_location is `100` and the last node (by index) of this route is exactly the AFS node `patch2+nbClients` (lines 75–103):** i.e. the AFS is the *last* stop before depot (post-leg empty). Test whether removing the AFS keeps the route feasible on distance:
           - `a = Route_related(patch2,6) + Route_related(patch2,7)` (total route distance).
           - `b = a − dist(pre_of_AFS → AFS) − dist(AFS → succ_of_AFS) + dist(pre_of_AFS → succ_of_AFS)` (distance if AFS spliced out). Uses distance-table indices with the AFS mapping.
           - If `b <= vrp.V_Dmax` (route fits in one leg without the AFS): **remove the AFS**:
             - Set `update_pc = 1`.
             - Rebuild `chromR` from Node_related via `get_chromR`; take route `patch2`, drop its max element (the AFS), giving `route_now` (customer-only).
             - `[patch_d,~,patch_t,~] = get_pd_pt(vrp, route_now, T_Customer, V_speed)` — recompute one-leg distance & time.
             - If `Route_related(patch2,3) ~= max(0, patch_t − T_max_V)` set `update_pt = 1` (time penalty changed).
             - Overwrite `Route_related(patch2,:) = [0, 0, max(0,patch_t−T_max_V), patch_t, patch_t, patch_d, 0]` (pd_pre=0,pd_su=0 assuming feasible distance, pt, time_afs=patch_t, time_v=patch_t, dist_pre=patch_d, dist_su=0).
             - Splice the AFS node out of the linked list (fix its pre/succ), then delete the AFS node record `patch2+nbClients` from Node_related, set all remaining nodes of route `patch2` to `node_location=−1`.
             - If `patch2` is not the last route, **rotate** this route to the end (move its Route_related row to the bottom and shift the rest up; renumber routeIDs accordingly; renumber AFS abstract IDs `>patch2+nbClients` down by 1), and decrement `patch2` so the shifted-in route is re-examined.
         - **Case B — the route's max node_location `~= 100` (no AFS) (lines 104–113):** if `patch2 <= (count of AFS nodes)`, this no-AFS route is positioned *before* AFS routes; rotate it to the end (so AFS routes are grouped first / no-AFS routes last), renumber, decrement `patch2`.
         - **Case C — has AFS but AFS is not the last node (line 115–117):** record in `sort_ind(patch2,:) = [lastNodeIndexOfRoute, patch2, lastNodeIndex − nbClients]` for a later re-sort of AFS records.
       - **AFS-record re-sort (lines 121–130):** compact `sort_ind`; if nonempty, sort the AFS rows of Node_related (`rows nbClients+1..end`) by routeID (col 3), then remap the abstract AFS IDs used inside predecessors/successor so they stay consistent with the new AFS row order. (This keeps AFS node `nbClients+k` aligned with route `k`.)
       - **Recompute changed penalties (lines 131–138):**
         - If `update_pt==1`: `Penalty_all(2,1) = Penalty_all(1,1) * sum(Route_related(:,3))` (re-sum time penalty × wt).
         - If `update_pc==1`: `[pc_now, afs_time_delay] = AFSdelay_new(Route_related(1:(N−nbClients),4), Route_related(1:(N−nbClients),5), T_max_V, T_Afs, C_Afs, afs_time_delay)` and `Penalty_all(2,2) = Penalty_all(1,2)*pc_now`.
       - **Degenerate-AFS flag (lines 139–143):** for each AFS record row `patch1`, if its predecessor is `0` and `Route_related(patch1−nbClients,4) ~= 2`, set `Route_related(patch1−nbClients,4) = 2` (marks that route's AFS as detached).

       **C++ NOTE:** because the port fully re-evaluates every candidate rather than doing incremental deltas, you may replace this entire patch block with a single canonicalization routine that, after any accepted move: (i) drops empty routes, (ii) removes AFS nodes whose post-leg is empty and whose route then fits `V_Dmax`, (iii) re-groups routes so AFS-routes precede no-AFS routes (only if you must preserve identical route numbering; otherwise route order is immaterial to cost), (iv) renumbers AFS node IDs to `nbClients+k`, and (v) recomputes PT and PC from scratch. The *semantics* above (feasibility test `b<=V_Dmax`, the `[0,0,pt,t,t,d,0]` row form) must be preserved; the exact rotation order only matters for reproducing identical tie-break sequences.

     - **(4b) OPERATOR CASCADE (lines 146–193).** Guard: run the cascade only if `loopID==0` OR `max(whenLastModified(routeID(nodeU)), whenLastModified(routeID(nodeV))) > lastTestRINodeU` (i.e. one of the two involved routes was modified after U was last tested — the classic HGS "don't retest unchanged pairs" filter). If the guard passes, try operators **in this exact order**, and on the FIRST that returns `isSuccess==true`, `continue` (skip the rest, move to next V):
       `m1, m2, m3, m4, m5, m6, m7, m8, m9`, then **only if `Node_related(nodeV,1) == 0`** (nodeV is a route-start, predecessor = depot): `Depot_m1, Depot_m2, Depot_m3, Depot_m8, Depot_m9`.
       Each operator has signature (all state passed & returned):
       `[isSuccess, nbMoves, searchCompleted, …, whenLastModified, Penalty_all, Route_related, Node_related, afs_time_delay] = mX(nodeU, nodeV, vrp, nbMoves, searchCompleted, tspid, par_hgs, whenLastModified, Penalty_all, Route_related, Node_related, afs_time_delay)`.
       See Section 4 for the operator semantics.

     - **(4c) NEW-ROUTE operators (lines 195–205).** After the V-loop, i.e. when `nodeV == correlatedU(end)` AND `loopID ~= 1`, try `NewRoute_m1, NewRoute_m2, NewRoute_m3` in order; `continue` on first success. (These create a brand-new route to place `nodeU` into. `loopID~=1` means: skip these on pass index 1 but allow on pass 0 and passes ≥2. Combined with the two-pass guard this is why "at least two loops" are needed — empty-route moves are only tried on the first pass.)

5. **End of U-loop → `loopID = loopID + 1`**, back to `while`.

**(D) Finalize & write-back (lines 211–238):**
- Rebuild `chromR = get_chromR(predecessors, successor, routeID, nbClients)`.
- `distance_Total = sum(sum(Route_related(:,[6,7])))`.
- `cost_Total = distance_Total + sum(Penalty_all(2,:))` (travel distance + all four current penalties).
- `IsFeasible = (sum(Penalty_all(2,:)) == 0) ? 1 : 0`.
- Store into `sol_table` at index `tspid`: `afs_time_delay`, `chromR_move` (=chromR), `distance_Total`, `cost_Total`, `predecessors`, `successor`, `routeID`, `node_location` (row), `distance_pre_su` (Route_related cols 6,7), `time_V` (col 5), `afs_time` (col 4), `penalty_D_v` (cols 1,2), `penalty_T_v` (col 3), scalar penalties `penalty_T/C/D/m` (=Penalty_all(2,1..4)), and `IsFeasible`.

---

### 2.2 `get_chromR` — reconstruct routes from linked list

**Signature:** `chromR = get_chromR(predecessors, successor, routeID, nb)`
(`nb` = nbClients, unused for math but part of interface.)

**Purpose:** given the predecessor/successor/routeID arrays (columns of Node_related), rebuild the ordered route arrays.

**Logic:**
```
chromR = cell(numberOfDistinctRouteIDs, 1)   // one slot per unique routeID
sort_candidate = 1 : max(max(predecessors), max(successor))   // all node IDs 1..maxNodeID
nn = numel(sort_candidate)
depot = all node IDs id in sort_candidate with predecessors[id] == 0   // route-start nodes
for each route-start node d = depot(i):
    route = [d]
    cur = d
    repeat up to nn-1 times:
        cur = successor[cur]
        if cur == 0: break     // reached depot end
        append cur to route
    nn = nn - length(route)    // (optimization only; safe to ignore in C++)
    chromR[ routeID[d] ] = route
remove empty cells from chromR
```
**Edge cases:** a route-start is any node whose predecessor is 0. Traversal stops when successor hits 0. `routeID[d]` places the route in the correct slot; unused/empty slots are stripped at the end. **C++:** iterate over all node IDs; those with `predecessors[id]==0 && routeID[id]!=-1` are starts; follow `successor` until 0; index the result by `routeID`.

---

### 2.3 `phrase_chromR` — build linked list from routes

**Signature:** `[predecessors, successor, routeID] = phrase_chromR(vrp, tspid, nbClients, chromR_move)`

**Purpose:** inverse of `get_chromR`. From the route arrays, fill the predecessor/successor/routeID arrays, sized to the largest node ID appearing.

**Logic:**
```
chromR = chromR_move
a = total number of elements across all routes   // sum of numel(chromR{i})
predecessors = successor = routeID = zeros(a, 1)   // NOTE: sized to 'a'
for j = 1 .. numberOfRoutes:
    if chromR{j} nonempty:
        for k = 2 .. numel(chromR{j}):
            predecessors[ chromR{j}(k) ]   = chromR{j}(k-1)
            successor[   chromR{j}(k-1) ]  = chromR{j}(k)
        routeID[ chromR{j} ] = j     // set routeID for every node in route j
```
**Notes / edge cases:**
- The arrays are indexed BY node ID (`chromR{j}(k)` is used as an index), so their logical length is `max node ID`, even though they are allocated to length `a` (= total element count). In practice for a valid solution `a >= maxNodeID`. The first node of each route keeps `predecessors=0` (route start) and the last keeps `successor=0` (route end) since only interior links are written.
- `vrp`, `tspid`, `nbClients` are unused in the body. **C++:** allocate the three arrays to `max node ID` (safer than `a`), zero-init, then wire interior links and set routeID.

---

### 2.4 `get_pd_pt` — leg distance & time for a route

**Signature:** `[d1, d2, t1, t2] = get_pd_pt(vrp, route_now, everTime, speed)`

**Inputs:** `route_now` = ordered abstract-ID sequence of ONE route (customers, and possibly the AFS node), NOT including depot; `everTime = vrp.T_Customer`; `speed = vrp.V_speed`.

**Outputs:** `d1` = pre-leg distance, `d2` = post-leg distance, `t1` = pre-leg time, `t2` = post-leg time. For a route with no AFS: `d1` = whole route distance, `t1` = whole route time, `d2=t2=0`.

**Logic:**
```
if route_now empty: return d1=d2=t1=t2=0

[a, n] = max(route_now)       // a = max ID, n = its position (1-based). The AFS has the largest ID.
if a > vrp.nb_customer:        // route CONTAINS an AFS at position n
    // Pre-leg = depot -> route_now(1) -> ... -> route_now(n) (the AFS)
    // Post-leg = route_now(n) -> ... -> route_now(end) -> depot
    D = 0; d1 = 0
    for i = 1 .. numel(route_now)-1:
        D += dist( route_now(i) -> route_now(i+1) )      // distance-table indexing
        if i == n-1:            // just closed the last hop of the pre-leg (…->AFS)
            d1 = D; D = 0       // freeze pre-leg accumulation, restart for post-leg
    d1 = d1 + dist( depot -> route_now(1) )              // add depot->first hop to pre-leg
    d2 = D + dist( route_now(end) -> depot )             // add last->depot hop to post-leg
    t1 = d1/speed + (n-1)*everTime      // pre-leg: (n-1) customers served before AFS (AFS itself no service time here)
    t2 = d2/speed + (numel(route_now) - n + 1)*everTime  // post-leg service count
else:                          // route has NO AFS
    D = 0
    for i = 1 .. numel(route_now)-1:
        D += dist( route_now(i) -> route_now(i+1) )
    d1 = D + dist(depot -> route_now(1)) + dist(route_now(end) -> depot)
    t1 = d1/speed + numel(route_now)*everTime
    d2 = 0; t2 = 0
```
**Distance indexing:** `vrp.distance_table(X+1, Y+1)` in MATLAB. For the AFS node (ID `nbClients+k`) the MATLAB code uses `route_now(i)+1` as the column index — BUT note `get_pd_pt` is called with `route_now` that already contains the AFS abstract ID `nbClients+k`, so `route_now(i)+1 = nbClients+k+1`. **This only maps to the correct AFS distance-table row `nbClients+2` when `k==1`.** In this single-AFS model the AFS record used in `get_pd_pt` re-evaluation contexts is always route-local and the caller (e.g. m1) constructs `routeU_now` with `numel(node_location)+1` as the AFS placeholder — you MUST map any node ID `> nbClients` to the single AFS distance-table index when porting. Concretely: `dist(X→Y) = distance_table[ toDistIdx(X) ][ toDistIdx(Y) ]` with `toDistIdx` from Section 0.

**Service-time semantics:** each customer visited adds `everTime`. In the AFS branch, the pre-leg counts `n-1` service events (the nodes strictly before the AFS position; the AFS at position `n` is not a customer so no service), and the post-leg counts `numel(route)-n+1` — note this INCLUDES position `n` in its count (the `+1`), i.e. the AFS position is billed to the post-leg service counter. **Reproduce these counts exactly** (`(n-1)` pre, `(len-n+1)` post) even though it looks asymmetric; it is how PT is computed.

---

### 2.5 `get_pc_now` — AFS capacity-overload penalty

**Signature:** `pc_now = get_pc_now(afs_time_delay, T_max, time_v, time_afs, T_Afs, C_Afs)`

**Purpose:** Given the (already conflict-resolved) refueling start times, compute the total "vehicle·time of over-capacity" at the AFS. This IS the raw PC (before multiplying by `wc`).

**Logic:**
```
afs_time_delay(afs_time_delay == 0) = []      // drop routes with no AFS
afs_time_end = afs_time_delay + T_Afs         // each vehicle occupies [start, start+T_Afs)
afs_time_end(afs_time_end == 0) = []          // (defensive; after removing zeros above)
c = sort( [afs_time_delay , afs_time_end] )   // all event times, ascending (2*m values, m = #AFS vehicles)
time_during = diff(c)                          // lengths of the (2m-1) elementary intervals

// For each interval, count how many vehicles are refueling:
d = afs_time_delay - c(1:end-1)                // (m × (2m-1)) broadcast: start_i minus interval-left-edge
overlap_matrix = (eps < d + T_Afs) & (d <= 0)  // vehicle i covers this interval iff start_i <= left AND start_i + T_Afs > left
nb_fueling = sum(overlap_matrix, 1)            // per interval: number of simultaneous vehicles
nb_fueling = max(nb_fueling - C_Afs, 0)        // amount over capacity
pc_now = max( nb_fueling * time_during' , 0)   // sum over intervals of (over-capacity * interval length)
```
**Meaning:** `pc_now` = ∫ max(concurrent_vehicles − C_Afs, 0) dt over the timeline — total vehicle-time of AFS congestion beyond capacity. `T_max`, `time_v`, `time_afs` are passed but not used in the body (kept for interface parity). `eps` is machine epsilon (use a tiny positive tolerance, e.g. `1e-12`).

**Vectorization detail for C++:** build the sorted event list of `2m` times; for each of the `2m-1` gaps `[c[j], c[j+1])`, count vehicles `i` with `afs_time_delay[i] <= c[j]` and `afs_time_delay[i] + T_Afs > c[j]` (strictly greater, with the `eps` tolerance); over-capacity `= max(count − C_Afs, 0)`; accumulate `over-capacity * (c[j+1]-c[j])`. Result clamped at ≥0.

---

## 3. OBJECTIVE & PENALTIES (exact formulas)

Let, per route k: `distance_pre(k)`, `distance_su(k)` (leg distances); `time_afs(k)` (pre-leg time), `time_su(k)` (post-leg time), `time_v(k)=time_afs(k)+time_su(k)`.

**Distance penalty PD (per leg over V_Dmax):**
```
pd_v(k,1) = max(distance_pre(k) - V_Dmax, 0)     // pre-leg
pd_v(k,2) = max(distance_su(k) - V_Dmax, 0)      // post-leg
PD (weighted) = pd = wd * sum_over_k( pd_v(k,1) + pd_v(k,2) )
```
So each *leg* (a leg is bounded by depot and AFS, or the whole route if no AFS) must fit within `V_Dmax`; overflow on either leg is penalized. This models the fuel/range constraint: a vehicle can travel at most `V_Dmax` between refuels.

**Time penalty PT (per route over T_max_V):**
```
pt_v(k) = max( time_afs(k) + time_su(k) - T_max_V , 0 )   // = max(time_v(k) - T_max_V, 0)
PT (weighted) = pt = wt * sum_over_k( pt_v(k) )
```
Route duration = travel time (`distance/speed`) + service time (`T_Customer` per customer). AFS refuel time `T_Afs` is NOT added into route duration here — it is handled separately by the capacity model. Overflow beyond `T_max_V` is penalized.

**AFS-capacity penalty PC:**
```
raw_pc = get_pc_now(afs_time_delay, …)             // Section 2.5 (uses conflict-resolved start times)
PC (weighted) = pc = wc * raw_pc
```
The conflict-resolved `afs_time_delay` comes from `AFSdelay_new`/`AFSdelay_recursion` (Section 3.1). PC = weighted total vehicle-time that the AFS is over its simultaneous-service capacity `C_Afs`.

**Vehicle-count / route penalty PM:**
```
pm = Penalty_all(2,4)      // maintained incrementally by operators
```
`pm` is carried in `Penalty_all(2,4)`. It changes only when a route is created or deleted. On AFS-route deletion (`deleteAFS_delta`): `pm_now = max(pm - wm, 0)` (one fewer route ⇒ subtract weight `wm`, floored at 0). New-route operators add `wm`. So `pm = wm * (numberOfRoutesBeyondSomeBaseline)` — effectively `wm ×` number of routes/vehicles used (a fixed cost per vehicle). **C++:** treat PM as `wm × (#routes − allowed)` or simply track additively as the operators do: +`wm` per created route, −`wm` (floored 0) per deleted route.

**Total cost / biased fitness used for acceptance:**
```
cost_Total = distance_Total + pt + pc + pd + pm
           = sum_over_legs(distance) + Penalty_all(2,1) + (2,2) + (2,3) + (2,4)
```
This penalized cost is what operators minimize. A move is accepted iff its `delta < -1e-6` where
```
delta = (change in travel distance)
      + (pm_now - pm) + (pt_now - pt) + (pc_now - pc) + (pd_now - pd)
      [ + AFS-open/close distance corrections, e.g.  -2*dAll(1, nbClients+2)  when an AFS-only route is deleted ]
```
`IsFeasible = (pt+pc+pd+pm == 0)`. The four weights `wt, wc, wd, wm` (row 1 of `Penalty_all`) are constant *within* one ELS call; they are adapted by the outer algorithm between calls via the targetFeasible mechanism (Section 5.3).

### 3.1 AFS delay model — `AFSdelay_new` / `AFSdelay_recursion` (semantics for PC)

`AFSdelay_new(time_afs, time_v, T_max, T_Afs, C_Afs, afs_time_delay)`:
```
if no AFS entries: pc_now = 0; return
time_V_shifting = max(T_max - time_v, 0)        // per route, slack: how much its AFS visit may be delayed without violating T_max
afs_time_delay = time_afs                        // initial refuel start = pre-leg arrival time at AFS
conflict_table = zeros(m,m)
[pc_now, afs_time_delay] = AFSdelay_recursion(time_afs, time_v, T_max, T_Afs, C_Afs, time_V_shifting, afs_time_delay, conflict_table)
```

`AFSdelay_recursion`: a scheduling/backtracking procedure that tries to *shift* AFS visit start times (each route's shift bounded by its slack `time_V_shifting`) so that no more than `C_Afs` vehicles overlap at the AFS, then returns the residual `pc_now` from `get_pc_now`. Algorithm:
```
pc_now = get_pc_now(afs_time_delay, …)          // baseline over-capacity
delay_over = true
while delay_over:
    delay_over = false
    for i = 1 .. m-1:
        for ii = i+1 .. m:
            // order the pair so t1 = earlier start, t2 = later start; q1,q2 = their slacks; t3 records which was which
            if afs_time_delay(i) >= afs_time_delay(ii): t1=delay(ii), t2=delay(i), q1=slack(ii), q2=slack(i), t3=1
            else:                                       t1=delay(i),  t2=delay(ii), q1=slack(i),  q2=slack(ii), t3=2
            if t1 + T_Afs <= t2: continue               // no overlap -> ok
            if (t1+T_Afs-t2 > q2) && (t2+T_Afs-t1 > q1): // neither can be delayed enough -> unresolvable pairwise
                conflict_table(i,ii)++; conflict_table(ii,i)++; continue
            if t1+T_Afs-t2 > q2:                         // only route with slack q1/later can move: delay t1 to t2+T_Afs
                delaytime = t2 + T_Afs - t1; t1 += delaytime
                apply delay to the appropriate route (per t3); reset i=ii=1 (restart scan); continue
            if t2+T_Afs-t1 > q1:                         // symmetric: delay t2 to t1+T_Afs
                delaytime = t1 + T_Afs - t2; t2 += delaytime
                apply; restart; continue
            // both delays individually feasible -> branch (backtrack): try delaying t1, recurse; try delaying t2, recurse; keep whichever yields smaller pc
            [afs_time_delay, time_V_shifting, pc_now] = backtrack(...)   // recursive, keeps min-pc branch
            delay_over = true; break
        if delay_over: break
pc_now = get_pc_now(afs_time_delay, …)          // final residual over-capacity
```
`backtrack` (nested): copies state; tries "delay vehicle-1 branch" (delay to `t2+T_Afs`), recurses `AFSdelay_recursion`; if resulting `pc_now` improves, keep it; then restores and tries "delay vehicle-2 branch" (delay to `t1+T_Afs`), recurses; keep the better of the two. Returns best `afs_time_delay`, `time_V_shifting`, `pc_now`.

**C++ semantics summary:** Given each AFS-using route's earliest AFS arrival time `time_afs(k)` and its slack `max(T_max - time_v(k), 0)`, find a set of *non-negative* start-time shifts (each ≤ its slack) minimizing the total over-capacity integral (`get_pc_now`). Greedily resolve pairwise overlaps by pushing one vehicle's start to just after the other finishes (when only one side has enough slack), and when both sides could resolve it, branch and keep the lower-penalty schedule. The returned `afs_time_delay` (final start times) and `pc_now` (residual) are what PC uses. **Note the known MATLAB bug at m1 line 232** (`AFSdelay_new(...)` called with only 5 args, missing `afs_time_delay`); the correct call passes `afs_time_delay` as the 6th argument — use the 6-arg form everywhere.

---

## 4. LOCAL-SEARCH OPERATORS — neighborhood, moves, feasibility

The cascade (m1..m9, then Depot_m1/2/3/8/9 when `nodeV` is a route-start, then NewRoute_m1/2/3 at end of V-loop) is **first-improvement**: each operator computes the single move it defines for the ordered pair `(nodeU, nodeV)`, accepts iff `delta < -1e-6`, and on acceptance mutates the state and returns `isSuccess=true`, which aborts the cascade for this pair.

Define, for the current pair, using the linked list:
- `preU = predecessors(nodeU)`, `nodeX = successor(nodeU)` (U's neighbors).
- `nodeY = successor(nodeV)` (V's successor).
- `routeU = routeID(nodeU)`, `routeV = routeID(nodeV)`.

**m1 — "Insert U after V" (relocate single node U to between V and Y).** Documented in full from `m1.m`; the other operators follow the same template with different move geometry. Its exact behavior:

- **Trivial reject:** if `nodeU == nodeY` (U already sits right after V) → fail.
- **Move geometry:** remove U from between `preU` and `nodeX` (splice `preU→nodeX`); insert U between `nodeV` and `nodeY` (`nodeV→U→nodeY`).
- **Cost of the splice:**
  - `costOne = −dist(preU→U) − dist(U→X) + dist(preU→X)` (distance saved by removing U from its route).
  - `costTwo = −dist(V→Y) + dist(U→Y) + dist(V→U)` (distance added by inserting U after V).
- **Same-route case (`routeU == routeV`):** apply `costOne`/`costTwo` to the affected leg of the single route (pre-leg if `node_location(U)==−1`, post-leg if `==1`; likewise for V's location choosing which leg V+Y edit lands on). Recompute that route's `pd_v`, `pt_v`, `time_afs`/`time_su`, then PC via `AFSdelay_new`. Time bookkeeping when moving across the service-count: removing U subtracts `everTime` from its leg time (`costOne/speed − everTime`), inserting U adds `everTime` (`costTwo/speed + everTime`).
- **Cross-route case (`routeU ≠ routeV`):** two subcases keyed on whether **routeV has an AFS** (`max(node_location of routeV)==100`):
  - **(a) routeV has an AFS:** plain relocation — apply costs to routeU's affected leg and routeV's affected leg, recompute both routes' `pd_v/pt_v/time`, recompute PC.
  - **(b) routeV has NO AFS — Conditional AFS Insertion (CRI rule):** inserting U into routeV may push routeV over `V_Dmax`, so the move inserts U **and a fresh AFS**: the new route becomes `[…V, U, AFS, Y…]`. It reconstructs routeV's full node order (`routeU_now`), appends `nodeU`, the new AFS placeholder (`numel(node_location)+1`), and the tail after V, then calls `get_pd_pt` to compute both legs (`d1,d2,t1,t2`) with the AFS present. `costTwo` includes the AFS edges: `−dist(V→Y) + dist(V→U) + dist(AFS→Y) + dist(U→AFS)`. A new AFS node record is appended (`predecessors=U, successor=Y, routeID=routeV, node_location=100`), and all nodes up to the AFS get `node_location=−1`, all after get `+1`.
- **Route deletion (`isdelete`):** if routeU had exactly 2 nodes and one was the AFS (`sum(routeID==routeU)==2 && max(node_location of routeU)==100`), then after removing U the route holds only the AFS → delete the whole route: `deleteAFS_delta` removes routeU's metric rows and does `pm_now = max(pm − wm, 0)`; `deleteAFS_node` removes the AFS node record and renumbers IDs/routeIDs. The delta then subtracts `2*dist(depot→AFS)` (the removed depot↔AFS round trip): `delta += −2*dAll(1, nbClients+2)`.
- **Acceptance:** `delta = costOne + costTwo + (pm_now−pm) + (pt_now−pt) + (pc_now−pc) + (pd_now−pd) [ −2*dAll(1,nbClients+2) if isdelete ]`. Accept iff `delta < −1e-6`. On accept: rewire predecessors/successor/routeID/node_location as above; rebuild `Route_related = [pd_v, pt_v, time_afs, time_v, distance_pre, distance_su]`; if a route was emptied to a single node, drop it and decrement higher routeIDs and drop its `whenLastModified`; set `Penalty_all(2,:) = [pt_now, pc_now, pd_now, pm_now]`; `nbMoves++`; `searchCompleted=false`; set `whenLastModified(routeU)=whenLastModified(routeV)=nbMoves` (or handle deletion: remove routeU's entry, set routeV's).

**Operator family semantics (for m2..m9, Depot_*, NewRoute_*).** They share the same evaluate-delta / accept-if-`<−1e-6` / rewire template; only the *move* differs. The standard HGS/METS move set these names imply (implement each as: build the resulting route(s), FULLY re-evaluate distance + PT + PD + PC + PM, accept iff penalized-cost delta `< −1e-6`, then rewire):

- **m1:** relocate U to after V (single node). [fully specified above]
- **m2:** relocate the pair `(U, X)` (U and its successor) to after V — insert `U,X` after V.
- **m3:** relocate the pair `(U, X)` reversed (`X,U`) to after V.
- **m4:** swap U and V (exchange two single nodes between their positions).
- **m5:** swap `(U,X)` with `V` (2-node segment for 1-node).
- **m6:** swap `(U,X)` with `(V,Y)` (2-node for 2-node).
- **m7:** 2-opt within a route or between routes — reverse the segment between U and V (intra-route edge exchange).
- **m8:** 2-opt* variant — exchange the tails of routeU and routeV after U and V (cross-route 2-opt joining `…U → Y…` and `…V → X…`).
- **m9:** the other 2-opt* tail-exchange orientation (`…U → V…` with reversed segment), complementary to m8.
- **Depot_m1/2/3/8/9:** the m1/m2/m3/m8/m9 moves specialized to the case where `nodeV` is a route-start (`predecessors(nodeV)==0`), i.e. inserting/joining at the *beginning* of routeV (before the first customer, adjacent to the depot). Only tried when `Node_related(nodeV,1)==0`.
- **NewRoute_m1/2/3:** create a NEW empty route and place U (m1: U alone; m2: `(U,X)`; m3: `(X,U)`) into it. Tried once per U at the end of its V-loop, only when `loopID ≠ 1`. These add a route so `pm_now = pm + wm` and add the depot↔first-node round-trip distance; accepted iff still `delta < −1e-6` (i.e. splitting off is worth the extra vehicle cost).

For EVERY operator, the AFS handling is the **Conditional AFS Insertion (CRI)** rule shown in m1 case (b): when a customer is added to a route that currently has no AFS and the new route would need refueling, insert exactly one AFS node so the route is split into a feasible pre-leg + post-leg; and the **AFS removal** rule (normalize pass / `deleteAFS`): when a route degenerates to only-AFS, or an AFS's post-leg becomes empty and the route then fits `V_Dmax` (test `b <= V_Dmax`), remove the AFS. Feasibility of any candidate is judged purely by whether the resulting per-leg distances exceed `V_Dmax` (→PD), per-route time exceeds `T_max_V` (→PT), and AFS concurrency exceeds `C_Afs` (→PC); these penalties enter the accept delta with weights `wd,wt,wc`, and `wm` for route count.

**C++ implementation guidance (per the task):** for each candidate move, construct the resulting route(s) as explicit node sequences, run `get_pd_pt` per route to get `d1,d2,t1,t2`, compute `pd_v`, `pt_v`, assemble `time_afs`/`time_v` for all routes, run `AFSdelay_new`+`get_pc_now` for PC, compute PM from route count, then `delta = newTotalPenalizedCost − oldTotalPenalizedCost` and accept iff `< −1e-6`. Do NOT port the O(1) incremental `costOne/costTwo/distance_pre±` bookkeeping; full re-evaluation reproduces the same acceptance decisions (the incremental math is exact-equal to full re-eval up to floating rounding, and the `−1e-6` tolerance absorbs that).

---

## 5. CONTROL FLOW SUMMARY

### 5.1 ELS termination
- Outer `while ~searchCompleted`. At the top of every pass after the first (`loopID>0`), `searchCompleted` is set true; any accepted move sets it back to false. Therefore ELS stops after the first *complete* pass (with `loopID≥1`) in which **no** move is accepted. Guarantees ≥2 passes because NewRoute/empty-route moves are skipped on pass `loopID==1` and empty-route checks are deferred to pass 0’s two-loop rule.
- Hard cutoff: if `tspid≠1` and elapsed time exceeds `sol_table.time(end-1) * 20 / numel(sol_table.ID)`, `break`.

### 5.2 Per-pass randomization
- Re-seed `rng(SEED + tspid)` each pass (reproducible).
- For each client, with probability ~`1/nbGranular` (test `mod(randi(999999999), nbGranular)==0`), randomly permute its `correlatedVertices` row. This perturbs neighbor exploration order between passes.

### 5.3 Penalty adaptation (targetFeasible — occurs OUTSIDE ELS, context for the port)
Within a single `ELS_mian` call the weights `wt, wc, wd, wm` are constant. The surrounding memetic/HGS driver adapts them between calls based on the fraction of recent solutions that were feasible vs a `targetFeasible` ratio: if too few feasible, multiply the relevant weight(s) up; if too many, multiply down (standard HGS multiplicative adaptation). `IsFeasible` (all four current penalties == 0) is the per-solution signal fed back. The C++ port should expose `Penalty_all(1,:)` (the weights) as inputs to the LS and adapt them in the driver, not inside the LS.

### 5.4 Repair probability & population/selection (context)
- `isrepair` toggles the input decoding (2.1 step A3) and whether the LS runs in repair mode; the driver decides per offspring whether to repair infeasible children (a repair probability governs this), running ELS again with `isrepair=true` and typically increased penalty weights so the result is feasible.
- Initial-population construction, parent selection (binary tournament on biased fitness = penalized cost + diversity contribution), and survivor selection (biased-fitness-based with clone/diversity pruning) live in the outer HGS driver (not in these five files). ELS is invoked on each new individual (`tspid` identifies it) to educate it; the educated individual is written back into `sol_table(tspid)`.

---

## 6. CONSTANTS / THRESHOLDS / TIE-BREAKS (checklist for the port)

- Acceptance threshold: `delta < −0.000001` (i.e. `−1e-6`). Moves with `delta > −1e-6` are rejected (note strict: exactly `−1e-6` is rejected because condition is `delta > −1e-6 ⇒ fail`).
- Empty-route detection tolerance: `1e-10` on `Route_related(:,6)` and `(:,7)`.
- AFS overlap epsilon in `get_pc_now`: `eps` (machine epsilon; use ~`1e-12`), condition `eps < d + T_Afs & d <= 0`.
- Neighbor-shuffle test: `mod(randi(999999999), par_hgs.nbGranular) == 0`.
- Time cutoff factor: `20 / numel(sol_table.ID)` times `sol_table.time(end-1)`.
- AFS-only-route deletion distance correction: `−2 * dist(depot ↔ AFS)` in the accept delta.
- Operator order is load-bearing (first-improvement): `m1→m2→…→m9→[Depot_m1,2,3,8,9 if V is route-start]→(end of V-loop, if loopID≠1) NewRoute_m1→m2→m3`.
- Service time `everTime = T_Customer` added once per customer visit; AFS visit adds no `everTime` in `get_pd_pt` but `T_Afs` occupancy is modeled only through PC.
- Leg definition: pre-leg = depot→…→AFS, post-leg = AFS→…→depot; no-AFS route = single pre-leg with `d2=t2=0`.
- Distance-table index of ANY AFS node = `nbClients+2` (single physical AFS), not `abstractID+1`.

---

### Files read to produce this spec (all absolute):
- `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/ELS_mian.m`
- `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/get_chromR.m`
- `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/phrase_chromR.m`
- `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/get_pd_pt.m`
- `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/get_pc_now.m`
- Supporting (read for correct semantics): `AFSdelay_new.m`, `AFSdelay_recursion.m`, `deleteAFS_node.m`, `deleteAFS_delta.m`, `m1.m`

Two source-code caveats the C++ port must NOT replicate: (1) m1.m line 232 calls `AFSdelay_new` with 5 args (missing `afs_time_delay`) — use the 6-arg form; (2) `AFSdelay_recursion`'s nested `backtrack` references `conflict_table` without receiving it as a parameter (a MATLAB closure/scope quirk) — in C++ pass `conflict_table` explicitly or drop it (it only accumulates unresolved-pair counts and does not affect the returned `pc_now`/`afs_time_delay`).