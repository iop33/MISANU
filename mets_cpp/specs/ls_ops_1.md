I have all the information I need. Note the distance-table indexing: `dAll(node+1, ...)` means the node ID stored in the arrays is the "node index" where depot=0's row is at MATLAB index 1, so a node value `n` addresses table row `n+1`. I have complete clarity on the AFS node = `nbClients+1` (value), addressing table column `nbClients+2`. Now I'll write the exhaustive spec.

Below is the complete specification.

---

# GrVRP-PCAFS — METS Local Search Operators m1–m5: C++ Implementation Specification

This spec is self-contained. It covers the five "efficient local search" move operators (`m1.m`–`m5.m`), all helpers they call (`AFSdelay_new`, `AFSdelay_recursion`/`backtrack`, `get_pc_now`, `get_pd_pt`, `deleteAFS_delta`, `deleteAFS_node`), the objective/penalty model, the data structures they read and write, and the driver loop (`ELS_mian.m`) that invokes them. Per the porting strategy, the incremental O(1) delta bookkeeping inside each operator is described, but the recommended C++ approach is to **fully re-evaluate** each candidate route's distance/time/penalties after applying the move; the incremental formulas below define the exact target values so you can validate a full re-evaluation against them.

---

## 0. GLOBAL CONVENTIONS AND INDEXING

MATLAB is 1-indexed. The problem uses a specific node numbering that MUST be preserved.

### 0.1 Node numbering (node "value", as stored in arrays)
- Let `nbClients = vrp.nb_customer` = number of customers.
- **Depot** has node value `0` in the predecessor/successor arrays (a predecessor of `0` means "attached to depot / start of route"). But in the distance table the depot is row/column `1`.
- **Customer** `c` (for `c = 1..nbClients`) has node value `c`.
- **AFS (refueling-station) nodes**: each route that visits an AFS has exactly ONE AFS node. AFS nodes are stored AFTER all customers in the node arrays. The i-th route's AFS node (if it exists) occupies array index `nbClients + i`. Its node *value* used for distance lookups is the "AFS super-node" — all AFS instances share one physical AFS location whose distance-table row/column is `nbClients + 2` (see 0.2). The AFS node's `node_location` marker value is `100`.

### 0.2 Distance table indexing — CRITICAL
`dAll = vrp.distance_table` is a square matrix. A node with stored value `n` addresses the table at row/column `n+1`:
- Depot (value 0) → table index `1`.
- Customer `c` (value c) → table index `c+1`.
- The shared AFS location → table index `nbClients + 2` (equivalently, the AFS "value" is `nbClients + 1`). In the operator code this appears both as `dAll(1, nbClients+2)` (depot→AFS) and `dAll(vrp.nb_customer+2, ...)` / `dAll(..., nbClients+2)`.

So in every operator, a distance between two nodes with values `a` and `b` is `dAll(a+1, b+1)`. When an AFS is freshly inserted, the code uses the AFS value `nbClients+1` for distance lookups (table index `nbClients+2`), NOT the per-route array index.

In C++: keep a `distance[i][j]` matrix (0-based) where index 0 = depot, index `c` = customer c, index `nbClients+1` = AFS location. Then "MATLAB `dAll(a+1,b+1)`" = C++ `distance[a][b]`.

### 0.3 node_location marker semantics
`node_location(k)` classifies node `k` relative to its route's AFS:
- `-1` : node is BEFORE the AFS on its route, OR the route has no AFS at all (whole route is the "pre" segment). Contributes to `distance_pre` / `time_afs`.
- `1`  : node is AFTER the AFS on its route. Contributes to `distance_su` / `time_su`.
- `100`: the node IS the AFS.
- `0`  : sometimes used transiently; a node value of `0` in the pred/succ arrays denotes the depot (route end/start).

Rationale: total route distance is split into `distance_pre` (depot → … → AFS) and `distance_su` (AFS → … → depot). If no AFS, everything is `distance_pre` and `distance_su = 0`.

### 0.4 Floating tolerance
Every acceptance test uses: a move is accepted iff `delta <= -1e-6` (written `delta > -0.000001 → reject`). Preserve this exact threshold.

---

## 1. DATA STRUCTURES

### 1.1 `vrp` (problem instance, read-only)
Fields used by m1–m5 and helpers:
- `vrp.nb_customer` (int): number of customers `nbClients`.
- `vrp.T_Customer` (double): fixed service time added per customer visit (`everTime`).
- `vrp.V_speed` (double): vehicle speed; travel time = distance / speed.
- `vrp.V_Dmax` (double): max distance a vehicle may travel on ONE fuel segment (pre-AFS OR post-AFS). Distance penalty is per-segment.
- `vrp.T_max_V` (double): max route duration (`T_max`).
- `vrp.T_Afs` (double): refueling duration at an AFS (occupies the AFS for `T_Afs` time).
- `vrp.C_Afs` (int): AFS capacity = number of vehicles that may refuel simultaneously.
- `vrp.distance_table` (matrix): see 0.2. Symmetric distance/cost matrix, 1-indexed with depot at 1 and AFS at `nbClients+2`.
- `vrp.correlatedVertices` (nbClients × K int matrix): for each customer `u`, the ordered granular neighbor list (candidate `V` nodes). Used only by the driver.
- `vrp.V_nb` (int): number of available vehicles (fleet size); excess beyond this is penalized (PM).

### 1.2 `par_hgs` (parameters, read-only in operators)
- `par_hgs.nbGranular` (int): granular neighborhood shuffle modulus (driver only).
- Other fields (`popSizeMu`, `popSizeLambda`, `targetFeasible`, `nbLast`, `nClosest`, `maxIter`, `maxIterNonProd`, `eliteNum`, etc.) belong to the outer genetic algorithm; not used inside m1–m5. `tspid` is the solution/individual ID (used only for RNG seeding and table storage).

### 1.3 `Node_related` (per-node state) — `N × 4` matrix
`N` = number of nodes currently in the solution = `nbClients + (number of routes that have an AFS)`. Columns:
1. `predecessors(k)`: node value preceding node `k` on its route (`0` = depot / route start).
2. `successor(k)`: node value following node `k` (`0` = depot / route end).
3. `routeID(k)`: which route node `k` belongs to (`1..R`). For a "free/unlinked" freshly created AFS node it is transiently `-1` before assignment.
4. `node_location(k)`: marker (see 0.3): `-1`, `1`, or `100`.

Row order: rows `1..nbClients` are the customers (row `c` = customer `c`). Rows `nbClients+1 .. N` are the AFS nodes, one per AFS-bearing route, stored in route order (the AFS node for route `i` is at row `nbClients + i` when routes are compacted). A newly inserted AFS is appended at row `N+1` (value `numel(node_location)+1`).

### 1.4 `Route_related` (per-route state) — `R × 7` matrix
`R` = number of routes. Columns:
1. `pd_v(i,1)` = `penalty_D_pre(i)` = `max(distance_pre(i) - V_Dmax, 0)`.
2. `pd_v(i,2)` = `penalty_D_su(i)` = `max(distance_su(i) - V_Dmax, 0)`.
3. `pt_v(i)`   = `penalty_T(i)` = `max(time_v(i) - T_max, 0)` (per-route time overage).
4. `time_afs(i)` = time at which the vehicle ARRIVES at its AFS on route `i` (no waiting), i.e. travel+service time of the pre-AFS segment. If no AFS, this is the full route time (see note in 1.6). Also used by AFS conflict logic. NOTE: value `2` is used as a sentinel in one driver patch branch (a route with `predecessors==0` AFS entry) — see driver §6; operators do not set that.
5. `time_v(i)`  = total route duration = `time_afs(i) + time_su(i)`.
6. `distance_pre(i)` = distance of the pre-AFS segment (depot → … → AFS), or full route distance if no AFS.
7. `distance_su(i)`  = distance of the post-AFS segment (AFS → … → depot), or `0` if no AFS.

Derived inside operators:
- `time_su(i) = time_v(i) - time_afs(i)`.

### 1.5 `Penalty_all` — `2 × 4` matrix
Row 1 = weights, Row 2 = current weighted penalty totals.
- `Penalty_all(1,1) = wt` (PT weight), `(1,2) = wc` (PC weight), `(1,3) = wd` (PD weight), `(1,4) = wm` (PM weight).
- `Penalty_all(2,1) = pt` = `wt * sum_i pt_v(i)` (total weighted time penalty).
- `Penalty_all(2,2) = pc` = `wc * penalty_C` (total weighted AFS-capacity penalty).
- `Penalty_all(2,3) = pd` = `wd * sum_i (pd_v(i,1)+pd_v(i,2))` (total weighted distance penalty).
- `Penalty_all(2,4) = pm` = `wm * max(R - V_nb, 0)` (total weighted vehicle-count penalty).

### 1.6 Distance/time split for a route (definition, from `chromR_detail_all.m`)
Given a route as a node sequence `[0(depot), n1, n2, …, nk, 0(depot)]` where exactly zero or one of the `n` are the AFS:
- If NO AFS: `distance_pre` = full cyclic distance depot→n1→…→nk→depot; `time_afs = distance_pre/speed + (k)*T_Customer` (service time counted for each of the k customers); `distance_su = 0`, `time_su = 0`; `time_v = time_afs`.  *(In `get_pd_pt`, the "no AFS" branch returns `t1 = d1/speed + numel(route_now)*everTime`, counting service for every node in the segment list, which are all customers.)*
- If AFS present at position `p` (route = pre-segment `[depot, …, AFS]` then post-segment `[AFS, …, depot]`):
  - `distance_pre` = depot→…→AFS distance; `distance_su` = AFS→…→depot distance.
  - `time_afs` = `distance_pre/speed + (#customers strictly before AFS)*T_Customer`. (Arrival time at AFS; refuel time `T_Afs` NOT included in `time_afs`.)
  - `time_su` = `distance_su/speed + (#customers after AFS)*T_Customer`. (Refuel time is NOT added to route duration in these formulas; it only matters for AFS-capacity conflicts via `get_pc_now`.)
  - `time_v = time_afs + time_su`.

The exact `get_pd_pt` implementation for a given route list `route_now` (containing customers, possibly the freshly-inserted AFS value `N+1`, with depot implicit at both ends) is in §2.6.

---

## 2. HELPER FUNCTIONS (implement first)

### 2.1 `get_pc_now(afs_time_delay, T_max, time_v, time_afs, T_Afs, C_Afs) → pc_now`
Computes the AFS-capacity penalty PC = total (over-capacity vehicles × overlap duration) at the single shared AFS.

Inputs: `afs_time_delay` = vector of (possibly delayed) AFS arrival times, one per route (0 if that route has no AFS or arrival time is 0). `T_Afs`, `C_Afs` as above. (`T_max`, `time_v`, `time_afs` are passed but unused in this function.)

Algorithm (exact):
```
afs_start = afs_time_delay with all entries == 0 removed
afs_end   = afs_start + T_Afs        (elementwise)   ; then remove any == 0 (none normally)
c = sort( concat(afs_start, afs_end) )   // ascending, length 2m for m AFS-visiting routes
time_during(k) = c(k+1) - c(k)  for k=1..len(c)-1     // interval widths
for each interval k (k=1..len(c)-1):
    // count vehicles refueling during interval [c(k), c(k+1)]
    d_i = afs_start_i - c(k)                     // for each start time i
    overlap_i = ( eps < d_i + T_Afs )  AND  ( d_i <= 0 )
        // i.e. start_i <= c(k)  and  start_i + T_Afs > c(k)  (strictly, with eps guard)
    nb_fueling(k) = count of i where overlap_i is true
nb_fueling = max(nb_fueling - C_Afs, 0)              // elementwise over intervals
pc_now = max( sum_k nb_fueling(k) * time_during(k), 0 )
```
`eps` = MATLAB machine epsilon ≈ 2.22e-16. If `afs_time_delay` is all zeros/empty, `pc_now = 0`.

Meaning: it lays out all refuel windows `[start, start+T_Afs)` on a timeline, and for each elementary interval counts how many windows cover it; any count above `C_Afs` incurs penalty proportional to interval width. This is the raw (unweighted) PC; the caller multiplies by `wc`.

### 2.2 `AFSdelay_new(time_afs, time_v, T_max, T_Afs, C_Afs, [afs_time_delay]) → (pc_now, afs_time_delay)`
Wrapper that (a) if `numel(time_afs)==0` returns `pc_now=0`; else (b) sets `time_V_shifting = max(T_max - time_v, 0)` (slack each route has before violating T_max), initializes `afs_time_delay = time_afs` (start from the un-delayed arrival times), and calls `AFSdelay_recursion`.

`time_afs`/`time_v` here are the per-route AFS-arrival and total-time vectors, restricted to the AFS-bearing routes (in the operators the caller passes `time_afs(1:numel(node_location)-nbClients)` = first `R_afs` entries; effectively the routes that have AFS nodes). `time_V_shifting(i)` is how much route `i`'s AFS visit can be delayed without pushing `time_v` past `T_max`.

### 2.3 `AFSdelay_recursion(...)` + `backtrack(...)` → resolves AFS scheduling conflicts, returns `pc_now`, `afs_time_delay`
This tries to DELAY some vehicles' AFS arrival times (within their slack `time_V_shifting`) to reduce simultaneous-refuel overlap, i.e. to minimize PC.

Algorithm:
```
pc_now = get_pc_now(afs_time_delay, ...)     // initial penalty with current delays
delay_over = true
while delay_over:
    delay_over = false
    for i = 1 .. m-1:               // m = number of AFS routes
      for ii = i+1 .. m:
        // order the pair so t1 = earlier arrival, t2 = later arrival
        if afs_time_delay(i) >= afs_time_delay(ii):
            t1=afs_time_delay(ii); t2=afs_time_delay(i); q1=shift(ii); q2=shift(i); t3=1
        else:
            t1=afs_time_delay(i);  t2=afs_time_delay(ii); q1=shift(i);  q2=shift(ii); t3=2
        if t1 + T_Afs <= t2: continue          // no overlap → no conflict
        // Overlap exists. Two "infeasible-to-separate" checks:
        if (t1+T_Afs-t2 > q2) AND (t2+T_Afs-t1 > q1):
            conflict_table(i,ii)++ ; conflict_table(ii,i)++ ; continue   // cannot separate; leave overlap
        if t1+T_Afs-t2 > q2:                   // can only fix by delaying t1
            delaytime = t2 + T_Afs - t1
            t1 += delaytime
            push t1 back into afs_time_delay of the later/earlier route per t3;
            reduce that route's shift by delaytime; reset i=ii=1; continue
        if t2+T_Afs-t1 > q1:                   // can only fix by delaying t2
            delaytime = t1 + T_Afs - t2
            t2 += delaytime
            (symmetric update); reset; continue
        // Both delays feasible → branch (backtracking) to try both and keep the better:
        [afs_time_delay, shift, pc_now] = backtrack(...)
        delay_over = true; break
      if delay_over: break
pc_now = get_pc_now(afs_time_delay, ...)      // final penalty
```
`backtrack` recursively tries "delay t1 by (t2+T_Afs−t1)" vs "delay t2 by (t1+T_Afs−t2)", recomputes via `AFSdelay_recursion`, and keeps whichever yields lower `pc_now`.

**Porting note:** This is a heuristic conflict-resolver whose only externally-observed output is the scalar `pc_now` (and the delayed arrival vector, stored back into `Route_related(:,4)` as `afs_time_delay` in the driver, but NOT re-read by operators — operators always recompute from `time_afs`). There is a known bug in the MATLAB (`i=1;ii=1` inside the for-loop does not actually restart the loops in MATLAB; and `backtrack` references an out-of-scope `conflict_table` and can error). For a **faithful port you may either** (a) replicate the exact loop including its quirks, or (b) since PC is re-derivable, reimplement a correct min-overlap scheduler. The initial `chromR_detail_all.m` version (lines 136–221) is a cleaner, self-contained variant that does the same pairwise separation with explicit tie-break `if q1 > q2 delay the later, else delay earlier` and no recursion; that version is the more reliable reference for what PC "should" be. **Recommendation:** implement `get_pc_now` exactly, and implement the delay-resolution using the `chromR_detail_all.m` pairwise algorithm (deterministic, no recursion), because operators fully recompute PC each call anyway.

### 2.4 `get_pd_pt(vrp, route_now, everTime, speed) → (d1, d2, t1, t2)`
Computes pre/post-AFS distance and time for a route given as a **node-value list `route_now`** (customers plus, if present, the AFS placeholder value). Depot is implicit at both ends (added via `distance_table(1, first+1)` and `distance_table(last+1, 1)`).

```
if route_now is empty: return (0,0,0,0)
[a, n] = max(route_now)              // a = largest value, n = its position (1-based)
if a > nb_customer:                  // an AFS placeholder is present (its value = N+1 > nb_customer, largest)
    // AFS is at position n; split route there
    D = 0; d1 = 0
    for i = 1 .. len-1:
        D += distance_table( route_now(i)+1, route_now(i+1)+1 )
        if i == n-1:                 // reached the edge into the AFS
            d1 = D; D = 0
    d1 += distance_table(1, route_now(1)+1)          // add depot → first
    d2  = D + distance_table( route_now(end)+1, 1 )  // add last → depot
    t1 = d1/speed + (n-1)*everTime          // (n-1) customers before AFS   (AFS at pos n)
    t2 = d2/speed + (len - n + 1)*everTime  // customers from AFS pos to end
    // NOTE: this service-time count is as written in the MATLAB; it counts positions,
    //       treating the AFS placeholder position boundaries per the (n-1)/(len-n+1) split.
else:                                // no AFS
    D = sum of consecutive distances along route_now
    d1 = D + distance_table(1, route_now(1)+1) + distance_table(route_now(end)+1, 1)
    t1 = d1/speed + len*everTime
    d2 = 0; t2 = 0
```
This is used only in the "insert into a route that has no AFS, so we must add one" branches (CRI rule), where `route_now` is rebuilt as `[pre-part…, movedNode(s), AFS_placeholder(=N+1), post-part…]`. The returned `(d1,d2,t1,t2)` become the new route's `distance_pre, distance_su, time_afs, time_su`.

### 2.5 `deleteAFS_delta(delete_idx, node_location, pd_v, pt_v, time_v, time_afs, distance_pre, distance_su, pm_now, pm, wm)` → updates penalty/route vectors after a route becomes empty
Triggered when a move empties route `delete_idx` down to only its AFS node (route must be deleted). It:
- Removes row `delete_idx` from `pd_v, pt_v, time_v, time_afs, distance_pre, distance_su`.
- Decrements the vehicle-count penalty: `pm_now = max(pm_now - wm, 0)` (one fewer route used).
- (`node_deleteAfs = delete_idx + sum(node_location~=100)` is computed for a sanity print only.)

### 2.6 `deleteAFS_node(routeID, delete_idx, node_location, predecessors, successor)` → renumbers node/route arrays after deletion
Removes the AFS node of the deleted route and compacts IDs:
- `node_deleteAfs = delete_idx + sum(node_location ~= 100)` — the array index of the AFS node to delete. (`sum(node_location~=100)` = number of non-AFS nodes = nbClients; so this is `delete_idx + nbClients`, the AFS row for route `delete_idx`.)
- `routeID(routeID > delete_idx) -= 1` (shift down route IDs above the deleted one).
- `predecessors(predecessors > node_deleteAfs) -= 1`; `successor(successor > node_deleteAfs) -= 1` (shift node values referencing rows above the removed AFS row).
- Delete row `node_deleteAfs` from `node_location, predecessors, successor, routeID`.

**Porting note:** In C++ prefer explicit route/node objects and rebuild; the exact `−1` renumbering just keeps the packed-array invariant.

---

## 3. OBJECTIVE, PENALTIES, AND ACCEPTANCE

### 3.1 Weighted penalties (current totals, kept in `Penalty_all` row 2)
```
PT (pt) = wt * Σ_i max(time_v(i) − T_max, 0)
PC (pc) = wc * pc_now                         // pc_now from AFSdelay_new
PD (pd) = wd * Σ_i [ max(distance_pre(i) − V_Dmax, 0) + max(distance_su(i) − V_Dmax, 0) ]
PM (pm) = wm * max(R − V_nb, 0)               // R = number of routes
```
Note PD penalizes EACH fuel segment (pre and post) separately against `V_Dmax`. Note refuel time `T_Afs` is NOT part of route duration for PT; it only affects PC via overlap.

### 3.2 Total cost of a solution
```
distance_Total = Σ_i (distance_pre(i) + distance_su(i))
cost_Total     = distance_Total + PT + PC + PD + PM
IsFeasible     = (PT + PC + PD + PM == 0)
```

### 3.3 Move acceptance (`delta`) — the general formula used by every operator
For a candidate move, compute the NEW weighted penalty totals `pt_now, pc_now, pd_now, pm_now` (using the same weights) after applying the move to a scratch copy of the affected routes, and the raw distance change `costOne + costTwo` (sum of edge-length deltas from removing node(s) from their old position and inserting at the new one). Then:
```
delta = costOne + costTwo
      + (pm_now − pm) + (pt_now − pt) + (pc_now − pc) + (pd_now − pd)
      [ − 2*dAll(1, nbClients+2)   if this move deletes a route that had an AFS ]
```
The extra `−2*dAll(1, nbClients+2)` term (twice depot↔AFS distance) accounts for the depot→AFS and AFS→depot edges disappearing when an AFS-only route is removed (that raw distance was part of `distance_Total` but is not captured by `costOne/costTwo`, which only track edges around the moved nodes).

**Acceptance rule:** apply the move iff `delta <= −1e-6`. Otherwise return `isSuccess = false` and leave the solution unchanged. This is a **first-improvement** local search (the first neighbor operator that yields an accepted move wins; the driver then `continue`s).

**Full-reevaluation strategy (recommended for C++):** For each candidate: (1) build the post-move route(s) as explicit node lists; (2) recompute `distance_pre/su`, `time_afs/su`, `time_v` for the affected routes via the §1.6 / §2.4 definitions; (3) recompute `pd_now, pt_now` by summing over ALL routes, `pc_now` via the AFS scheduler over all AFS routes, `pm_now = wm*max(R'−V_nb,0)` where `R'` is the new route count; (4) compute `costOne+costTwo` OR equivalently `newDistanceTotal − oldDistanceTotal` (they are equal, minus the AFS-edge correction which is naturally captured if you use full distance totals); (5) `delta = (newDistanceTotal + pt_now+pc_now+pd_now+pm_now) − (oldDistanceTotal + pt+pc+pd+pm)`; accept iff `delta <= −1e-6`. Using full totals AUTOMATICALLY includes the `−2*dAll(1,AFS)` correction, so you do not special-case it.

---

## 4. THE CONDITIONAL AFS INSERTION (CRI) RULE — shared logic

Several operators move customer(s) into a target route `routeV`. Whether the moved node(s) may go into `routeV` depends on whether `routeV` already has an AFS:

- **If `routeV` HAS an AFS** (`max(node_location(routeID==routeV)) == 100`): the moved node(s) are simply spliced in after `nodeV`, inheriting `nodeV`'s `node_location` (`-1` if inserted in the pre-segment, `1` if in the post-segment). No new AFS is added. The route may now violate `V_Dmax` on the affected segment (penalized, not forbidden).
- **If `routeV` has NO AFS** (`max(node_location(routeID==routeV)) ~= 100`): **CRI inserts a NEW AFS** into `routeV`, placed immediately AFTER the moved node(s) (pattern `[…, nodeV, movedNode(s), AFS, nodeY, …]`). Concretely:
  - A new AFS node is appended at array index `node_addAfs = numel(node_location)+1`, with `node_location=100`, `routeID=routeV`.
  - The route is rebuilt as `routeU_now = [pre-part(from nodeV back to route start), movedNode(s), node_addAfs, post-part(from nodeV forward to route end)]`, and `get_pd_pt` recomputes its pre/su distances and times (§2.4). Everything up to and including the new AFS gets `node_location=-1`; from the AFS to the end gets `1`; the AFS itself `100`.

The choice is made ONLY on the presence/absence of an AFS in `routeV`; the CRI always places the AFS directly after the just-inserted customer(s). There is no search over AFS insertion positions here.

### 4.1 Route-deletion side effect
When a customer (or pair) is moved OUT of its source route `routeU` and that leaves `routeU` containing ONLY its AFS node, the source route is deleted:
- m1 condition: `sum(routeID==routeU)==2 && route has an AFS` → `isdelete=1`.
- m2/m3 condition: `sum(routeID==routeU)==3 && route has an AFS` → `isdelete=1` (because two nodes leave).
- On `isdelete`: call `deleteAFS_delta` (drop the route's penalty rows, `pm_now = max(pm_now−wm,0)`), recompute PC/PT/PD over the reduced set, subtract `2*dAll(1,nbClients+2)` in `delta`, then `deleteAFS_node` to renumber. Also `whenLastModified(routeU)` is removed and `whenLastModified(routeV)=nbMoves`.

Only m1, m2, m3 (single-direction insertions) can trigger `isdelete`. m4 and m5 (swaps) never delete a route (node counts are preserved).

---

## 5. OPERATORS m1–m5 — EXACT SEMANTICS

Common preamble for ALL operators (read from the state arrays):
```
nbClients, everTime(=T_Customer), speed, V_Dmax, T_max, T_Afs, C_Afs   ← vrp
wt,wc,wd,wm ← Penalty_all(1,:) ;  pt,pc,pd,pm ← Penalty_all(2,:)
pd_v ← Route_related(:,[1,2]); pt_v ← Route_related(:,3)
time_v ← Route_related(:,5); time_afs ← Route_related(:,4); time_su = time_v − time_afs
distance_pre ← Route_related(:,6); distance_su ← Route_related(:,7)
predecessors,successor,routeID,node_location ← Node_related columns 1..4
For nodeU: nodeU_loc=node_location(nodeU); routeU=routeID(nodeU); preU=pred(nodeU); nodeX=succ(nodeU)
For nodeV: nodeV_loc=node_location(nodeV); routeV=routeID(nodeV); nodeY=succ(nodeV)
yu=routeU; yv=routeV  (returned for driver's whenLastModified logic)
```
`nodeX_loc`, `nodeY_loc`, `suX = succ(nodeX)`, `preV = pred(nodeV)` are read where needed. A value of `0` for a successor/predecessor means depot (route boundary). The incremental distance/time updates below are the MATLAB's O(1) deltas; in C++ you may instead fully recompute the two affected routes and derive `delta` per §3.3.

Throughout: when a node's `node_location == −1`, its incident edge-length change flows into `distance_pre`/`time_afs`; when `== 1`, into `distance_su`/`time_su`. The `everTime` adjustments (`± k*everTime`) account for customers moving between the pre and post segments (which changes the service-time bookkeeping of `time_afs` vs `time_su`).

---

### 5.1 m1 — "Insert u after v" (relocate single customer u to just after v)

**Move:** remove `nodeU` from between `preU` and `nodeX`; insert it between `nodeV` and `nodeY`. Resulting local structure: `… preU nodeX …` (source) and `… nodeV nodeU nodeY …` (target).

**Early rejects:**
- If `nodeU == nodeY` (u already directly after v): `isSuccess=false; return`.

**Case A — same route (`routeU == routeV`):**
- `costOne = −dAll(preU+1,nodeU+1) − dAll(nodeU+1,nodeX+1) + dAll(preU+1,nodeX+1)` (close the gap left by u).
- `costTwo = −dAll(nodeV+1,nodeY+1) + dAll(nodeU+1,nodeY+1) + dAll(nodeV+1,nodeU+1)` (open a slot for u after v).
- Apply the segment-distance/time deltas: for the REMOVAL, if `nodeU_loc==−1` add `costOne` to `distance_pre(routeU)`, `time_afs += costOne/speed − everTime` (u leaves pre-segment, so one fewer serviced customer there); if `nodeU_loc==1` add to `distance_su`, `time_su += costOne/speed − everTime`. For the INSERTION, using `nodeV_loc`: `+costTwo` to the corresponding segment, and `time += costTwo/speed + everTime` (u joins that segment). Recompute `pd_v(routeU,:)`, `pt_v(routeU)`, `pc_now` (via `AFSdelay_new`), `pm_now=pm`.
- `delta = costOne+costTwo + (pm_now−pm)+(pt_now−pt)+(pc_now−pc)+(pd_now−pd)`. Reject if `> −1e-6`.
- Apply pointers: `pred(nodeU)=nodeV; succ(nodeU)=nodeY; succ(nodeV)=nodeU; routeID(nodeU)=routeV; node_location(nodeU)=nodeV_loc; if preU: succ(preU)=nodeX; if nodeX: pred(nodeX)=preU; if nodeY: pred(nodeY)=nodeU`.

**Case B — different routes (`routeU != routeV`):**
- Compute `isdelete`: if source route has exactly 2 nodes AND has an AFS (`sum(routeID==routeU)==2 && route has AFS`), deleting u empties it → `isdelete=1`.
- **B1 — routeV HAS an AFS** (`max(node_location(routeID==routeV))==100`):
  - Same `costOne`, `costTwo` as Case A. Removal deltas applied to `routeU` segment; insertion deltas applied to `routeV` segment (per `nodeV_loc`). Recompute `pd_v` for both `routeU` and `routeV`, `pt_v` both, `pc_now`, `pm_now`.
  - If `isdelete`: `deleteAFS_delta(routeU,…)`, recompute PC/PT/PD, `pm_now=max(pm−wm,0)`.
  - `delta = costOne+costTwo+Δpenalties [ −2*dAll(1,nbClients+2) if isdelete ]`. Reject if `>−1e-6`.
  - Apply pointers (relocate u into routeV as in Case A but with cross-route `routeID(nodeU)=routeV`), then if `isdelete` call `deleteAFS_node`.
- **B2 — routeV has NO AFS** (CRI insert): result is `[…, nodeV, nodeU, AFS_new, nodeY, …]`.
  - `costTwo = −dAll(nodeV+1,nodeY+1) + dAll(nodeV+1,nodeU+1) + dAll(nbClients+2, nodeY+1) + dAll(nodeU+1, nbClients+2)` (open slot for u AND for a new AFS between u and nodeY; distances to/from the shared AFS index `nbClients+2`).
  - Rebuild `routeU_now` = walk `routeV` from `nodeV` forward (successors) → `routeU_NOW`, walk backward (predecessors) → reversed prefix; then `routeU_now = [prefix, nodeU, (N+1), routeU_NOW]` where `N+1 = numel(node_location)+1` is the AFS placeholder. Call `get_pd_pt` → `(d1,d2,t1,t2)` = new `distance_pre/su`, `time_afs/su` of `routeV`.
  - Apply removal delta to `routeU`; set `routeV`'s `distance_pre=d1, distance_su=d2, time_afs=t1, time_su=t2`. Recompute `pd_v`, `pt_v` for both routes; `pc_now`; `pm_now`.
  - If `isdelete`: as in B1.
  - `delta = costOne+costTwo+Δpenalties [ −2*dAll(1,nbClients+2) if isdelete ]`. Reject if `>−1e-6`.
  - Append the new AFS node at `node_addAfs = numel(node_location)+1` (`node_location=100, pred=succ=−1, routeID=−1` initially). Splice pointers: `pred(nodeU)=nodeV; succ(nodeU)=node_addAfs; succ(nodeV)=nodeU; pred(node_addAfs)=nodeU; succ(node_addAfs)=nodeY; if nodeY: pred(nodeY)=node_addAfs; routeID(nodeU)=routeV; routeID(node_addAfs)=routeV`. Fix source pointers (`succ(preU)=nodeX` — but if `isdelete` and `preU` is the source route's AFS, skip that repoint), `pred(nodeX)=preU`. Set `node_location` along `routeU_now`: everything up to & incl. `node_addAfs` = `−1`, from `node_addAfs` to end = `1`, the AFS itself = `100`. If `isdelete`: `deleteAFS_node`.

**Postamble (all cases):** rebuild `Route_related=[pd_v, pt_v, time_afs, time_v, distance_pre, distance_su]`; if the source route ended up with a single node (`sum(Node_related col3==routeU)==1`) drop that empty route row and decrement higher routeIDs; transpose `node_location`; rebuild `Node_related`; write `pm_now,pt_now,pc_now,pd_now` into `Penalty_all(2,:)`; `nbMoves++`; `searchCompleted=false`; update `whenLastModified` (`isdelete`: remove `routeU`, set `routeV=nbMoves`; else set both `routeU,routeV=nbMoves`); `isSuccess=true`.

---

### 5.2 m2 — "Insert u,x after v" (relocate the PAIR (u, its successor x) to after v, keeping order u then x)

**Move:** remove the consecutive pair `nodeU, nodeX` (where `nodeX = succ(nodeU)`, `suX = succ(nodeX)`) from the source route; insert as `… nodeV nodeU nodeX nodeY …`.

**Early rejects:** `if nodeU==nodeY or nodeV==nodeX → false`. `if nodeX==0 → false` (u must have a successor). If cross-route and `nodeX > nbClients` (x is an AFS) → false (cannot move an AFS as part of the pair). Same-route with `nodeX_loc==100` (x is the AFS) → false.

**Case A — same route:**
- `costOne = −dAll(nodeU+1,nodeX+1) − dAll(preU+1,nodeU+1) − dAll(nodeX+1,suX+1) + dAll(preU+1,suX+1)` (excise the pair, reconnect `preU→suX`; note the internal u–x edge is removed and re-added, so it nets out).
- `costTwo = +dAll(nodeU+1,nodeX+1) − dAll(nodeV+1,nodeY+1) + dAll(nodeV+1,nodeU+1) + dAll(nodeX+1,nodeY+1)` (insert `v→u→x→y`, restoring the internal u–x edge).
- Segment deltas use `∓2*everTime` (TWO customers move between segments). Removal: `time += costOne/speed − 2*everTime`; insertion: `time += costTwo/speed + 2*everTime`. Recompute penalties for `routeU` (single route). `pm_now=pm`.
- `delta` as usual; reject if `>−1e-6`.
- Pointers: `succ(preU)=suX; pred(nodeU)=nodeV; succ(nodeX)=nodeY; pred(suX)=preU; succ(nodeV)=nodeU; pred(nodeY)=nodeX; node_location(nodeU)=node_location(nodeX)=nodeV_loc`. (The u→x link is unchanged.)

**Case B — different routes (`nodeX <= nbClients` required):**
- `isdelete = 1` iff source route has exactly 3 nodes and has an AFS (`sum==3 && AFS`).
- **B1 — routeV HAS AFS:** same `costOne/costTwo` as Case A; apply removal to `routeU`, insertion to `routeV`; recompute both routes' penalties; handle `isdelete` (deleteAFS_delta, −2·dAll depot-AFS in delta). Pointers as in Case A but with `routeID(nodeU)=routeID(nodeX)=routeV`; if `isdelete`, `deleteAFS_node`.
- **B2 — routeV NO AFS (CRI):** result `[…, nodeV, nodeU, nodeX, AFS_new, nodeY, …]`.
  - `costTwo = +dAll(nodeU+1,nodeX+1) − dAll(nodeV+1,nodeY+1) + dAll(nodeV+1,nodeU+1) + dAll(nodeX+1,nbClients+2) + dAll(nbClients+2,nodeY+1)` (insert pair then new AFS before y).
  - Build `routeU_now = [prefix(from nodeV back), nodeU, nodeX, (N+1), suffix(from nodeV forward)]`; `get_pd_pt` → new routeV pre/su distances & times.
  - Recompute penalties; handle `isdelete`; `delta` (with −2·dAll depot-AFS if isdelete).
  - Append AFS at `node_addAfs=numel(node_location)+1`. Pointers: `succ(preU)=suX; pred(nodeU)=nodeV; succ(nodeU)=nodeX; succ(nodeX)=node_addAfs; pred(suX)=preU; succ(nodeV)=nodeU; pred(nodeY)=node_addAfs; pred(node_addAfs)=nodeX; succ(node_addAfs)=nodeY; routeID(nodeU)=routeID(nodeX)=routeID(node_addAfs)=routeV`. Set `node_location` along `routeU_now` (≤AFS→−1, ≥AFS→1, AFS→100). If `isdelete`, `deleteAFS_node`.

**Postamble:** same shape as m1 (note m2/m3 do NOT include the "drop single-node source route" block that m1's postamble has; they rely on the driver's patch and `isdelete`). Update `Penalty_all`, `nbMoves`, `searchCompleted`, `whenLastModified`, `isSuccess=true`.

---

### 5.3 m3 — "Insert x,u after v" (relocate the pair but REVERSE order: x then u)

Identical setup to m2 (`nodeX=succ(nodeU)`, `suX=succ(nodeX)`), but the pair is inserted **reversed**: `… nodeV nodeX nodeU nodeY …`.

**Early rejects:** `if nodeX==0 or nodeX==nodeV or nodeU==nodeY → false`. Cross-route with `nodeX>nbClients` → false. Same-route `nodeX_loc==100` → false.

**Case A — same route:**
- `costOne` identical to m2 (excise pair, reconnect `preU→suX`).
- `costTwo = +dAll(nodeU+1,nodeX+1) − dAll(nodeV+1,nodeY+1) + dAll(nodeV+1,nodeX+1) + dAll(nodeU+1,nodeY+1)` (insert `v→x→u→y`).
- `∓2*everTime` in time deltas as in m2. Recompute `routeU` penalties.
- Pointers: `succ(preU)=suX; pred(nodeU)=nodeX; succ(nodeU)=nodeY; pred(nodeX)=nodeV; succ(nodeX)=nodeU; pred(suX)=preU; succ(nodeV)=nodeX; pred(nodeY)=nodeU; node_location(nodeU)=node_location(nodeX)=nodeV_loc`. (Order reversed: v→x→u→y.)

**Case B — different routes (`nodeX<=nbClients`):**
- `isdelete = 1` iff `sum(routeID==routeU)==3 && AFS`.
- **B1 — routeV HAS AFS:** costs as Case A; apply, recompute both routes; handle isdelete; pointers as Case A with cross-route routeID; deleteAFS_node if isdelete.
- **B2 — routeV NO AFS (CRI):** result `[…, nodeV, nodeX, nodeU, AFS_new, nodeY, …]`.
  - `costTwo = +dAll(nodeU+1,nodeX+1) − dAll(nodeV+1,nodeY+1) + dAll(nodeV+1,nodeX+1) + dAll(nodeU+1,nbClients+2) + dAll(nbClients+2,nodeY+1)`.
  - `routeU_now = [prefix, nodeX, nodeU, (N+1), suffix]`; `get_pd_pt` → new routeV metrics.
  - Recompute penalties; isdelete handling; delta.
  - Append AFS; pointers: `succ(preU)=suX; pred(nodeU)=nodeX; succ(nodeU)=node_addAfs; pred(nodeX)=nodeV; succ(nodeX)=nodeU; pred(suX)=preU; succ(nodeV)=nodeX; pred(nodeY)=node_addAfs; pred(node_addAfs)=nodeU; succ(node_addAfs)=nodeY; routeID(nodeU)=routeID(nodeX)=routeID(node_addAfs)=routeV`; set `node_location` along `routeU_now`; deleteAFS_node if isdelete.

**Postamble:** same as m2.

---

### 5.4 m4 — "Swap u and v" (exchange the two single customers u and v)

**Move:** swap positions of `nodeU` and `nodeV` (each a single node). `preU,nodeX` around u; `preV,nodeY` around v.

**Early rejects:** `if nodeU==preV or nodeU==nodeY or nodeU>nodeV → false`. (The `nodeU>nodeV` guard makes each unordered pair considered once, avoiding duplicate/self swaps.)

**Cost:**
- `costOne = −dAll(preU+1,nodeU+1) − dAll(nodeU+1,nodeX+1) + dAll(preU+1,nodeV+1) + dAll(nodeV+1,nodeX+1)` (u's slot now holds v).
- `costTwo = −dAll(preV+1,nodeV+1) − dAll(nodeV+1,nodeY+1) + dAll(preV+1,nodeU+1) + dAll(nodeU+1,nodeY+1)` (v's slot now holds u).
- Segment deltas: for u's route, use `nodeU_loc` → `distance/time` of `routeU += costOne`; NOTE m4 uses **no `everTime` term** (a swap does not change how many customers are in each segment). For v's route, `nodeV_loc` → `+= costTwo`. This works whether same or different routes (if same route, both deltas hit `routeU=routeV`).
- Recompute `pd_v`, `pt_v` for `routeU` and `routeV`; `pc_now`; `pm_now=pm`.
- `delta = costOne+costTwo+Δpenalties`. Reject if `>−1e-6`. **No route deletion possible.**

**Pointers:** `succ(preU)=nodeV; pred(nodeU)=preV; succ(nodeU)=nodeY; pred(nodeX)=nodeV; succ(preV)=nodeU; pred(nodeV)=preU; succ(nodeV)=nodeX; pred(nodeY)=nodeU; routeID(nodeU)=routeV; routeID(nodeV)=routeU; node_location(nodeU)=nodeV_loc; node_location(nodeV)=nodeU_loc`.

**Postamble:** update `Route_related`, `Node_related`, `Penalty_all`, `nbMoves`, `searchCompleted`, `whenLastModified(routeU)=whenLastModified(routeV)=nbMoves`, `isSuccess=true`. (No route-deletion branch, no CRI.)

---

### 5.5 m5 — "Swap (u,x) pair with single v" (exchange the pair u,x for the single node v)

**Move:** the consecutive pair `nodeU, nodeX` (with `suX=succ(nodeX)`) trades places with the single node `nodeV`. Result: u's route position now holds `v`; v's position now holds `u,x`. `preU … nodeV … suX` on one side; `preV … nodeU nodeX … nodeY` on the other.

**Early rejects:** `if nodeU==preV or nodeX==preV or nodeU==nodeY or nodeX==0 → false`. Cross-route with `nodeX>nbClients` → false. Same-route with `nodeX_loc==100` (x is AFS) → false.

**Cost (all variants):**
- `costOne = −dAll(nodeU+1,nodeX+1) − dAll(preU+1,nodeU+1) − dAll(nodeX+1,suX+1) + dAll(preU+1,nodeV+1) + dAll(nodeV+1,suX+1)` (remove the pair, drop v into `preU … v … suX`).
- `costTwo (no new AFS) = +dAll(nodeU+1,nodeX+1) − dAll(preV+1,nodeV+1) − dAll(nodeV+1,nodeY+1) + dAll(preV+1,nodeU+1) + dAll(nodeX+1,nodeY+1)` (place `preV → u → x → y`, internal u–x edge restored).
- Segment time deltas use `∓everTime` (net one customer moves between the two routes/segments: two leave one side, one enters). Removal side (u's segment, per `nodeU_loc`): `time += costOne/speed − everTime`. Insertion side (v's route/segment, per `nodeV_loc`): `time += costTwo/speed + everTime`.

**Case A — same route (`routeU==routeV`), `nodeX_loc != 100`:** apply both deltas to `routeU`, recompute its penalties. Pointers: `succ(preU)=nodeV; pred(nodeU)=preV; succ(nodeX)=nodeY; pred(suX)=nodeV; succ(preV)=nodeU; pred(nodeV)=preU; succ(nodeV)=suX; pred(nodeY)=nodeX; routeID(nodeX)=routeID(nodeU)=routeV; routeID(nodeV)=routeU; node_location(nodeU)=node_location(nodeX)=nodeV_loc; node_location(nodeV)=nodeU_loc`.

**Case B — different routes (`nodeX<=nbClients`):**
- **B1 — routeV HAS AFS:** apply removal deltas to `routeU`, insertion deltas to `routeV`; recompute both; `delta` (NO route-deletion, NO −2·dAll term — m5 never deletes). Pointers as Case A (cross-route). `node_location` updates as Case A.
- **B2 — routeV NO AFS (CRI):** result: v goes into routeU's segment; `u,x,AFS_new` go into routeV.
  - `costTwo = +dAll(nodeU+1,nodeX+1) − dAll(preV+1,nodeV+1) − dAll(nodeV+1,nodeY+1) + dAll(preV+1,nodeU+1) + dAll(nodeX+1,nbClients+2) + dAll(nbClients+2,nodeY+1)`.
  - Rebuild `routeV_now`: collect routeV's node sequence (forward + backward from nodeV), then: set positions equal to `nodeX`→0, `nodeU`→0 (they are leaving via the pair-move? — actually in m5 u,x are ENTERING routeV; this bookkeeping removes any stale references and marks `nodeV`→−1 as the split anchor), then reinsert: `routeV_now = [ part up to the −1 anchor, nodeU, nodeX, (N+1), rest ]`, remove the −1 marker. Call `get_pd_pt` → new routeV `d1,d2,t1,t2`.
    *(Faithful C++: simplest is to construct the post-move routeV node list explicitly as `[…preV, nodeU, nodeX, AFS_new, nodeY…]` inserted at v's old location within routeV, then compute via §2.4. The MATLAB's `0/−1` marker juggling is just building that same list.)*
  - Apply removal deltas to `routeU` (v inserted there); set routeV metrics from `get_pd_pt`; recompute penalties; `delta` (no deletion term).
  - Append AFS at `node_addAfs=numel(node_location)+1`. Pointers: `succ(preU)=nodeV; pred(nodeU)=preV; succ(nodeX)=node_addAfs; pred(suX)=nodeV; succ(preV)=nodeU; pred(nodeV)=preU; succ(nodeV)=suX; pred(nodeY)=node_addAfs; pred(node_addAfs)=nodeX; succ(node_addAfs)=nodeY; routeID(nodeV)=routeU; routeID(nodeU)=routeID(nodeX)=routeID(node_addAfs)=routeV; node_location(nodeV)=nodeU_loc`; then set `node_location` along `routeV_now` (≤AFS→−1, ≥AFS→1, AFS→100).

**Postamble:** update all matrices, `Penalty_all`, `nbMoves`, `searchCompleted=false`, `whenLastModified(routeU)=whenLastModified(routeV)=nbMoves`, `isSuccess=true`. No route deletion.

---

## 6. DRIVER / CONTROL FLOW (`ELS_mian.m`)

The local search is a **first-improvement, granular, restart-until-stable** loop over customer pairs `(U, V)`.

### 6.1 Setup
- `afs_time_delay = Route_related(:,4)` (initial AFS arrival times).
- `correlatedVertices = vrp.correlatedVertices` (per-customer candidate `V` lists).
- Build initial `Node_related` from the chromosome routes `chromR` via `phrase_chromR` (§ below) plus the passed-in `node_location`. If not a repair pass, first strip depots from `chromR` by subtracting 1 and dropping zeros (`chromR_move{i} = chromR{i}−1; remove 0s`).
- `whenLastTestedRI = zeros(nbClients,1)`; `whenLastModified = zeros(numRoutes,1)`; `nbMoves=0`; `loopID=0`; `searchCompleted=false`.

`phrase_chromR(vrp,tspid,nbClients,chromR_move)`: from route lists, fill `predecessors/successor` (consecutive links within each route, boundary links to 0) and `routeID` (route index j for every node in `chromR{j}`).

### 6.2 Outer loop `while ~searchCompleted`
- If `loopID > 0`, set `searchCompleted = true` at the top (so it runs at least twice; a move setting `searchCompleted=false` forces another pass — this yields the standard "keep looping while improvements are found, with a mandatory second pass for empty-route moves").
- Re-seed RNG: `rng(SEED + tspid)`. For each customer `i`, with probability `1/nbGranular` (`mod(randi(...), nbGranular)==0`) randomly permute its `correlatedVertices(i,:)` row (granular neighbor shuffle).
- Time cutoff (only when `tspid != 1`): if `toc − sol_table.time(end-1) > sol_table.time(end-1)*20/numel(sol_table.ID)` then `break` (per-individual time budget).

### 6.3 Inner double loop over (U, V)
For `ii = 1..nbClients`: `nodeU = ii`; `correlatedU = correlatedVertices(nodeU,:)`; record `lastTestRINodeU = whenLastTestedRI(nodeU)`; set `whenLastTestedRI(nodeU) = nbMoves`.
  For each `jj` over `correlatedU`: `nodeV = correlatedU(jj)`.
  - **Consistency patch** (runs when `whenLastModified` changed since last check, i.e. after any accepted move): a maintenance block that (a) removes any route row whose `distance_pre` and `distance_su` are both ≈0 (dead route), (b) for each single-AFS route checks whether removing the AFS and directly connecting its neighbors keeps distance ≤ `V_Dmax`; if so it **drops the AFS** and recomputes that route's `pt/pd` via `get_pd_pt` and re-derives PC/PT, reordering AFS node rows and updating `routeID`/pred/succ renumbering. This is a periodic "conditional AFS REMOVAL" cleanup (the counterpart to CRI insertion). Port it as: after each accepted move, scan routes; for any AFS route where the shortcut distance (skip the AFS) ≤ V_Dmax, remove the AFS, recompute that route, and rebuild PC/PT. (Faithful behavior; the exact array-renumbering is bookkeeping.)
  - **Move gate:** attempt moves only if `loopID==0` OR `max(whenLastModified(routeID(nodeU)), whenLastModified(routeID(nodeV))) > lastTestRINodeU` (i.e. one of the involved routes changed since U was last tested — the "don't-re-test-unchanged" acceleration).
  - **Move sequence (first-improvement):** call in order `m1, m2, m3, m4, m5, m6, m7, m8, m9`; after each, `if isSuccess: continue` to the next `(U,V)`. (m6–m9 are additional operators not in scope here; keep the order.)
  - **Depot moves:** if `predecessors(nodeV)==0` (V is first on its route, adjacent to depot), additionally try `Depot_m1, Depot_m2, Depot_m3, Depot_m8, Depot_m9` (each with `if isSuccess: continue`).
  - **New-route moves:** after the full `V` loop for a given `U` (`nodeV == correlatedU(end)`) and if `loopID != 1`, try `NewRoute_m1, NewRoute_m2, NewRoute_m3` (open a brand-new route for U).
- After the `U` loop, `loopID++`.

### 6.4 Finalization
- Rebuild `chromR` from pred/succ/routeID via `get_chromR`.
- `distance_Total = Σ Route_related(:,[6,7])`; `cost_Total = distance_Total + Σ Penalty_all(2,:)`; `IsFeasible = (Σ Penalty_all(2,:)==0)`.
- Store all fields back into `sol_table` for this `tspid` (chromR, distances, times, per-route/total penalties, feasibility, etc.).

`get_chromR(pred,succ,routeID,nb)`: reconstruct each route by starting at every node whose predecessor is 0 (depot) and following successors until 0; assign the traversed sequence to `chromR{routeID(start)}`; drop empty routes.

### 6.5 Termination
The loop ends when a full pass produces no accepted move (`searchCompleted` stays `true`) — i.e. a local optimum under the granular neighborhood — or the per-individual time cutoff `break`s. This is deterministic given `SEED+tspid`.

---

## 7. PENALTY WEIGHT ADAPTATION (context; happens OUTSIDE m1–m5)
Operators read weights `wt,wc,wd,wm` from `Penalty_all(1,:)` and never change them. The genetic outer loop (population management, not in these files) adapts weights based on the fraction of feasible offspring vs `par_hgs.targetFeasible` (standard HGS penalty adaptation: raise a weight if too few solutions satisfy that constraint, lower it if too many). For a faithful port of just the local search, treat weights as constant inputs; expose the four weights so the outer loop can adjust them between LS invocations. The initial `pm` is `wm*max(R−V_nb,0)`.

---

## 8. IMPLEMENTATION CHECKLIST / GOTCHAS
1. Distance lookups: `dAll(a+1,b+1)` → `distance[a][b]`; AFS location value = `nbClients+1` (table index `nbClients+2`); depot value 0 (table index 1).
2. Acceptance: `delta <= −1e-6` exactly; first-improvement (stop at first accepting operator).
3. PD penalizes pre- and post-AFS segments SEPARATELY vs `V_Dmax`.
4. PT uses per-route `time_v` vs `T_max`; refuel time `T_Afs` is NOT in `time_v`.
5. PC = over-capacity overlap area at the single shared AFS, computed by `get_pc_now` after the delay-scheduler resolves conflicts; recompute for the WHOLE solution (all AFS routes) on every candidate.
6. PM = `wm*max(numRoutes − V_nb, 0)`; only route creation/deletion changes it.
7. CRI: inserting into an AFS-less target route ALWAYS adds one new AFS immediately after the inserted customer(s); inserting into an AFS-bearing route adds none.
8. Route deletion (only m1/m2/m3, cross-route, source route reduced to just-its-AFS): remove the route, `pm -= wm` (floored at 0), and add `−2*dAll(1,nbClients+2)` to `delta`. Using full distance totals for `delta` captures this automatically.
9. m4/m5 never delete a route and m4 has no CRI branch; m5 has a CRI branch but no deletion.
10. `everTime` bookkeeping in segment-time deltas (`±everTime` for single-node moves m1/m5-cross, `±2*everTime` for pair moves m2/m3, none for m4) tracks how many customers cross between the pre/post segments — but a full re-evaluation via §1.6/§2.4 makes these adjustments unnecessary; they will fall out correctly.
11. The `whenLastModified` / `whenLastTestedRI` machinery is a pruning optimization; you may implement it faithfully for identical search trajectories, or re-test all pairs (slower but same local optimum set under first-improvement is not guaranteed identical, so replicate it if bit-exact behavior matters).
12. `node_addAfs = numel(node_location)+1` — new AFS nodes are appended; the driver's consistency patch later re-sorts AFS rows into route order. In C++ prefer per-route AFS ownership objects and rebuild the flat arrays only at boundaries.

---

Source files (all under `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm`):
- Operators: `Efficient local search/m1.m`, `m2.m`, `m3.m`, `m4.m`, `m5.m`
- Helpers: `Efficient local search/AFSdelay_new.m`, `AFSdelay_recursion.m`, `get_pc_now.m`, `get_pd_pt.m`, `deleteAFS_delta.m`, `deleteAFS_node.m`, `get_chromR.m`, `phrase_chromR.m`
- Driver: `Efficient local search/ELS_mian.m`
- Init/objective (data-structure construction, penalties, node_location): `Load/chromR_detail_all.m`, `Load/get_vrp.m`