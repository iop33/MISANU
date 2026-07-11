I now have a complete and detailed understanding of the entire algorithm. I have everything required to write the exhaustive spec covering all requested sections, with m6–m9 as the primary focus. Here is the specification.

---

# SPECIFICATION: METS local-search operators m6–m9 (GrVRP-PCAFS)

This spec is implementation-ready. In C++ we FULLY RE-EVALUATE each candidate move's cost/feasibility from the resulting route(s) — you do NOT need to replicate the MATLAB incremental delta bookkeeping. The delta formulas below are given only so you can reproduce the exact acceptance decision if you choose to (recommended: they are what MATLAB uses). Where they matter for correctness (tie-break, penalty adaptation), they are called out.

Path context: all MATLAB under `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/`. Operators live in `Efficient local search/m6.m … m9.m`; helpers `AFSdelay_new.m`, `AFSdelay_recursion.m`, `get_pc_now.m`, `get_pd_pt.m`, `deleteAFS_delta.m`, `deleteAFS_node.m`; driver `ELS_mian.m`; builder `Load/chromR_detail_all.m`; main `Main_METS.m`.

============================================================
## 1. DATA STRUCTURES
============================================================

### 1.1 Indexing convention (CRITICAL — MATLAB 1-based)
There are three distinct numbering spaces. Keep them separate in C++.

- **Node IDs used in linked-list arrays** (`predecessors`, `successor`, `routeID`, `node_location` indices, and the values of `nodeU`, `nodeV`, `nodeX`, etc.):
  - Depot is represented as **0** (a node whose predecessor/successor value is 0 means "depot / route end").
  - Customer `c` (1..nbClients) has node ID `c`.
  - AFS (refueling-station) nodes have IDs `nbClients+1, nbClients+2, …`. Each active route that uses an AFS owns exactly one AFS node; AFS node IDs are appended AFTER all customers in the arrays.
  - So the arrays `predecessors/successor/routeID/node_location` are indexed 1..(nbClients + nbAFSnodes). Row `k` for `k<=nbClients` = customer k; row `k>nbClients` = an AFS node.

- **Distance-table indices** (`dAll = vrp.distance_table`, `dAll(a+1,b+1)`):
  - The distance table is 1-based with **row/col 1 = depot**. A node with node-ID `n` maps to distance-table index `n+1`. Thus `dAll(nodeU+1, nodeX+1)` = distance from node U to node X; `dAll(1, k+1)` = depot→node k; `dAll(1, nbClients+2)` = depot→(the "canonical" AFS at node-ID `nbClients+1`, whose table index is `nbClients+2`). The value `dAll(1,nbClients+2)` is the depot↔AFS one-way distance used in AFS-add/remove deltas.
  - **The distance table is symmetric** and holds real (unrounded) distances during the search. (Main_METS rounds to 2 decimals only at the very end, after the run.)

- **chromR representation** (`chromR{i}`): a 1-based route as a row vector `[1 … 1]` where **1 = depot** and other entries are `nodeID+1` (customer c is `c+1`, AFS node n is `n+1`). i.e. chromR uses distance-table indexing. `chromR_move{i}` is the same route with `-1` applied and depot entries (value 0 after -1) removed, i.e. it lists node IDs (customers as c, AFS as `nbClients+…`) WITHOUT depot. Conversion between the two spaces happens in `get_chromR`/`phrase_chromR`/`chromR_detail_all`.

### 1.2 Route / node arrays (the state m6–m9 operate on)
Passed in packed matrices; unpacked at the top of every operator.

**`Node_related`** = `[predecessors, successor, routeID, node_location]`, size `(nbNodes × 4)`, indexed by node ID:
- `predecessors(n)`: node ID of n's predecessor in its route, or `0` if n is first (depot is the implicit start).
- `successor(n)`: node ID of n's successor, or `0` if n is last (route returns to depot).
- `routeID(n)`: 1-based index of the route containing n. Routes are numbered `1..nbRoutes`. Row `r` of `Route_related` describes route `r`.
- `node_location(n)`: position tag relative to the route's AFS:
  - `-1` = node is **before** the AFS (or route has no AFS at all → whole route is "pre").
  - `1`  = node is **after** the AFS.
  - `100` = node **is** the AFS itself.
  Meaning: distance/time of a route is split into a "pre-AFS" part (d1) and "post-AFS" part (d2); `node_location` says which part a node belongs to. A route with no AFS has all nodes at `-1` and only d1 nonzero.

**`Route_related`** = `[pd_v(:,1), pd_v(:,2), pt_v, time_afs, time_v, distance_pre, distance_su]`, size `(nbRoutes × 7)`, row = route:
- col1 `pd_v(:,1)` = distance penalty of the pre-AFS leg = `max(distance_pre - V_Dmax, 0)`.
- col2 `pd_v(:,2)` = distance penalty of the post-AFS leg = `max(distance_su - V_Dmax, 0)`.
- col3 `pt_v` = time penalty of the route = `max(time_afs + time_su - T_max_V, 0)` = `max(time_v - T_max_V, 0)`.
- col4 `time_afs` = **arrival time at the AFS** (i.e. cumulative travel+service time from depot to the AFS, excluding waiting). Also serves as the duration of the pre-AFS leg (d1/speed + service). For a route with no AFS this equals the whole route time.
- col5 `time_v` = **total route duration** = `time_afs + time_su`.
- col6 `distance_pre` (d1) = travel distance of the pre-AFS leg (depot → AFS), or whole-route distance if no AFS.
- col7 `distance_su` (d2) = travel distance of the post-AFS leg (AFS → depot); `0` if no AFS.
- Derived in every operator: `time_su = time_v - time_afs` (post-AFS leg duration).

**`Penalty_all`** = `[wt, wc, wd, wm ; pt, pc, pd, pm]`, size `2×4`:
- Row 1 = penalty **weights** (coefficients): `wt` (time/Tmax), `wc` (AFS capacity), `wd` (distance/Dmax), `wm` (vehicle count).
- Row 2 = current **weighted penalty totals**: `pt`, `pc`, `pd`, `pm`.

**`afs_time_delay`** = column vector, size `(nbRoutes × 1)` (one entry per route's AFS; `0` if the route has no AFS). Holds the AFS arrival time possibly *shifted later* to resolve station-capacity conflicts (see §3.2). Recomputed by `AFSdelay_new`.

**`whenLastModified`** = `(nbRoutes × 1)`: last `nbMoves` value at which each route changed (used to skip untouched pairs).

**Scalars threaded through:** `nbMoves` (move counter), `searchCompleted` (bool), `isSuccess` (bool return), `yu`/`yv` (= routeU/routeV of the attempted move, always returned), `isdelete` (m8 only: 1 if a route got emptied and deleted).

### 1.3 Instance data (`vrp`, read-only)
- `vrp.nb_customer` = number of customers (`nbClients`).
- `vrp.T_Customer` = service time per customer (`everTime`).
- `vrp.V_speed` = travel speed (`speed`); time = distance/speed + service.
- `vrp.V_Dmax` = max travel distance per **leg** (pre or post AFS) (`V_Dmax`).
- `vrp.T_max_V` = max route duration (`T_max`).
- `vrp.T_Afs` = refueling duration at an AFS (`T_Afs`).
- `vrp.C_Afs` = AFS capacity = number of vehicles that can refuel simultaneously (`C_Afs`).
- `vrp.V_nb` = free (non-penalized) vehicle count.
- `vrp.distance_table` = symmetric matrix, index 1 = depot (`dAll`).
- `vrp.correlatedVertices` = `(nbClients × nbGranular)` neighbor lists (granular search order).
- `vrp.last_customer`, `vrp.last_F_location`, `vrp.T_Start` = used by split/broken-pair distance.

============================================================
## 2. OBJECTIVE & PENALTIES (exact formulas)
============================================================

### 2.1 Per-route quantities
For a route with an AFS, split at the AFS:
- `distance_pre` = sum of edge distances depot→…→AFS.
- `distance_su` = sum of edge distances AFS→…→depot.
- `time_afs` = `distance_pre/speed + (#customers strictly before AFS)*everTime`  (arrival at AFS; AFS itself and the depot add no service time).
- `time_su` = `distance_su/speed + (#customers after AFS)*everTime`.
- `time_v` = `time_afs + time_su`. **Refueling time `T_Afs` and any AFS waiting time are NOT added into `time_v` here** — station conflicts are handled separately by the PC penalty (§3.2), not by inflating route time.
For a route WITHOUT an AFS: `distance_pre` = whole route distance, `distance_su = 0`, `time_afs = time_v = distance_pre/speed + (#customers)*everTime`, `time_su = 0`.

(Reference implementation: `get_pd_pt(vrp, route_now, everTime, speed)` returns `[d1,d2,t1,t2]`. `route_now` is a node-ID list WITHOUT depot; if its max element `> nb_customer` it contains an AFS at position `n` = index of the max element, and it splits there; else no AFS. See §4.6.)

### 2.2 The four penalties (weighted)
- **PD (distance / Dmax):** per leg, `pd_v(r,1)=max(distance_pre(r)-V_Dmax,0)`, `pd_v(r,2)=max(distance_su(r)-V_Dmax,0)`. Total `pd = wd * sum(sum(pd_v))` over all routes and both legs.
- **PT (time / Tmax):** `pt_v(r)=max(time_afs(r)+time_su(r)-T_max_V,0)`. Total `pt = wt * sum(pt_v)`.
- **PM (vehicle count):** `penalty_m = max(nbRoutes - vrp.V_nb, 0)`; `pm = wm * penalty_m`. (`wm` = `Penalty_all(1,4)`, initialized to 0 in Main_METS, so PM is effectively 0 unless set; operators still carry it.)
- **PC (AFS capacity):** `pc = wc * pc_now`, where `pc_now` = total (over time) vehicle-overlap beyond capacity at the shared AFS resource; computed by `AFSdelay_new` (§3.2).

### 2.3 Total cost & feasibility
- `distance_Total = sum(distance_pre + distance_su)` over all routes.
- `cost_Total = distance_Total + pt + pc + pd + pm`  (weighted penalties added).
- **Feasible** iff `pt + pc + pd + pm == 0` (all weighted penalty totals zero). `IsFeasible=1` else `0`.

### 2.4 Move acceptance (all operators identical rule)
Every candidate computes a `delta` = (change in raw distance) + (change in each weighted penalty), where the "change" of a penalty is `new_weighted_total - old_weighted_total`:
```
delta = costOne + costTwo
      + (pm_now - pm) + (pt_now - pt) + (pc_now - pc) + (pd_now - pd)
      [ - AFS-removal distance credit, if a route/AFS was deleted ]  // see per-operator
Accept iff delta <= -1e-6  (MATLAB: `delta > -0.000001` → reject).
```
`costOne`, `costTwo` are the raw edge-distance deltas of the move (negative for removed edges, positive for added). `pt`,`pc`,`pd`,`pm` are the pre-move weighted totals from `Penalty_all(2,:)`; `pt_now`,… are recomputed post-move weighted totals. On accept, the operator commits linked-list changes and writes the new `Penalty_all(2,:) = [pt_now,pc_now,pd_now,pm_now]`, increments `nbMoves`, sets `searchCompleted=false`, stamps `whenLastModified`.

**C++ FULL RE-EVALUATION strategy (recommended):** Rather than tracking d1/d2/time deltas, for each candidate build the resulting route(s) as node-ID lists, call the equivalent of `get_pd_pt` on each changed route to get `(d1,d2,t1,t2)`, rebuild `pd_v/pt_v/time_afs/time_v/distance_pre/distance_su` for the changed routes, recompute `pc_now` via the AFS-delay routine over ALL routes, form `pt_now/pd_now/pm_now` totals, then compute `delta` with the SAME formula and the SAME `-1e-6` threshold. **You must keep costOne/costTwo (raw distance delta) in the delta even though distance_Total also changes — MATLAB's `delta` uses the raw edge delta for the distance component and the penalty *differences* for the penalties; they are consistent because `pd` is a penalty not raw distance.** Concretely: `delta = (new_distance_Total - old_distance_Total) + (pt_now-pt)+(pc_now-pc)+(pd_now-pd)+(pm_now-pm)`. This equals the MATLAB `costOne+costTwo+...` because `costOne+costTwo` is exactly the raw distance change. Use the deletion credit terms below when an AFS/route is removed.

### 2.5 AFS insertion/removal distance accounting (Conditional AFS Insertion)
When a move ADDS a new AFS to a route (needed because a segment moved into an AFS-less route that must be split), the added edges go via node-ID `nbClients+1` (table index `nbClients+2`); the pre/post split distances are recomputed by `get_pd_pt`. When a move REMOVES an AFS/route entirely, MATLAB subtracts a distance credit in the delta:
- Remove one AFS-only route (both depot↔AFS edges gone): `delta -= 2*dAll(1,nbClients+2)`.
- Merge that spares one depot↔AFS edge: `delta -= dAll(1,nbClients+2)` (m8 branch "xx≠0, both non-AFS" with isdelete).
`deleteAFS_delta` also does `pm_now = max(pm_now - wm, 0)` (one fewer vehicle) and drops that route's rows; `deleteAFS_node` renumbers `routeID`/`predecessors`/`successor` and removes the AFS node (see §4.7–4.8).

============================================================
## 3. SHARED HELPERS (needed by m6–m9)
============================================================

### 3.1 `get_pd_pt(vrp, route_now, everTime, speed) -> [d1,d2,t1,t2]` (§4.6)
Recomputes a single route's split distances/times from a node-ID list (no depot). This is the "full re-evaluation" primitive; used by m7 and m9 when a reversal makes incremental math impossible.

### 3.2 AFS capacity penalty — `AFSdelay_new` / `AFSdelay_recursion` / `get_pc_now`
Signature: `[pc_now, afs_time_delay] = AFSdelay_new(time_afs, time_v, T_max, T_Afs, C_Afs, afs_time_delay)`.
- `time_afs`, `time_v` are the per-route vectors (only the AFS-owning routes; callers pass `time_afs(1:nbRoutesWithAFS)` where relevant — in m6/m7/m9 they pass the full vectors, in m8/m1 sometimes `time_afs(1:numel(node_location)-nbClients)` = the AFS rows).
- If `numel(time_afs)==0` → `pc_now=0`.
- `time_V_shifting(r) = max(T_max - time_v(r), 0)` = slack by which route r's AFS visit may be delayed without violating Tmax.
- `afs_time_delay := time_afs` initially, then a **greedy backtracking scheduler** (`AFSdelay_recursion`) delays overlapping AFS visits within available slack to minimize overlap-capacity penalty. Semantics:
  - Two vehicles i, ii conflict if their `[afs_time_delay, afs_time_delay+T_Afs]` intervals overlap (`t1+T_Afs > t2` where t1≤t2 are the sorted delays).
  - If neither can be delayed enough (`t1+T_Afs-t2 > slack_of_2` AND `t2+T_Afs-t1 > slack_of_1`) → unavoidable conflict (backtracking tries both delay directions, keeps the lower `pc_now`).
  - Else delay the one that resolves it, consuming its slack; restart the scan.
- Final `pc_now` from `get_pc_now(afs_time_delay, T_max, time_v, time_afs, T_Afs, C_Afs)`:
  - Drop zero entries from `afs_time_delay`; `afs_time_end = afs_time_delay + T_Afs`.
  - Sort all start/end times into `c`; `time_during = diff(c)` (interval lengths).
  - For each interval, count vehicles refueling: `overlap = (afs_time_delay - c(1:end-1))` with `eps < d+T_Afs & d <= 0`; `nb_fueling = sum(overlap)`.
  - `nb_fueling = max(nb_fueling - C_Afs, 0)` (overload beyond capacity).
  - `pc_now = max(nb_fueling · time_during', 0)` = total overloaded vehicle-time.
- **C++ note:** this is a self-contained scheduling+overlap routine. Re-implement faithfully; it is called once per candidate to get `pc_now`. The exact backtracking is complex but only affects PC; for a faithful port replicate `AFSdelay_recursion` + `get_pc_now`. The initial (non-recursive) conflict resolution used at solution-build time is in `chromR_detail_all` lines 136–221 and matches this logic (with the extra `q1>q2` tie-break branch when both can be delayed).

### 3.3 `deleteAFS_delta` / `deleteAFS_node` (§4.7–4.8) — route deletion bookkeeping.

============================================================
## 4. EACH FUNCTION — STEP-BY-STEP
============================================================

Common prologue for m6–m9 (identical unpacking):
```
nbClients=vrp.nb_customer; everTime=vrp.T_Customer; speed=vrp.V_speed;
V_Dmax=vrp.V_Dmax; T_max=vrp.T_max_V; T_Afs=vrp.T_Afs; C_Afs=vrp.C_Afs;
wt,wc,wd,wm = Penalty_all(1,1..4);   pt,pc,pd,pm = Penalty_all(2,1),(2,2),(2,3),(2,4);
pd_v=Route_related(:,[1,2]); pt_v=Route_related(:,3);
dAll=vrp.distance_table;
distance_pre=Route_related(:,6); distance_su=Route_related(:,7);
time_v=Route_related(:,5); time_afs=Route_related(:,4); time_su=time_v-time_afs;
predecessors,successor,routeID = Node_related(:,1..3); node_location=Node_related(:,4)';
// nodeU info:
nodeU_loc=node_location(nodeU); routeU=routeID(nodeU); preU=predecessors(nodeU); nodeX=successor(nodeU);
if nodeX: suX=successor(nodeX); nodeX_loc=node_location(nodeX);
// nodeV info:
nodeV_loc=node_location(nodeV); routeV=routeID(nodeV); nodeY=successor(nodeV); preV=predecessors(nodeV);
if nodeY: suY=successor(nodeY); nodeY_loc=node_location(nodeY);
yu=routeU; yv=routeV;
```
`nodeX` = node after U; `nodeY` = node after V; `suX`/`suY` = nodes after X/Y; `preU`/`preV` = nodes before U/V. Value 0 = depot/none.

Common epilogue (all operators, only reached on ACCEPT):
```
Route_related = [pd_v, pt_v, time_afs, time_v, distance_pre, distance_su];   // repack
Node_related  = [predecessors, successor, routeID, node_location'];
Penalty_all(2,:) = [pt_now, pc_now, pd_now, pm_now];   // note col order: (2,1)=pt,(2,2)=pc,(2,3)=pd,(2,4)=pm
nbMoves++; searchCompleted=false;
whenLastModified(routeU)=nbMoves; whenLastModified(routeV)=nbMoves;   // (m8 differs if isdelete)
isSuccess=true;
```

---
### 4.1 m6 — SWAP (exchange node U and node V, each with its successor)
"Swap u,x and v,y": swap the pair (U, its successor X) with the pair (V, its successor Y) between/within routes. Effectively U↔V positions AND X↔Y positions.

**Neighborhood / candidate generation:** driven by ELS main loop: `nodeU` iterates 1..nbClients, `nodeV` iterates over `correlatedVertices(nodeU,:)`. m6 is tried after m1–m5 fail.

**Guard rejections (return isSuccess=false):**
1. Any of: `nodeX==0` OR `nodeY==0` OR `nodeY==preU` OR `nodeU==nodeY` OR `nodeX==nodeV` OR `nodeV==suX`.
2. If both `nodeX_loc==100` and `nodeY_loc==100` (both successors are AFS) → reject.
3. If **exactly one** of `nodeX_loc==100 / nodeY_loc==100` (the `else` at line 130) → reject. So m6 only proceeds when **neither X nor Y is an AFS** (`nodeX_loc~=100 && nodeY_loc~=100`).

**Move semantics (when accepted):** Replace edges around U..X and V..Y so that V takes U's slot and U takes V's slot; the successor chains splice X after Y's old position and Y after X's old position. Precisely the pointer updates (lines 102–129):
```
successor(preU)=nodeV; predecessors(nodeU)=preV;
successor(nodeX)=suY;  predecessors(suX)=nodeY;
successor(preV)=nodeU; predecessors(nodeV)=preU;
successor(nodeY)=suX;  predecessors(suY)=nodeX;
routeID(nodeX)=routeV; routeID(nodeU)=routeV; routeID(nodeV)=routeU; routeID(nodeY)=routeU;
node_location(nodeU)=nodeV_loc; node_location(nodeX)=nodeV_loc;
node_location(nodeV)=nodeU_loc; node_location(nodeY)=nodeU_loc;
```
Net effect: U (and X) move into route V at V's location tag; V (and Y) move into route U at U's location tag. Works intra- or inter-route.

**Cost (raw distance delta):**
```
costOne = -d(U,X) -d(preU,U) -d(X,suX) + d(V,Y) + d(preU,V) + d(Y,suX)
costTwo = -d(V,Y) -d(preV,V) -d(Y,suY) + d(U,X) + d(preV,U) + d(X,suY)
```
(all `d(a,b)=dAll(a+1,b+1)`.)

**Feasibility/cost re-eval:** Update pre/post distance & time of routeU and routeV depending on `nodeU_loc`/`nodeV_loc` (∈{-1,1}) — costOne applied to routeU's pre-leg if `nodeU_loc==-1` else su-leg; costTwo to routeV similarly. Recompute `pd_v` (both legs, both routes), `pt_v` (both routes), `time_v=time_afs+time_su`, then `pc_now=wc*AFSdelay_new(time_afs,time_v,…)`, `pm_now=pm`. Compute delta (§2.4). Accept iff `delta<=-1e-6`.

**AFS:** m6 never inserts/removes an AFS (guaranteed by rejecting AFS successors); node_location tags simply follow U/V.

---
### 4.2 m7 — 2-opt within one route: `(u,x)(v,y) → (u,v)(x,y)`
Reverse the sub-path between X and V (inclusive) inside a single route. U..V..Y with segment X…V reversed.

**Guards:**
- `routeU != routeV` → reject (line 8; both U,V must be same route).
- Walk predecessors from U back through the route; if V is reachable going backward (V is before U) → reject (lines 65–75). So V must be *after* U in the route.
- `nodeX==nodeV` → reject (adjacent, nothing to reverse).

**Segment `xxvv`:** the node-ID list from U forward up to and including V, then U dropped: `xxvv = [X, …, V]` (the sub-path to be reversed), built by following successors from U until V (lines 80–89).

**Two cases by AFS presence in the reversed segment:**
- **Case A — no AFS in `xxvv`** (`max(node_location(xxvv)) ~= 100`):
  - `costOne = -d(U,X) -d(V,Y)`, `costTwo = +d(U,V) + d(X,Y)`.
  - Apply `costOne+costTwo` to routeU's pre-leg if `nodeX_loc==-1` else su-leg if `nodeU_loc==1` (lines 94–100). Recompute pd_v/pt_v/pc/time. Delta, accept iff ≤ −1e-6.
  - Commit: reverse the segment pointers — `successor(U)=V`; `predecessors(Y)=X` (if Y); reverse predecessors/successors along `xxvv` (lines 118–127): `predecessors(xxvv(end))=U`; for i=1..end-1 `predecessors(xxvv(i))=xxvv(i+1)`; `successor(xxvv(1))=Y`; for i=2..end `successor(xxvv(i))=xxvv(i-1)`.
- **Case B — AFS inside `xxvv`** (`max==100`): the reversal crosses the AFS, so pre/post split must be fully recomputed.
  - Build the full ordered route `routeU_now` (walk successors from V forward → `routeU_NOW`; walk predecessors from V backward → `routeU_now`, flip, concatenate). Then reverse the sub-array from `nodeX` to `nodeV` within it (lines 132–153).
  - `[d1,d2,t1,t2] = get_pd_pt(vrp, routeU_now, everTime, speed)`; set `distance_pre/su(routeU)=d1/d2`, `time_afs/su(routeU)=t1/t2`. Recompute pd_v/pt_v/pc. Delta, accept.
  - Commit: same segment-reversal pointer updates as Case A, PLUS retag node_location around the (now relocated) AFS: everything from the new AFS position onward = `-1`? — precisely lines 185–187: `node_location(xxvv(pos(maxRoute):end)) = -1`; `node_location(xxvv(1:pos(maxRoute))) = 1`; `node_location(max(routeU_now)) = 100` where `max(routeU_now)` = the AFS node ID.

**C++ note:** In both cases you can simply rebuild `routeU_now` (reversed segment) as a node list and call `get_pd_pt` for the single route, then recompute penalties over all routes. That reproduces Case B exactly and is also correct for Case A.

---
### 4.3 m8 — 2-opt* (inter-route): `(u,x)(v,y) → (u,v)(x,y)`
Cross two DIFFERENT routes: tail of routeU after U is swapped with tail of routeV up to V. `xx` = nodes after U in routeU (the tail to move), `vv` = nodes from V backward to route start (the head of routeV to move).

**Guards:** `routeV==routeU` → reject (must be different routes).

**Segments:**
- `xx` = successors from U forward (tail of routeU after U); `d_xx` = its internal edge-length sum. `xx` excludes U (drops first).
- `vv` = predecessors from V backward incl. V (`[V, preV, …, start]`); `d_vv` = its internal edge sum.

**Branch structure (by whether `xx` empty, and AFS presence). All produce a valid crossed solution or reject.** The AFS-presence checks enforce that at most one AFS ends up per resulting route (Conditional AFS Insertion / consistency):

- **`xx` empty (U is last real node of routeU):**
  - Both routeU-region and `vv` AFS-free, or both have AFS → reject.
  - Exactly routeU has AFS, `vv` AFS-free (lines 95–160): possibly `isdelete=1` if routeV becomes AFS-only after the move (`nodeY && vv==all of routeV minus AFS && node_location(nodeY)==100`). Costs:
    ```
    costOne = -d_vv - d(depot, vv_end) - d(V,Y) + d(depot,Y)
    costTwo = +d_vv + d(U,V) + d(vv_end,depot) - d(U,depot)
    ```
    Apply costOne to routeV pre-leg (adjust time_afs by `-numel(vv)*everTime`), costTwo to routeU su-leg (`+numel(vv)*everTime`). Recompute pd/pt/pc. If isdelete: call `deleteAFS_delta` (drop routeV rows, `pm_now=max(pm_now-wm,0)`), recompute pc/pt/pd, and `delta -= 2*dAll(1,nbClients+2)`. Accept iff ≤ −1e-6. Commit pointer splice (lines 135–160) moving `vv` into routeU and reversing pointer order; retag `node_location(vv)=nodeU_loc`; if isdelete call `deleteAFS_node`.
  - routeU AFS-free & `vv` has AFS → reject.
- **`xx` non-empty:**
  - `xx` AFS-free & `vv` has AFS → reject; and `xx` has AFS & `vv` AFS-free → reject (asymmetry not allowed).
  - **Both `xx` and `vv` have AFS** (lines 172–274): may `isdelete` if `nodeY==0 && numel(xx)==1 && node_location(xx)==100`. Costs:
    ```
    costOne = -d_vv -d(depot,vv_end) -d(V,Y) + d_xx + d(X,Y) + d(depot,xx_end)
    costTwo = +d_vv +d(U,V) + d(vv_end,depot) - d_xx - d(U,X) - d(xx_end,depot)
    ```
    Full recompute of both routes' split distances/times using the saved d1/d2/t1/t2 of both routes and the segment sizes (lines 181–208), then pd/pt; `pc_now` from `AFSdelay_new(time_afs(1:nbAFSrows), time_v(1:nbAFSrows), …)`. isdelete path subtracts `2*dAll(1,nbClients+2)`. Commit: splice `vv`→routeU, `xx`→routeV, flip node_location signs (`node_location(vv)=-node_location(vv)` etc., then `-100→100`), and **swap the two AFS nodes' linked-list positions** (`xxafs=max(xx)`, `vvafs=max(vv)`; exchange their pre/su/routeID — lines 256–271) so each AFS stays with its route.
  - **Both `xx` and `vv` AFS-free** (lines 275–340): simplest cross.
    ```
    costOne (same as above), costTwo (same as above)
    ```
    Apply costTwo to routeU pre/su-leg (by `nodeX_loc==-1`/`nodeU_loc==1`) with service-time correction `(numel(vv)-numel(xx))*everTime`; costOne to routeV pre/su-leg with `(numel(xx)-numel(vv))*everTime`. Recompute pd/pt/pc. Delta (isdelete would subtract `dAll(1,nbClients+2)` but isdelete stays 0 here). Commit splice: `vv`→routeU, `xx`→routeV; `node_location(vv)=nodeU_loc`, `node_location(xx)=nodeV_loc`.

**Deletion / whenLastModified (m8 epilogue lines 358–370):** if `isdelete==1` and `delete_idx==routeV`: `whenLastModified(routeU)=nbMoves`, remove `whenLastModified(routeV)`; if `delete_idx==routeU`: remove `whenLastModified(routeU)`, set `whenLastModified(routeV)`. Else both set.

---
### 4.4 m9 — 2-opt* (inter-route, other orientation): `(u,x)(v,y) → (u,y)(x,v)`
Also crosses two different routes, but connects U→Y (V's successor tail) and X→V. Uses three segments: `xx` (routeU tail after U), `vv` (routeV head, predecessors from V back), `yy` (routeV tail after V, successors from V forward).

**Guards:** `routeV==routeU` → reject.

**Segments:** `xx` (successors from U, drop U; `d_xx`), `vv` (predecessors from V incl. V; `d_vv`), `yy` (successors from V, drop V; `d_yy`).

**Four top-level cases by (`xx` empty?, `yy` empty?), each with AFS sub-cases; each computes costOne/costTwo, re-evaluates penalties, checks `delta<=-1e-6`, and on accept splices pointers + retags node_location.** Key structure:
- **`xx` empty & `yy` empty:** only valid if `vv` has AFS (`max(node_location(vv))==100`); it simply *reverses routeV's AFS split* (swap distance_pre↔su, adjust times by ±everTime — lines 103–108) because U now attaches to Y (=0/depot) and V-chain flips. costOne=costTwo=0. Else reject.
- **`xx`≠0 & `yy` empty:** valid when exactly one of `xx`/`vv` has AFS.
  - `xx` has AFS, `vv` AFS-free (lines 162–229): `costOne = -d_xx -d(U,X) -d(xx_end,depot) + d(U,depot)`; `costTwo = +d_xx + d(X,V) + d(depot,xx_end) - d(V,depot)`. Full recompute of routeU (becomes AFS-free, `distance_su=0`,`time_su=0`) and routeV (absorbs `xx` reversed after V). node_location: `vv→1`, `xx` sign-flipped, `-100→100`.
  - `xx` AFS-free, `vv` has AFS (lines 230–301): symmetric; costs same formula; routeV split reversed; `node_location(xx)=-1`, `vv` flipped.
- **`xx` empty & `yy`≠0:** valid when exactly one of (routeU region)/`yy` has AFS.
  - routeU has AFS, `yy` AFS-free (lines 309–387): `costOne = -d(U,depot) + d_yy + d(U,Y) + d(yy_end,depot)`; `costTwo = -d_yy -d(V,Y) -d(yy_end,depot) + d(V,depot)`. Attach `yy` after U (into routeU); routeV keeps `vv`. Branch on whether `vv` has AFS for how routeV's split updates (lines 328–337). `node_location(yy)=nodeU_loc`.
- **`xx`≠0 & `yy`≠0:** both tails non-empty.
  - **Both AFS present** (`max(node_location(xx))==100 && max(node_location(yy))==100`, lines 393–482): rebuild both routes fully. `uu` = predecessors from U back (routeU head incl. U), flipped. `a=[uu, yy]` (new routeU), `b=[fliplr(xx), vv]` (new routeV). `[d1u,d2u,t1u,t2u]=get_pd_pt(vrp,a,…)`, likewise for `b`. Set both routes' distance_pre/su/time_afs/su. costOne/costTwo per lines 395–396. On accept, splice pointers and **swap AFS nodes** between routes (`xxafs=max(xx)`, `vvafs=max(yy)`, exchange pre/su/routeID — lines 467–482).
  - **Both AFS-free** (lines 483–559): same `uu`/`a`/`b` rebuild via `get_pd_pt`; node_location retag: `yy→nodeU_loc`, `xx→node_location(nodeV)`, `vv` flipped only if it had AFS.
  - else → reject.

**C++ note for m7/m8/m9:** the cleanest faithful port is: for each candidate, construct the two (or one) resulting routes as ordered node-ID lists exactly as the commit section would produce, run `get_pd_pt` on each changed route to get `(d1,d2,t1,t2)`, rebuild that route's `Route_related` row, recompute `pc_now` over all AFS routes with `AFSdelay_new`, form the penalty totals, and apply the SAME delta+threshold. The many MATLAB sub-branches exist only to do this incrementally; full re-evaluation collapses them, provided you keep the AFS-consistency REJECTION guards (an operator must reject moves that would put two AFS in one leg / leave an AFS-only fragment unless it explicitly deletes it) and the AFS add/remove distance credits.

---
### 4.5 `get_pd_pt` (route split re-evaluation) — full logic
```
function [d1,d2,t1,t2] = get_pd_pt(vrp, route_now, everTime, speed)   // route_now: node-ID list, NO depot
  if empty(route_now): return 0,0,0,0
  [a,n] = max(route_now)                       // a=max ID, n=its 1-based position
  if a > nb_customer:                          // route contains an AFS at position n
     d1 = sum dist route_now(1..n) inter-edges  + dist(depot, route_now(1))
        // (accumulate D over i=1..n-1, capture as d1 at i==n-1, reset)
     d2 = sum dist route_now(n..end) inter-edges + dist(route_now(end), depot)
     t1 = d1/speed + (n-1)*everTime
     t2 = d2/speed + (numel(route_now) - n + 1)*everTime
  else:                                         // no AFS
     d1 = sum inter-edges + dist(depot,first) + dist(last,depot)
     t1 = d1/speed + numel(route_now)*everTime
     d2 = 0; t2 = 0
```
Note the service-time counts: with AFS, pre-leg counts `n-1` customers (positions before AFS), post-leg counts `numel-n+1` (includes AFS position offset). Distances always use `dAll(id+1, id2+1)`, depot = index 1.

### 4.6 `AFSdelay_new` / `AFSdelay_recursion` / `get_pc_now` — see §3.2. (Re-implement verbatim for PC.)

### 4.7 `deleteAFS_delta(delete_idx, node_location, pd_v,pt_v,time_v,time_afs,distance_pre,distance_su, pm_now,pm, wm)`:
Removes route `delete_idx`'s rows from all Route_related vectors (`pd_v(delete_idx,:)=[]`, etc.) and `pm_now = max(pm_now - wm, 0)`. Called when the route holds only its AFS after the move.

### 4.8 `deleteAFS_node(routeID, delete_idx, node_location, predecessors, successor)`:
`node_deleteAfs = delete_idx + sum(node_location~=100)` = the AFS node-ID of that route (customers count = nodes not tagged 100). Then: decrement all `routeID > delete_idx` by 1; decrement all `predecessors/successor > node_deleteAfs` by 1; remove row `node_deleteAfs` from node_location/predecessors/successor/routeID. This keeps AFS node IDs contiguous after deletion.

============================================================
## 5. CONTROL FLOW
============================================================

### 5.1 ELS_mian (local-search driver that calls m1..m9)
Inputs: current solution's `chromR`, `node_location`, `Penalty_all`, `Route_related`, plus `isrepair`. Builds `Node_related` via `phrase_chromR`. Then:
```
whenLastTestedRI = zeros(nbClients,1); whenLastModified = zeros(nbRoutes,1);
nbMoves=0; loopID=0; searchCompleted=false;
while ~searchCompleted:
    if loopID>0: searchCompleted=true      // guarantees ≥2 passes
    rng(SEED+tspid)                          // reseed each pass
    for i=1..nbClients: with prob (randi mod nbGranular==0) shuffle correlatedVertices(i,:)
    if tspid~=1 and (toc - prev_time) > prev_time*20/nSol : break   // time cut for non-first
    for ii=1..nbClients:
        nodeU=ii; correlatedU=correlatedVertices(nodeU,:);
        lastTestRINodeU=whenLastTestedRI(nodeU); whenLastTestedRI(nodeU)=nbMoves;
        for jj over correlatedU:
            nodeV=correlatedU(jj);
            // (a "patch" block re-normalizes Route_related/Node_related when a route was emptied — lines 60–144; it deletes zero-distance routes, re-sorts AFS rows, recomputes PT/PC if needed. Replicate to keep indexing consistent after deletions.)
            if loopID==0 OR max(whenLastModified(routeID(nodeU)), whenLastModified(routeID(nodeV))) > lastTestRINodeU:
                try m1; if success continue
                try m2; … try m9  (in this exact order; first accepting move wins → "first-improvement")
                if predecessors(nodeV)==0:   // V is first in its route (depot-adjacent)
                    try Depot_m1, Depot_m2, Depot_m3, Depot_m8, Depot_m9
        if loopID~=1 and nodeV==correlatedU(end):
            try NewRoute_m1, NewRoute_m2, NewRoute_m3   // open a new route with U
    loopID++
// finalize: rebuild chromR, compute distance_Total, cost_Total, IsFeasible; write sol_table row.
```
- **Acceptance = first-improvement:** operators are tried in fixed order m1…m9 (then Depot_*, then NewRoute_* once per U); the first with `isSuccess` commits and moves to next V (`continue`).
- **Move ordering / neighborhoods** for m6–m9 are exactly (U over customers) × (V over U's granular neighbors), plus the operator-internal segment definitions above.
- **Termination:** at least two full passes; stop when a whole pass makes no accepted move (`searchCompleted` stays true).

### 5.2 Main_METS — memetic (HGS-style) driver
Parameters (exact): `split_prob=0.5`, `PT=527, PC=195, PD=430`, `penaltyScaleFactor=1.2`, `penaltyDecreaseFactor=0.85`, `popSizeMu=154`, `popSizeLambda=68`, `targetFeasible=0.2`, `nbLast=20`, `maxIterNonProd=300`, `maxIter=2000`, `timeLimit=100000`, `el=0.5` → `eliteNum=floor(0.5*154)=77`, `nc=0.2` → `nClosest=floor(0.2*154)=30`, `nbGranular=20`. `Penalty_all = [PT,PC,PD,0; 0,0,0,0]` (wm=0).

**Initialization (population of `popSizeMu*4 = 616` individuals, capped by maxIter/timeLimit):**
- `tsp_all` = `popSizeMu*4` random permutations of `1..nbClients` (seeded `rng(SEED+1)`).
- For each i: split the permutation into routes — with prob `split_prob` use `split_Dmax` else `split_Tmax` (§5.3); build details via `chromR_detail_all`; run `ELS_mian`; add to population via `PopManagement`; update best; add to `Last100Sol`.
- **Repair:** if the result is infeasible and `rand-0.5>1e-6` (≈50% prob), call `Repair_sol` = multiply every currently-positive penalty weight by `WP=10` and re-run `ELS_mian` once; if it becomes feasible, add to pop + update best.

**Main loop (`tspid = popSizeMu*4+1 … maxIter`):**
```
while nbIterNonProd <= maxIterNonProd and toc <= timeLimit:
    p1 = selectparents(...); p2 = selectparents(...)   // binary tournament (§5.4)
    offspring_tsp = Crossover(p1,p2,vrp)               // OX-style (§5.5)
    rng(SEED+tspid)
    split (Dmax or Tmax, prob 0.5)  [NOTE: MATLAB passes stale `tsp` here, not offspring_tsp — a bug, but reproduce as-is for fidelity: chromR is split from `tsp` while details use `offspring_tsp`]
    chromR_detail_all(offspring_tsp) ; ELS_mian ; PopManagement ; update best ; Last100Sol
    repair if infeasible & coin
    nbIterNonProd = 1 if new best else +1
    // penalty adaptation every nbLast iters once tspid>=100 (§5.6)
```

### 5.3 Split (initial route construction)
- **`split_Tmax`:** greedily append customers to a route while the route stays within `T_max_V` (checked by `T_able`); when the next customer would exceed Tmax, close the route and start a new one. Produces `chromR{i}=[1, custs+1, 1]`.
- **`split_Dmax`:** greedily append while pre-AFS distance `≤ V_Dmax`; when depot→…→customer plus depot→AFS exceeds Dmax, insert an AFS (node `nbClients+j+1`) splitting the route, and continue tracking the post-AFS leg against `V_Dmax`; close route when post-leg exceeds Dmax. Single-customer routes get `[1 a i+nb_customer+1 1]` (customer + its AFS).

### 5.4 Parent selection (`selectparents`): binary tournament — pick two random indices over feasible∪infeasible pool, keep the one with smaller `Fitness`; do this twice (p1, p2).

### 5.5 Crossover (`Crossover`): order-based (OX). Take parent gene sequences (customers only, AFS stripped), pick two cut points `point1<point2`; child inherits `[point1:point2]` from parent1, fills the rest from parent2 in cyclic order skipping already-present customers. Returns offspring customer permutation.

### 5.6 Survivor selection & biased fitness (`PopManagement`/`add2Pop`):
- Two subpopulations: feasible / infeasible. New solution inserted sorted by `cost_Total`.
- **Biased fitness** `Fitness = fitRank + (1 - eliteNum/popSize) * divRank`, where `fitRank` = normalized rank by cost (0..1), `divRank` = normalized rank by diversity. **Diversity** = `avgBrokenDist` = negative mean of the `nClosest` smallest **broken-pair distances** to other individuals; broken-pair distance between two solutions = (# positions where successor/predecessor pairs differ) / `vrp.last_customer` (AFS all mapped to `nbClients+1` before comparison).
- **Survivor removal:** when a subpop exceeds `popSizeMu+popSizeLambda`, repeatedly delete the individual with the largest (worst) biased fitness, preferring clones (broken-pair distance ≈ 0) first, down to `popSizeMu`.

### 5.7 Penalty adaptation (targetFeasible logic, Main_METS lines 171–239):
Every `nbLast=20` iterations, once `tspid>=100`, over the last `nbLast` recorded solutions (`Last100Sol`):
```
for each of T, C, D:
  fractionFeasible = (#solutions with that penalty == 0) / count
  if fractionFeasible <= targetFeasible - 0.05:   weight *= penaltyScaleFactor,   capped at 100000  (too few feasible → harsher)
  elif fractionFeasible >= targetFeasible + 0.05:  weight *= penaltyDecreaseFactor, floored at 0.1    (too many feasible → softer)
```
After changing a weight, rescale every infeasible individual's stored penalty by `(oldPenalty/oldWeight)*newWeight`, recompute their `cost_Total = penalty_D+penalty_C+penalty_T+distance_Total`, and re-sort/re-rank the infeasible pop (`infeasiblePop_updateBiasedFitnesses`). `targetFeasible=0.2`, band ±0.05.

### 5.8 Best-solution tracking / termination / output:
- `bestSolOverall` starts `IsFeasible=0, cost_Total=999999`; updated whenever a strictly-better (and feasible, per uti_updateBestSol) solution appears; `nbIterNonProd` resets on new best.
- Terminate when `nbIterNonProd > maxIterNonProd(300)` or `toc > timeLimit`.
- Final result: if `bestSolOverall.IsFeasible==0` → `Result=99999`; else recompute total distance from `bestSolOverall.chromR_move` using the **rounded** (2-decimal, floor) distance table and return it plus `Resulttime`.

============================================================
## 6. IMPLEMENTATION NOTES / EDGE CASES (do not lose)
============================================================
- **Acceptance threshold is `delta <= -1e-6`** (reject if `> -0.000001`). Use exactly this to reproduce tie behavior.
- **node_location tags {-1,1,100}** must be maintained on every commit; they drive which leg (pre/post) distances/times update and which routes may accept a moved AFS.
- **AFS per leg invariant:** an operator rejects any move that would place two AFS in one leg or leave an AFS-stranded fragment, UNLESS it explicitly deletes the emptied route (`isdelete`) with the correct distance credit (`-2*dAll(1,nbClients+2)` for full AFS-route removal, `-dAll(1,nbClients+2)` for a spared depot↔AFS edge).
- **AFS distance uses `dAll(1, nbClients+2)`** = depot↔canonical-AFS one-way distance; all AFS share this depot distance in the delta credits.
- **Service time `everTime`** is added per customer only (not depot, not AFS). Correct the service-time count whenever nodes move between legs/routes (the `±numel(...)*everTime` terms).
- **PC is global:** `AFSdelay_new` runs over ALL AFS-owning routes (the shared station resource), so any move changing any route's `time_afs`/`time_v` requires recomputing PC over the whole solution.
- **Distances are symmetric and unrounded during search;** rounding happens only in the final output.
- **1-indexing everywhere:** node ID `n` → distance index `n+1`; depot = node 0 = distance index 1. Port to 0-based carefully and consistently.

---

Source files (all absolute):
- Operators: `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/m6.m`, `.../m7.m`, `.../m8.m`, `.../m9.m` (and `.../m1.m` for the AFS-insertion pattern).
- Helpers: `.../Efficient local search/{AFSdelay_new.m, AFSdelay_recursion.m, get_pc_now.m, get_pd_pt.m, deleteAFS_delta.m, deleteAFS_node.m, get_chromR.m, phrase_chromR.m, Repair_sol.m, ELS_mian.m}`.
- Build/objective: `.../Load/chromR_detail_all.m`, `.../Load/get_vrp.m`.
- Driver/pop: `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Main_METS.m`, `.../Population management/{PopManagement.m, add2Pop.m, infeasiblePop_updateBiasedFitnesses.m}`, `.../Crossover selectparents/{selectparents.m, Crossover.m}`, `.../Two novel segmentation/{split_Dmax.m, split_Tmax.m}`.