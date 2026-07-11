I have now read all 12 files in full. Let me produce the exhaustive specification. Note the `get_pc_now` helper is called but not provided in the file set; I will document it precisely from context and flag it as a dependency to be sourced separately.

I have all the information needed. Here is the complete specification.

---

# GrVRP-PCAFS — METS Local Search: Exhaustive C++ Port Specification

This spec covers 12 MATLAB files: the 8 "depot-relative" / "new-route" move operators (`Depot_m1/2/3/8/9`, `NewRoute_m1/2/3`), the AFS-capacity penalty engine (`AFSdelay_new`, `AFSdelay_recursion`), and the route-deletion helpers (`deleteAFS_delta`, `deleteAFS_node`). All 8 operators share one common data layout and one common tail; understand Section 1–3 fully before the per-operator sections.

IMPORTANT PORTING NOTE FROM THE USER: In C++ you will FULLY RE-EVALUATE each candidate solution from the linked-list structure rather than reproduce the O(1) incremental `distance_pre`/`time_afs`/... bookkeeping. Therefore the many lines that patch `distance_pre(routeU) = distance_pre(routeU) + costOne` etc. are the MATLAB delta trick; you may treat them as "recompute route U and route V's distance/time split from scratch after the move." What you MUST preserve exactly: (a) which move is generated, (b) the exact resulting linked-list topology, (c) the feasibility guards that reject a move (`isSuccess=false; return`), (d) the acceptance rule `delta > -0.000001 → reject`, (e) the AFS penalty semantics, (f) the route-deletion trigger conditions, (g) `whenLastModified` bookkeeping. Below, each operator gives the MOVE + RESULTING TOPOLOGY + COST/FEASIBILITY semantics; the incremental arithmetic is documented too so you can cross-check a full re-evaluation.

---

## 1. DATA STRUCTURES

### 1.1 Node indexing (MATLAB 1-indexed; keep or offset carefully)
- Node 1 = depot.
- Customer `c` (1..nbClients) occupies node index `c+1` in distance table `dAll`. BUT in the linked-list arrays (`predecessors`, `successor`, `routeID`, `node_location`), a customer `c` is stored at array index `c` directly (the depot node 1 is NOT stored as a linked-list entry — the depot is represented implicitly by predecessor/successor value `0`).
- Careful: there are TWO indexing conventions in play:
  - `dAll(a+1, b+1)` = distance between node `a` and node `b`, where `a,b` are "logical node ids" with depot = 0, customer = its number, AFS = a number `> nbClients`. So `dAll` is 1-indexed and every logical id is offset by +1 when indexing `dAll`. `dAll(1, k+1)` = distance depot→node k.
  - The linked-list arrays (`predecessors` etc.) are indexed directly by the logical id (customer `c` at row `c`; an AFS instance at a row `> nbClients`). Value `0` in `predecessors`/`successor` means "the depot" (route endpoint).
- `dAll(1, nbClients+2)`: this is the distance depot → (logical node `nbClients+1`), i.e. depot→AFS-station-representative. It is used as the "cost of one AFS out-and-back stub" reference: a route consisting solely of an AFS contributes `2*dAll(1,nbClients+2)` of distance (depot→AFS→depot). Used when a route collapses to only-AFS and gets deleted.

### 1.2 `node_location` (per node, one value)
Encodes where a node sits relative to the route's single AFS refuel point:
- `-1` : node is in the **pre-AFS segment** (between depot start and the AFS).
- `+1` : node is in the **post-AFS segment** (between the AFS and depot end).
- `100`: the node **is the AFS** (refuel station node) of its route.
Each route has exactly one AFS node (value 100). A route reads depot → [nodes with loc −1] → AFS(100) → [nodes with loc +1] → depot. `max(node_location(routeID==r)) == 100` tests "route r contains an AFS" (100 is larger than ±1). Some operators temporarily set loc to `-100` as a sentinel then map `-100 → 100`.

### 1.3 Linked-list arrays (column vectors, one row per node)
`Node_related = [predecessors, successor, routeID, node_location]` (N×4), where N = total number of stored nodes = nbClients + (number of AFS instances currently in the solution). Columns:
- `predecessors(k)`: logical id of node before k in its route; `0` = depot (k is first in route).
- `successor(k)`: logical id of node after k; `0` = depot (k is last).
- `routeID(k)`: which route (vehicle) k belongs to, 1..R.
- `node_location(k)`: as in 1.2. (Internally handled as a row vector `node_location = Node_related(:,4)'` then transposed back on output.)

Route traversal: start from any node, follow `predecessors` back until it hits 0 (that's route start = first pre-AFS node, unless AFS is first), follow `successor` forward until 0 (route end). AFS instances live at rows `> nbClients`.

### 1.4 `Route_related` (R×7, one row per route/vehicle)
`Route_related = [pd_v(:,1), pd_v(:,2), pt_v, time_afs, time_v, distance_pre, distance_su]`. Columns, per route r:
1. `pd_v(r,1)` = distance-penalty of the **pre-AFS** half = `max(distance_pre(r) − V_Dmax, 0)`.
2. `pd_v(r,2)` = distance-penalty of the **post-AFS** half = `max(distance_su(r) − V_Dmax, 0)`.
3. `pt_v(r)`   = time-penalty of route = `max(time_afs(r) + time_su(r) − T_max_V, 0)`.
4. `time_afs(r)` = elapsed time from depot to (and including refuel at) the AFS = travel time of pre-AFS distance / speed + service times of pre-AFS customers + (AFS handling). Precisely: the "time up to end of AFS visit."
5. `time_v(r)` = total route time = `time_afs(r) + time_su(r)`.
6. `distance_pre(r)` = travel distance of pre-AFS half (depot → … → AFS).
7. `distance_su(r)`  = travel distance of post-AFS half (AFS → … → depot).
Derived inside every operator: `time_su = time_v − time_afs` (post-AFS half time). `V_Dmax` is the per-half max distance (fuel range between refuels); note the vehicle refuels at the AFS, so BOTH halves are separately range-limited.

### 1.5 `Penalty_all` (2×4)
Row 1 = **weights** (adaptive multipliers): `[wt, wc, wd, wm]` = `Penalty_all(1,1..4)`.
- `wt` = time-penalty weight, `wc` = AFS-capacity-penalty weight, `wd` = distance-penalty weight, `wm` = vehicle-count-penalty weight.

Row 2 = **current penalty totals**: `Penalty_all(2,1..4) = [pt, pc, pd, pm]`.
- `pt` = total time penalty (already weighted: `wt * sum(pt_v)`).
- `pc` = total AFS-capacity penalty (already weighted: `wc * rawPc`).
- `pd` = total distance penalty (already weighted: `wd * sum(sum(pd_v))`).
- `pm` = total vehicle-count penalty (already weighted).

Inside operators, `pt,pc,pd,pm` (lowercase) hold the OLD totals; `pt_now,pc_now,pd_now,pm_now` hold the recomputed NEW totals. `delta` compares them.

### 1.6 `afs_time_delay` (vector, length = #AFS/#routes-with-AFS)
Per-route "delayed arrival time at its AFS," output of the capacity-penalty solver. Passed through operators as state, recomputed by `AFSdelay_new`. Initialized to `time_afs` inside `AFSdelay_new`.

### 1.7 `vrp` struct fields used
- `vrp.nb_customer` (`nbClients`): number of customers.
- `vrp.T_Customer` (`everTime`): fixed service time per customer visit.
- `vrp.V_speed` (`speed`): travel speed (distance/time).
- `vrp.V_Dmax`: max distance per refuel half.
- `vrp.T_max_V` (`T_max`): max route duration.
- `vrp.T_Afs`: refuel handling duration at an AFS (also the "occupancy window" width for capacity conflicts).
- `vrp.C_Afs`: AFS capacity = number of vehicles that may refuel simultaneously in one `T_Afs` window.
- `vrp.V_nb`: number of available vehicles (baseline; extra routes beyond this incur `pm`).
- `vrp.distance_table` (`dAll`): full distance matrix, 1-indexed with +1 offset (see 1.1).

### 1.8 Bookkeeping scalars / flags (operator I/O)
- `nodeU`, `nodeV`: the two anchor nodes for the move (logical ids). `nodeX = successor(nodeU)`, `preU = predecessors(nodeU)`, `suX = successor(nodeX)`, `nodeY = successor(nodeV)`, `preV = predecessors(nodeV)`.
- `nodeU_loc/nodeX_loc/nodeV_loc/nodeY_loc`: their node_location.
- `routeU = routeID(nodeU)`, `routeV = routeID(nodeV)`.
- `nbMoves`: running accepted-move counter; `+1` on success.
- `searchCompleted`: bool; set `false` on any success (signals LS not converged).
- `isdelete`: flag; `0`/`1` meaning "route collapsed to only-AFS and must be deleted"; `-1` used as "not applicable" init in m2/m3 variants.
- `yu`, `yv`: output route ids touched (for the caller's dirty-route tracking). Set to `routeU`, `routeV` (or `routeAFS` in NewRoute\_\*).
- `whenLastModified`: vector, length = #routes; stamps `nbMoves` on touched routes; an entry is DELETED (`[]`) when its route is removed.
- `tspid`, `par_hgs`: passed through, unused in these operators.
- `isSuccess`: bool return; `true` if move accepted+applied, `false` if rejected (nothing changed — caller must discard all other outputs).

---

## 2. OBJECTIVE, PENALTIES, ACCEPTANCE RULE

### 2.1 Penalty component formulas (recompute after a move)
Given the post-move route split arrays:
- **Distance penalty (PD)**: for each route r, `pd_v(r,1)=max(distance_pre(r)−V_Dmax,0)`, `pd_v(r,2)=max(distance_su(r)−V_Dmax,0)`. Total `pd_now = wd * sum(sum(pd_v))` (sum over all routes and both halves).
- **Time penalty (PT)**: for each route r, `pt_v(r)=max(time_afs(r)+time_su(r)−T_max_V,0)` = `max(time_v(r)−T_max_V,0)`. Total `pt_now = wt * sum(pt_v)`.
- **AFS-capacity penalty (PC)**: `[rawPc, afs_time_delay] = AFSdelay_new(time_afs, time_v, T_max, T_Afs, C_Afs, afs_time_delay)`; then `pc_now = wc * rawPc`. (See Section 5 for the full solver.)
- **Vehicle-count penalty (PM)**: `pm_now = max(0, (R − V_nb) * wm)` where R = current number of routes (in `Depot_*` operators the number of routes does not change, so `pm_now = pm`; in `NewRoute_*` a route is added so `pm_now = max(0,(routeAFS − V_nb)*wm)` with `routeAFS = max(routeID)+1`).

### 2.2 Move delta and acceptance rule (identical in all 8 operators)
```
delta = costOne + costTwo
      + (pm_now - pm) + (pt_now - pt) + (pc_now - pc) + (pd_now - pd)
      [ - 2*dAll(1, nbClients+2)   IF isdelete==1 ]     % route removal refunds an AFS out-back stub
ACCEPT iff delta <= -0.000001  (i.e. strictly improving beyond tolerance 1e-6)
If delta > -0.000001 → isSuccess=false; return (reject, change nothing).
```
`costOne`/`costTwo` are the raw travel-distance deltas of the move (removal cost + insertion cost, per operator). The `-2*dAll(1,nbClients+2)` term is applied ONLY when the move empties a route down to only its AFS, which is then deleted — refunding that route's depot↔AFS out-and-back distance.

For a full re-evaluation port: compute `objective(solution) = totalTravelDistance + pt + pc + pd + pm` (all with current weights), and accept the candidate iff `obj(candidate) − obj(current) <= -1e-6`. This is equivalent to the delta test above (the AFS-stub term is already included in the candidate's total distance).

### 2.3 Weight adaptation (targetFeasible logic)
The weights `wt,wc,wd,wm` are ADAPTIVE and set by the caller/driver (not in these 12 files) based on the fraction of recent solutions that were feasible vs a `targetFeasible` ratio: if too few feasible, raise the relevant weight(s); if too many, lower them. These 12 files only READ the weights from `Penalty_all(1,:)`. (Standard HGS adaptation; implement in the driver, out of scope of these files but must exist.)

---

## 3. COMMON OPERATOR TAIL (executed on every accepted move)
After a move is accepted and the linked-list mutated, every operator does:
```
Route_related = [pd_v(:,1), pd_v(:,2), pt_v, time_afs, time_v, distance_pre, distance_su]
Node_related  = [predecessors, successor, routeID, node_location(as column)]
Penalty_all(2,:) = [pt_now, pc_now, pd_now, pm_now]   % (2,1)=pt,(2,2)=pc,(2,3)=pd,(2,4)=pm
nbMoves += 1
searchCompleted = false
% whenLastModified:
if isdelete==1:
    whenLastModified(routeU) = []      % remove the deleted route's stamp
    whenLastModified(routeV) = nbMoves
else:
    whenLastModified(routeU) = nbMoves
    whenLastModified(routeV) = nbMoves   % (NewRoute_* use routeAFS instead of routeV)
isSuccess = true
```
Note the order in the delete case: `whenLastModified(routeU)=[]` shifts indices; MATLAB deletes element `routeU` from the vector, so subsequent route ids shift down by one — consistent with `deleteAFS_node` decrementing `routeID > delete_idx`. In C++: erase the entry at index routeU and keep whenLastModified aligned with the (now renumbered) routes.

---

## 4. THE 8 MOVE OPERATORS

Every operator's signature is:
```
[isSuccess, nbMoves, searchCompleted, (isdelete,) yu, yv, whenLastModified,
 Penalty_all, Route_related, Node_related, afs_time_delay] =
   OP(nodeU, nodeV, vrp, nbMoves, searchCompleted, tspid, par_hgs,
      whenLastModified, Penalty_all, Route_related, Node_related, afs_time_delay)
```
(`Depot_m8` omits `isdelete` from its output list; all others include it.)

All operators begin by unpacking `vrp`, `Penalty_all`, `Route_related`, `Node_related` into local scalars/vectors exactly as listed in Section 1, and compute `time_su = time_v − time_afs`, and the node neighborhood (`preU,nodeX,suX,nodeV_loc,nodeY,...`).

The moves are **relocate / swap of segments across the AFS structure**, always inserting the moved segment **immediately after the depot start of route V** ("after v(routeV(1))" / "after NewRouteV(1)"). `nodeV` is assumed to be the FIRST node of route V (the routes are enumerated so V is a route-head anchor). Insertion places U (and/or X) as the new head of route V's relevant half, right after the depot.

---

### 4.1 `Depot_m1` — relocate single node U to head of route V

MOVE: remove `nodeU` from its position in route U (bridging `preU → nodeX`); insert `nodeU` as the first node after the depot in route V, i.e. new order `depot → nodeU → nodeV → …`. `nodeU` adopts `node_location = nodeV_loc` (it joins V's same half as V's head).

Guards / early rejects:
- `nodeU == nodeV` → reject.
- Case A `routeU == routeV` (intra-route relocate): always proceed to cost calc (no AFS-presence test).
- Case B `routeU != routeV` (inter-route):
  - Determine `isdelete`: `isdelete=1` iff route U has exactly 2 nodes AND one of them is the AFS (`sum(routeID==routeU)==2 && any(node_location(routeID==routeU)==100)`) — i.e. after U leaves, only the AFS remains → delete route U.
  - REQUIRE route V contains an AFS: `max(node_location(routeID==routeV))==100`. If not (`~=100`) → reject.

Cost:
- `costOne = −dAll(preU+1,nodeU+1) − dAll(nodeU+1,nodeX+1) + dAll(preU+1,nodeX+1)` (bridge out U from route U).
- `costTwo = −dAll(1,nodeV+1) + dAll(1,nodeU+1) + dAll(nodeU+1,nodeV+1)` (splice U between depot and V).

Distance/time reassignment (delta trick — for full re-eval, just recompute both routes' halves):
- Removal side (route U half `nodeU_loc`): if `nodeU_loc==−1` adjust `distance_pre(routeU)+=costOne`, `time_afs(routeU)+=costOne/speed − everTime`; if `==1` adjust `distance_su`/`time_su` with `−everTime` (U's service time leaves).
- Insertion side (route's half `nodeV_loc`): if `nodeV_loc==−1` adjust pre/`time_afs += costTwo/speed + everTime`; if `==1` adjust su/`time_su += costTwo/speed + everTime` (U's service time added).
- Intra-route (Case A): both adjustments hit `routeU`. Inter-route (Case B): removal hits routeU, insertion hits routeV.

Then recompute `pd_v`, `pt_v` for the touched route(s), `pd_now/pt_now`, `time_v=time_afs+time_su`, `[pc_now,afs_time_delay]=AFSdelay_new(...)`, `pc_now*=wc`, `pm_now=pm`.

If `isdelete==1`: call `deleteAFS_delta(routeU, …)` to drop route U's row and reduce `pm_now` by `wm` (floored at 0); then recompute pc/pt/pd from the shrunk arrays. `delta` then includes `−2*dAll(1,nbClients+2)`.

Accept test as Section 2.2. On accept, mutate links:
- `predecessors(nodeU)=0; successor(nodeU)=nodeV; predecessors(nodeV)=nodeU; node_location(nodeU)=nodeV_loc;`
- `if preU: successor(preU)=nodeX; if nodeX: predecessors(nodeX)=preU;`
- Case B also: `routeID(nodeU)=routeV;` and if `isdelete==1` call `deleteAFS_node(routeID, routeU, …)`.

Tail per Section 3. `yu=routeU, yv=routeV`.

---

### 4.2 `Depot_m2` — relocate the PAIR (U, X) to head of route V, keeping order U→X

`isdelete` init `-1`. `nodeYY=nodeV`.
MOVE: remove consecutive pair `nodeU, nodeX` (X = successor of U) from route U (bridge `preU → suX`); insert as `depot → nodeU → nodeX → nodeV → …` at head of route V. Both U and X take `node_location = nodeV_loc`.

Guards:
- Reject if `nodeX==0` (U has no successor) or `nodeX==nodeV` or `nodeU==nodeV`.
- Case A `routeU==routeV`: require `nodeX_loc ~= 100` (X must not be the AFS). If `nodeX_loc==100` → reject.
- Case B `routeU!=routeV`:
  - If `nodeX > nbClients` (X is an AFS, not a customer) → reject.
  - Else (`nodeX<=nbClients`): `isdelete=1` iff `sum(routeID==routeU)==3 && any(node_location(routeID==routeU)==100)` (route U had U, X, and the AFS → after pair leaves only AFS remains). Require route V has AFS (`max(node_location(routeID==routeV))==100`); else reject.

Cost:
- `costOne = −dAll(nodeU+1,nodeX+1) − dAll(preU+1,nodeU+1) − dAll(nodeX+1,suX+1) + dAll(preU+1,suX+1)`.
- `costTwo = +dAll(nodeU+1,nodeX+1) − dAll(1,nodeV+1) + dAll(1,nodeU+1) + dAll(nodeX+1,nodeV+1)`.
- Service-time deltas use `2*everTime` (two customers move). Same pattern as m1 but with factor 2.

`isdelete`/pc/pt/pd/delta and accept identical structure to m1 (delete refunds `−2*dAll(1,nbClients+2)`).

Link mutation on accept:
- `if preU: successor(preU)=suX; predecessors(nodeU)=0; if nodeX: successor(nodeX)=nodeV; if suX: predecessors(suX)=preU; predecessors(nodeV)=nodeX; node_location(nodeU)=nodeV_loc; node_location(nodeX)=nodeV_loc;` (Note: U becomes head, X between U and V; predecessors(nodeV)=nodeX, and U's predecessor is 0. successor(nodeU) is left implicitly = nodeX from before, which is preserved.)
- Case B: `routeID(nodeU)=routeV; routeID(nodeX)=routeV;` and if delete, `deleteAFS_node(...)`.

Tail per Section 3.

---

### 4.3 `Depot_m3` — relocate the pair but REVERSED (X, U) to head of route V

Identical to `Depot_m2` EXCEPT the inserted order is `depot → nodeX → nodeU → nodeV → …` (segment reversed).

Cost difference (only `costTwo` changes):
- `costOne` same as m2.
- `costTwo = +dAll(nodeU+1,nodeX+1) − dAll(1,nodeV+1) + dAll(1,nodeX+1) + dAll(nodeU+1,nodeV+1)` (depot→X, then U→V).

Link mutation on accept (reversed):
- `if preU: successor(preU)=suX; predecessors(nodeU)=nodeX; successor(nodeU)=nodeV; if nodeX: predecessors(nodeX)=0, successor(nodeX)=nodeU; if suX: predecessors(suX)=preU; predecessors(nodeV)=nodeU; node_location(nodeU)=nodeV_loc; node_location(nodeX)=nodeV_loc;`
- Case B: `routeID(nodeU)=routeV; routeID(nodeX)=routeV;` and delete handling.

All guards, isdelete logic, delta, accept, tail identical to m2.

---

### 4.4 `Depot_m8` — move the ENTIRE tail of route U (from U to end) to head of route V

`afs_swap=0; isdelete=0`. Requires `routeV != routeU` and `nodeX != 0` else reject.

Build the moved segment `xx`: starting at `nodeU`, walk successors collecting `xx = [nodeU, successor, successor², …]` up to route length, stop at 0. Then remove the FIRST element (`xx(1)=[]`) so `xx` = everything AFTER U (i.e. `nodeX` onward to end of route U). Also accumulate `d_xx` = sum of consecutive intra-`xx` distances. If `xx` empty → reject.

Two sub-cases by whether the moved tail contains the AFS (`max(node_location(xx))`):

**Sub-case 8a: `max(node_location(xx)) != 100`** (tail carries NO AFS — pure customers):
- MOVE: detach tail `xx` after U (U becomes new route-U end: `successor(nodeU)=0`), and prepend the tail (REVERSED) to route V's head: `depot → xx-reversed → nodeV → …`. `xx` nodes join `routeID=routeV`, `node_location(xx)=nodeV_loc`.
- Cost: `costOne = −d_xx − dAll(nodeU+1,nodeX+1) − dAll(xx(end)+1,1) + dAll(nodeU+1,1)` (route U loses its tail, U now goes straight to depot). `costTwo = +d_xx − dAll(1,nodeV+1) + dAll(xx(end)+1,1) + dAll(nodeX+1,nodeV+1)`.
- Distance/time: routeU half chosen by `nodeX_loc==−1` (pre) else `nodeU_loc==1` (su); service time removed = `numel(xx)*everTime`. routeV half by `nodeV_loc`; service added `numel(xx)*everTime`.
- PC solver called on the sub-vector `time_afs(1:numel(node_location)-nbClients)` / `time_v(...)` — i.e. only over routes that have AFS rows (the first `#routes` entries correspond to AFS-bearing routes; this slices to AFS-count length).
- Link mutation: `successor(nodeU)=0; predecessors(nodeV)=nodeX;` then rebuild `xx` as a reversed chain in front of V: set `predecessors(xx(end))=0`, for i=1..len-1 `predecessors(xx(i))=xx(i+1)`; `successor(xx(1))=nodeV`, for i=2..len `successor(xx(i))=xx(i-1)`. `routeID(xx)=routeV; node_location(xx)=nodeV_loc`.

**Sub-case 8b: `max(node_location(xx)) == 100`** (moved tail CONTAINS the AFS):
- Extra guard: require route V has NO AFS (`max(node_location(routeID==routeV)) ~= 100`); if V already has AFS → reject (can't have two AFS).
- This performs an AFS-carrying tail move plus a swap of route roles: route U keeps only its pre-AFS prefix (its post-AFS distance goes to 0), and route V absorbs the AFS-bearing tail. The code recomputes both routes' pre/su distances and times explicitly (this is a from-scratch recompute already — replicate its formulas or simply re-evaluate the two routes after applying the topology change):
  - `distance_pre(routeU) = d1u + d2u + costOne; distance_su(routeU)=0; time_su(routeU)=0; time_afs(routeU) = distance_pre(routeU)/speed + (nodesInRouteU − numel(xx))*everTime`.
  - `distance_pre(routeV)=d2u; distance_su(routeV)=d1v+d2v+costTwo − distance_pre(routeV); time_afs(routeV)=t2u − everTime; time_v(routeV)=(d1v+d2v+costTwo)/speed + numel(xx)*everTime + nodesInRouteV*everTime; time_su(routeV)=time_v(routeV)−time_afs(routeV)`.
- PC solver called on `time_afs`/`time_v` with route U's entry excluded (`[1:routeU-1, routeU+1:end]`) because U now has no AFS.
- Link mutation same reversed-prepend as 8a plus set `afs_swap=1` and: `node_location(xx) = −node_location(xx); node_location(routeID==routeV)=1; node_location(node_location==−100)=100;` — i.e. within the moved segment, pre/su flip sign, everything now in V is marked post-AFS (`+1`), and the AFS sentinel `−100` maps back to `100`.
- AFTER building `Node_related`/`Route_related`, because `afs_swap==1`, SWAP route ids U↔V globally (`routeID==routeV → routeU`, `routeID==routeU → routeV`, via a `-1` temp) and swap `Route_related` rows U and V. This makes route V (now carrying the AFS-tail) take U's slot ordering. Net effect: the AFS-bearing composite route ends up correctly numbered.

Both sub-cases: standard delta/accept (`isdelete` stays 0 → no `-2*dAll` term). Tail per Section 3, but `whenLastModified(routeU)=nbMoves; whenLastModified(routeV)=nbMoves;` (no deletion branch). `yu=routeU, yv=routeV`.

For a full re-evaluation port: implement 8 as "cut route U at U (keep depot..U), take the remainder, reverse it, and prepend after route V's depot; if the remainder contains the AFS, additionally require V had no AFS and let the composite route be re-split so U's post-half is empty." Then recompute both routes fully and test delta. The `afs_swap` row/id swap is purely a renumbering to keep the AFS route in a canonical slot; in C++ just place nodes in the correct route object.

---

### 4.5 `Depot_m9` — swap the tail of route U with the tail of route V (2-opt* style exchange)

`isdelete=0`. Require `routeV != routeU`.

Build `xx` = tail of route U after U (same construction as m8: collect successors from nodeU, drop first, drop zeros) with `d_xx`. Reject if (`numel(xx)==1 && node_location(xx)==100`) — can't move a lone AFS this way.
Build `yy` = route V's segment starting AT `nodeV` and walking successors to the end, WITH zeros dropped (note: `yy` KEEPS `nodeV` as its first element, unlike xx). Accumulate `d_yy`. `routeV_now = yy`. Reject if `numel(yy)==0` or `numel(xx)==0`.

MOVE: exchange the two tails between routes — route U's post-U tail (`xx`) goes to route V, route V's from-nodeV tail (`yy`) goes to route U, spliced so `depot(U)…nodeU → nodeV(head of yy) …` and `depot(V) … → xx`. Both segments are re-linked reversed as in m8; the AFS node of each segment is repositioned (see below).

Three sub-cases by AFS membership of the two tails:

**9a: both tails contain the AFS** (`max(node_location(xx))==100 && max(node_location(yy))==100`):
- If `numel(xx)==1 && node_location(xx)==100` set `isdelete=1` (route U's tail is just the AFS → after swap route V collapses to only-AFS and gets deleted; `delete_idx=routeV`).
- Costs:
  - `costOne = −d_xx − dAll(nodeU+1,nodeX+1) − dAll(xx(end)+1,1) + d_yy + dAll(nodeU+1,nodeV+1) + dAll(yy(end)+1,1)`.
  - `costTwo = −d_yy − dAll(1,nodeV+1) − dAll(yy(end)+1,1) + d_xx + dAll(1,xx(end)+1) + dAll(nodeX+1,1)`.
- Full explicit recompute of both routes' pre/su distances and afs/su times (replicate or re-evaluate). Key formulas: `distance_su(routeU)=d2v; distance_pre(routeV)=d2u; distance_su(routeV)=d_xx+dAll(xx(end)+1,1)−d2u+dAll(nodeX+1,1); distance_pre(routeU)=d1u+d2u+costOne−distance_su(routeU); time_su(routeU)=t2v; time_afs(routeV)=t2u−everTime; time_v(routeU)=(d1u+d2u+costOne)/speed + nodesU*everTime + (numel(yy)−numel(xx))*everTime; time_v(routeV)=(d1v+d2v+costTwo)/speed + nodesV*everTime − (numel(yy)−numel(xx))*everTime; time_afs(routeU)=time_v(routeU)−time_su(routeU); time_su(routeV)=time_v(routeV)−time_afs(routeV)`.
- PC over AFS-bearing prefix (`1:numel(node_location)−nbClients`). If delete, `deleteAFS_delta(routeV,…)` then recompute; delta gets `−2*dAll(1,nbClients+2)`.
- Link mutation: complex reconnection swapping the two tails and RE-POSITIONING each tail's AFS node. Precisely (on accept):
  - `successor(nodeU)=nodeV;` rebuild xx's internal reversed links (`predecessors(xx(end))=0; predecessors(xx(i))=xx(i+1); successor(nodeX)=0; successor(xx(i))=xx(i-1)`); `predecessors(nodeV)=nodeU; routeID(xx)=routeV; routeID(yy)=routeU; node_location(xx)=−node_location(xx); node_location(node_location==−100)=100;`
  - AFS reposition: `xxafs=max(xx); vvafs=max(yy)` (the AFS instances are the max-indexed rows in each segment). Restore `routeID(xxafs)=routeU; routeID(vvafs)=routeV`. Then splice each AFS node into the OTHER's neighborhood: capture pre/su of both, cross-link so xxafs takes vvafs's neighbors and vice versa (the block lines 182–189). If `isdelete==1`, `deleteAFS_node(routeV,…)`.

**9b: NEITHER tail contains the AFS** (`max(node_location(xx))~=100 && max(...yy)~=100`): pure customer-tail exchange.
- Same `costOne`/`costTwo` formulas.
- Distance/time via half chosen by `nodeU_loc` (route U) and `nodeV_loc` (route V); service adjustments `(numel(yy)−numel(xx))*everTime` on U and `(numel(xx)−numel(yy))*everTime` on V.
- Link mutation: `successor(nodeU)=nodeV;` reverse-link xx; `successor(nodeX)=0` (route V's old continuation cut); `predecessors(nodeV)=nodeU; routeID(xx)=routeV; routeID(yy)=routeU; node_location(yy)=nodeU_loc; node_location(xx)=−1; node_location(node_location==−100)=100;`. No delete.

**9c: exactly one tail has the AFS** → `else: isSuccess=false; return` (reject; mixing an AFS tail with a non-AFS tail is disallowed).

Standard delta/accept. Tail per Section 3 with the delete branch when `isdelete==1` (delete removes `routeV`'s stamp: note here the deleted route is routeV, but the tail still executes `whenLastModified(routeU)=nbMoves; whenLastModified(routeV)=nbMoves` in m9 — m9's tail does NOT have the special delete branch; it unconditionally stamps both. The row removal is handled inside `deleteAFS_delta`/`deleteAFS_node`. Match the source: m9 always does `whenLastModified(routeU)=nbMoves; whenLastModified(routeV)=nbMoves;` and relies on the delete helpers for row removal — but note this leaves a length mismatch that the caller must reconcile; replicate exactly).

`yu=routeU, yv=routeV`.

---

### 4.6 `NewRoute_m1` — relocate single node U into a BRAND-NEW route (with a fresh AFS)

`isdelete=0`. Creates a new route consisting of `depot → nodeU → AFS → depot`.
Key indices: `routeAFS = max(routeID)+1` (new route id), `node_addAfs = numel(routeID)+1` (new AFS node row appended at end). `yu=routeU; yv=routeAFS`. Reject if `routeU==routeAFS`.

`isdelete=1` iff route U had exactly 2 nodes and one is the AFS (`sum(routeID==routeU)==2 && any(node_location(routeID==routeU)==100)`): moving U out empties route U to only-AFS → delete it.

MOVE/topology: remove U from route U (`preU → nodeX` bridge); create new route `routeAFS` with U in its PRE-AFS half (`node_location(nodeU)=−1`) and a new AFS node `node_addAfs` (`node_location=100`) as the route tail. New route reads `depot → U → AFS → depot`.

Costs:
- `costOne = −dAll(preU+1,nodeU+1) − dAll(nodeU+1,nodeX+1) + dAll(preU+1,nodeX+1)` (remove U).
- `costTwo = +dAll(1,nodeU+1) + dAll(nodeU+1,node_addAfs+1) + dAll(node_addAfs+1,1)` (depot→U→AFS→depot).
- New route split: `distance_pre(routeAFS)=dAll(1,nodeU+1)+dAll(nodeU+1,node_addAfs)`, `distance_su(routeAFS)=dAll(node_addAfs,1)`, `time_afs(routeAFS)=distance_pre/speed+everTime`, `time_su(routeAFS)=distance_su/speed+everTime`.
- `pm_now = max(0,(routeAFS − V_nb)*wm)` (a route was ADDED).
- Removal-side route U half via `nodeU_loc` (`−everTime`).

If `isdelete==1`: `deleteAFS_delta(routeU,…)`, recompute, delta gets `−2*dAll(1,nbClients+2)`.

Link mutation on accept:
- `predecessors(node_addAfs)=nodeU; successor(node_addAfs)=0; if preU: successor(preU)=nodeX; predecessors(nodeU)=0; successor(nodeU)=node_addAfs; if nodeX: predecessors(nodeX)=preU; routeID(nodeU)=routeAFS; routeID(node_addAfs)=routeAFS; node_location(node_addAfs)=100; node_location(nodeU)=−1;`
- If delete: `deleteAFS_node(routeU,…)`.

Tail: `whenLastModified(routeU)=nbMoves; whenLastModified(routeAFS)=nbMoves;` (no delete branch in the tail of NewRoute_m1 — but note if isdelete, deleteAFS already shifted ids; replicate exactly: NewRoute_m1's tail unconditionally stamps routeU and routeAFS).

---

### 4.7 `NewRoute_m2` — relocate pair (U, X) into a new route, order U→X

`isdelete=-1`. New route `depot → U → X → AFS → depot`. `routeAFS=max(routeID)+1`, `node_addAfs=numel(routeID)+1`.
Guards: reject if `nodeX==0`, or `nodeX>nbClients` (X is an AFS). If `nodeX<=nbClients`: if `sum(routeID==routeU)==3 && any(node_location(routeID==routeU)==100)` set `isdelete=1` AND immediately `isSuccess=false; return` (this specific collapse is DISALLOWED for m2 — you may not empty a 3-node route by moving a customer pair into a new route). Else `isdelete=0`.

Costs:
- `costOne = −dAll(preU+1,nodeU+1) − dAll(nodeU+1,nodeX+1) − dAll(nodeX+1,suX+1) + dAll(preU+1,suX+1)`.
- `costTwo = +dAll(1,nodeU+1)+dAll(nodeU+1,nodeX+1)+dAll(nodeX+1,node_addAfs+1)+dAll(node_addAfs+1,1)`.
- New route: `distance_pre(routeAFS)=dAll(1,nodeU+1)+dAll(nodeU+1,nodeX+1)+dAll(nodeX+1,node_addAfs+1)`, `distance_su(routeAFS)=dAll(node_addAfs,1)`, `time_afs(routeAFS)=distance_pre/speed+2*everTime`, `time_su=distance_su/speed+everTime`. Removal service `−2*everTime`.
- `pm_now=max(0,(routeAFS−V_nb)*wm)`. Since isdelete is forced to abort, the `−2*dAll` branch is effectively unreachable.

Link mutation: `predecessors(node_addAfs)=nodeX; successor(node_addAfs)=0; if preU: successor(preU)=suX; predecessors(nodeU)=0; if nodeX: successor(nodeX)=node_addAfs; if suX: predecessors(suX)=preU; routeID(nodeU/nodeX/node_addAfs)=routeAFS; node_location(nodeU)=−1; node_location(nodeX)=−1; node_location(node_addAfs)=100;`. (U head, then X, then AFS.)

Tail: stamps routeU and routeAFS.

---

### 4.8 `NewRoute_m3` — relocate pair into new route, REVERSED order X→U

Identical to NewRoute_m2 except new route is `depot → X → U → AFS → depot`. Only `costTwo` and link order differ:
- `costTwo = +dAll(1,nodeX+1)+dAll(nodeX+1,nodeU+1)+dAll(nodeU+1,node_addAfs+1)+dAll(node_addAfs+1,1)`.
- New route split: `distance_pre(routeAFS)=dAll(1,nodeX+1)+dAll(nodeX+1,nodeU+1)+dAll(nodeU+1,node_addAfs+1)`, `distance_su=dAll(node_addAfs+1,1)`.
- Link mutation: `predecessors(node_addAfs)=nodeU; successor(node_addAfs)=0; if preU: successor(preU)=suX; predecessors(nodeU)=nodeX; successor(nodeU)=node_addAfs; if nodeX: predecessors(nodeX)=0, successor(nodeX)=nodeU; if suX: predecessors(suX)=preU;` route ids and locations as m2.
- Same `isdelete==1 → reject` guard, same `pm_now`, same tail.

---

## 5. AFS-CAPACITY PENALTY SOLVER (PC)

The AFS is a *private capacitated* station: at most `C_Afs` vehicles can be refueling within any `T_Afs`-wide time window. Each AFS-bearing route arrives at its AFS at time `time_afs(r)` and occupies the station for `T_Afs`. If more than `C_Afs` routes overlap, that's a capacity violation; vehicles may be DELAYED (shifting their refuel later) to deconflict, but a route only has `time_V_shifting(r)=max(T_max − time_v(r),0)` slack before it violates `T_max`. The solver tries to deconflict by delaying; residual unavoidable overlap becomes `pc_now`.

### 5.1 `AFSdelay_new(time_afs, time_v, T_max, T_Afs, C_Afs, afs_time_delay) → [pc_now, afs_time_delay]`
```
if numel(time_afs)==0: pc_now=0; return          % no AFS routes → no capacity penalty
time_V_shifting = max(T_max - time_v, 0)          % per-route slack before T_max violated
afs_time_delay  = time_afs                        % initialize delayed arrival = actual arrival
conflict_table  = zeros(n, n)                     % n = numel(time_afs)
[pc_now, afs_time_delay] = AFSdelay_recursion(time_afs, time_v, T_max, T_Afs, C_Afs,
                                              time_V_shifting, afs_time_delay, conflict_table)
```

### 5.2 `AFSdelay_recursion(...) → [pc_now, afs_time_delay]`
Purpose: iteratively resolve pairwise AFS occupancy overlaps by delaying whichever route can absorb the delay within its slack; when a pair is unresolvable one way, branch (backtrack) trying both delay directions and keep the lower `pc_now`.

```
pc_now = get_pc_now(afs_time_delay, T_max, time_v, time_afs, T_Afs, C_Afs)   % (see 5.4)
delay_over = 1
while delay_over:
    delay_over = 0
    for i = 1 .. n-1:
      for ii = i+1 .. n:
        % order the pair so t2 >= t1 (t2 = later/greater delayed arrival)
        if afs_time_delay(i) >= afs_time_delay(ii):
            t1=afs_time_delay(ii); t2=afs_time_delay(i); q1=time_V_shifting(ii); q2=time_V_shifting(i); t3=1
        else:
            t1=afs_time_delay(i);  t2=afs_time_delay(ii); q1=time_V_shifting(i);  q2=time_V_shifting(ii); t3=2
        if t1 + T_Afs <= t2:  continue        % no overlap → no conflict
        % both directions infeasible within slack → record conflict, leave as violation
        if t1+T_Afs-t2 > q2 AND t2+T_Afs-t1 > q1:
            conflict_table(i,ii)+=1; conflict_table(ii,i)+=1; continue
        % only earlier one can be resolved by delaying t1 to just after t2's window:
        if t1+T_Afs-t2 > q2:
            delaytime = t2 + T_Afs - t1;  t1 += delaytime
            (delay the earlier route: if t3==1 delay index ii else index i),
            reduce that route's time_V_shifting by delaytime; reset i=1,ii=1; continue
        % only later one resolvable by delaying t2:
        if t2+T_Afs-t1 > q1:
            delaytime = t1 + T_Afs - t2;  t2 += delaytime
            (delay the later route accordingly), reduce its slack; reset i=1,ii=1; continue
        % (this final guard is logically unreachable given the two above, but present:)
        if t1+T_Afs-t2 > q2 AND t2+T_Afs-t1 > q1:
            [afs_time_delay, time_V_shifting, pc_now] = backtrack(...)   % branch both ways
            delay_over = 1; break
      if delay_over: break
pc_now = get_pc_now(afs_time_delay, T_max, time_v, time_afs, T_Afs, C_Afs)
```
Notes: `t3` records which of (i,ii) is the later one so the correct index gets the delay. The `i=1;ii=1` reassignments inside the loop body attempt to restart scanning (in MATLAB these do NOT actually reset the `for` loop counters — the `for` header controls them — so effectively they are no-ops for loop control and the `continue` proceeds to the next `ii`; the OUTER `while delay_over` provides the real restart only when a backtrack fires). Replicate the observable behavior: keep scanning pairs; when a resolvable overlap is found, apply the delay and continue scanning; loop the whole pass again only when `backtrack` sets `delay_over`.

### 5.3 `backtrack(...)` (nested helper)
Tries BOTH deconfliction directions for an unresolvable-by-single-move pair and keeps whichever yields the lower `pc_now`. It:
1. saves current state; delays route "1" by `delaytime=t2+T_Afs-t1`, recurses `AFSdelay_recursion`; if result `pc_now < pc_original`, adopt it.
2. restores; delays route "2" by `delaytime=t1+T_Afs-t2`, recurses; if better, adopt.
Returns best `afs_time_delay`, `time_V_shifting`, `pc_now`.
CAVEAT (bug to replicate faithfully OR guard against): as written, `backtrack` recurses with an undefined `conflict_table` and `original_afs_time_delay` is referenced before assignment in some paths; the branch that fires it is guarded by a condition already handled above, so in practice `backtrack` is effectively dead code / rarely entered. Safe port: implement the two-way branch-and-keep-min semantics, initializing `conflict_table` to zeros(n,n) on each recursion and seeding `original_afs_time_delay = afs_time_delay` before use. The intended semantics: explore delaying either conflicting route, recompute pc, keep the minimum.

### 5.4 `get_pc_now(afs_time_delay, T_max, time_v, time_afs, T_Afs, C_Afs)` — DEPENDENCY NOT IN FILE SET
This helper is CALLED by `AFSdelay_recursion` but its source is NOT among the 12 files provided. It computes the raw (unweighted) AFS-capacity penalty given the (possibly delayed) arrival times. Its required semantics (infer/obtain from the METS codebase before implementing): for the set of AFS occupancy intervals `[afs_time_delay(r), afs_time_delay(r)+T_Afs)`, count the maximum simultaneous overlap; the penalty reflects the amount by which concurrent refuels exceed `C_Afs` (over-capacity), plus any `T_max` violation induced by the applied delays (`time_v` shifted by delay). ACTION FOR PORTER: locate `get_pc_now.m` in the METS repo and spec it separately — do NOT invent the formula. The interface is fixed: input the delayed arrival times + capacity params, output a nonnegative scalar `pc_now`; `AFSdelay_new` multiplies by `wc`.

---

## 6. ROUTE-DELETION HELPERS

Triggered when a move empties a route to just its AFS. Two parts: `deleteAFS_delta` updates the R-indexed route arrays and `pm`; `deleteAFS_node` updates the N-indexed node arrays. Both compute `node_deleteAfs = delete_idx + sum(node_location~=100)`: the row index of the AFS node to remove. Meaning: AFS nodes are stored AFTER all customers; `sum(node_location~=100)` = number of non-AFS (customer + in-transit) nodes = base offset; adding `delete_idx` (the route number) gives the specific AFS node's row. (Both assert this resolves to a single index.)

### 6.1 `deleteAFS_delta(delete_idx, node_location, pd_v, pt_v, time_v, time_afs, distance_pre, distance_su, pm_now, pm, wm)`
```
node_deleteAfs = delete_idx + sum(node_location != 100)   % (assert scalar)
pd_v(delete_idx,:) = []      % remove route row from pd_v (both halves)
pt_v(delete_idx)   = []
time_v(delete_idx) = []
time_afs(delete_idx) = []
distance_pre(delete_idx) = []
distance_su(delete_idx)  = []
pm_now = max(pm_now - wm, 0)   % one fewer vehicle → reduce vehicle penalty
return updated arrays + pm_now (pm unchanged)
```

### 6.2 `deleteAFS_node(routeID, delete_idx, node_location, predecessors, successor)`
```
node_deleteAfs = delete_idx + sum(node_location != 100)   % (assert scalar)
% renumber routes above the deleted one:
routeID(routeID > delete_idx) -= 1
% renumber node rows above the deleted AFS node:
predecessors(predecessors > node_deleteAfs) -= 1
successor(successor > node_deleteAfs)       -= 1
% remove the AFS node row from all node arrays:
node_location(node_deleteAfs) = []
predecessors(node_deleteAfs)  = []
successor(node_deleteAfs)     = []
routeID(node_deleteAfs)       = []
return updated arrays
```
C++ equivalent: erase route object `delete_idx` (all higher route ids shift down by 1) and erase AFS node row `node_deleteAfs` (all node references above it decrement). Keep `whenLastModified` aligned (the operator tail erases its entry).

---

## 7. CONTROL FLOW CONTEXT (how these operators are driven)

These 12 files are the LOCAL-SEARCH move library called from an outer HGS/memetic driver (not in this file set). The contract the driver must honor:
- Enumerate anchor pairs `(nodeU, nodeV)` where `nodeV` is a route-head, apply each applicable operator; the operator itself gates feasibility and improvement.
- Acceptance is **first-improvement with tolerance 1e-6** (each operator applies the move iff strictly improving); `searchCompleted` flips to `false` whenever ANY move succeeds — the driver loops the neighborhood until a full pass yields no success (then `searchCompleted` stays true → LS converged).
- `nbMoves` monotonically counts accepted moves; `whenLastModified[r]` records the last move touching route r, enabling the driver to skip re-examining unchanged route pairs (don't-look bits at route granularity).
- On rejection (`isSuccess=false`) the operator returns WITHOUT mutating any output; the driver must discard the returned matrices and keep the prior state.
- `NewRoute_*` grow the solution by one route (a fresh AFS route); `Depot_m8/m9` never delete; `Depot_m1/2/3` and `NewRoute_m1` may delete a route via the AFS helpers, adjusting `pm` and shifting all ids.
- Penalty weights (`Penalty_all(1,:)`) are adapted OUTSIDE these files by the driver based on feasibility vs `targetFeasible`; these operators only read them. Objective used for population ranking / biased fitness, repair probability, parent/survivor selection all live in the driver and are OUT OF SCOPE of these 12 files (they are not referenced here).

---

## 8. PORTER CHECKLIST / GOTCHAS
1. `dAll` uses +1 offset for every logical id; the AFS-station reference distance is `dAll(1, nbClients+2)`. A collapsed-to-AFS route's stub = `2*dAll(1,nbClients+2)`.
2. `node_location`: −1 pre-AFS, +1 post-AFS, 100 AFS, transient −100 sentinel → remap to 100. `max(loc(routeID==r))==100` ⇔ route r has an AFS.
3. Acceptance tolerance is exactly `-0.000001`; reject if `delta > -0.000001`.
4. Delete triggers: m1/NewRoute_m1 when route had 2 nodes incl. AFS; m2/m3 when route had 3 nodes incl. AFS (but NewRoute_m2/m3 ABORT in that case instead of deleting); m9 when U's tail is a lone AFS. m8 never deletes.
5. In m2/m3/8/9, `2*everTime` / `numel*everTime` service-time bookkeeping scales with #moved customers.
6. `deleteAFS_*` renumber both route ids (`>delete_idx`) and node rows (`>node_deleteAfs`); keep `whenLastModified` and all arrays consistent after deletion.
7. `get_pc_now` is an EXTERNAL dependency not in this file set — source it from the METS repo; do not fabricate its formula.
8. The `backtrack`/`i=1;ii=1` internals of `AFSdelay_recursion` are partly dead/buggy; port the OBSERVABLE semantics (pairwise delay-to-deconflict within slack; unresolvable overlaps become penalty), not the literal broken recursion.
9. For a full-re-evaluation port (the intended approach): ignore the incremental `distance_pre += cost` arithmetic; instead, after applying each operator's exact link topology, recompute each affected route's `distance_pre/su`, `time_afs/su`, `time_v`, `pd_v`, `pt_v`, run `AFSdelay_new` for PC, and compute `pm` from route count; accept iff total objective improves by ≥1e-6.

---

FILE PATHS (all absolute, for reference):
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/Depot_m1.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/Depot_m2.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/Depot_m3.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/Depot_m8.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/Depot_m9.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/NewRoute_m1.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/NewRoute_m2.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/NewRoute_m3.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/AFSdelay_recursion.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/AFSdelay_new.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/deleteAFS_delta.m
- /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Efficient local search/deleteAFS_node.m

DEPENDENCY MISSING FROM FILE SET (must be sourced separately before the port compiles): `get_pc_now.m` (called by `AFSdelay_recursion`), and the outer HGS driver (population construction, parent/survivor selection, biased fitness, repair probability, weight adaptation vs `targetFeasible`, neighborhood enumeration order) — none of which are in these 12 files.