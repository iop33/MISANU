I have everything needed. Here is the exhaustive, implementation-ready specification.

---

# SPECIFICATION: `split_Tmax` and `split_Dmax` (GrVRP-PCAFS "Two Novel Segmentation" routines)

This spec covers the two MATLAB Split routines that convert a giant TSP-style permutation of customers (`tsp`) into a set of vehicle routes (`chromR`), one splitting on the time bound `Tmax`, the other on the distance/range bound `Dmax`. Everything below is derived directly from the MATLAB source plus the surrounding model files that define the data conventions. A C++ programmer can implement from this without seeing the MATLAB.

**CRITICAL PORTING NOTE — indexing.** All MATLAB code is 1-indexed. In a C++ port you have two choices; pick one and be consistent. The spec below states BOTH the mathematical node id (used in `distance_table` lookups) and the 1-based array positions. Where the MATLAB writes `+1` on a node value, that is converting a 0-based "logical" node id into the 1-based `distance_table` row/column index — this is a real arithmetic operation you must reproduce, not a language artifact.

---

## 1. DATA STRUCTURES

### 1.1 Node numbering convention (the single most important thing to get right)

There are three logical categories of node. Two parallel numbering systems are used, and the code moves between them.

**"Logical" (0-based) ids — the values stored in `tsp`:**
- `tsp` contains **customer indices `0 .. nb_customer-1`** (0-based customer ids). This is why every place a `tsp` value is put into a route it gets `+1` added.

**"distance_table" (1-based) ids — the values stored inside `chromR` routes and used to index `distance_table`:**
- `1` = the **depot** (node 1).
- `2 .. nb_customer+1` = the **customers**. Customer with logical id `c` (0-based, `0..nb_customer-1`) has distance-table id `c+2`. Equivalently: a `tsp` entry `t` becomes route node `t+1` (because `tsp` holds `c`, and MATLAB writes `tsp(...)+1`; note the depot occupies id 1, so customer logical id `c` maps to table id `c+1` when `c` is already 1-based... see the concrete rule below to avoid confusion).
- `nb_customer+2 .. ` = the **AFS (refueling / alternative-fuel-station) nodes**.

**Concrete, unambiguous mapping rules actually used in the code:**
- A route node value `v` satisfies `v == 1` ⟺ depot.
- A route node value `v` with `2 <= v <= nb_customer+1` ⟺ a customer (this is confirmed by `chromR{i}(chromR{i} > nb_customer + 1)` being used to select "AFS nodes only", and `time_window` treating `aaa(j+1) ~= 1` (i.e. not depot) as "add customer service time").
- A route node value `v` with `v > nb_customer + 1` ⟺ an AFS node.
- The **generic / first AFS** node id is `nb_customer + 2` (used in `split_Dmax` line 43 as the AFS to test-detour to).
- In `split_Dmax`, a **per-route private AFS** for route number `j` is given id `nb_customer + j + 1`. So route 1's private AFS = `nb_customer+2`, route 2's = `nb_customer+3`, etc. (This matches the generic AFS being `nb_customer+2` for the first route.)
- `distance_table` is a full square symmetric matrix of size `(nb_customer + 1 + nb_AFS) × (nb_customer + 1 + nb_AFS)`, 1-indexed, giving distance between any two nodes. `distance_table(a,b)` = distance from node `a` to node `b`.

### 1.2 The `vrp` struct (input) — fields consumed by these two functions

All are read-only inputs. Shapes/meanings:

| Field | Shape | Meaning |
|---|---|---|
| `vrp.nb_customer` | scalar int | Number of customers, `N`. |
| `vrp.distance_table` | `M×M` double, 1-indexed | Pairwise distances between all nodes (depot, customers, AFS). See 1.1 for id layout. Note: in `Main_METS` distances are scaled by `floor(100*d)/100` before use; that scaling is applied once at load time and is already baked into `vrp.distance_table` by the time Split runs. |
| `vrp.V_speed` | scalar double | Vehicle travel speed. Travel **time** on an edge = `distance / V_speed`. |
| `vrp.V_Dmax` | scalar double | Maximum driving **distance** (range) per fuel tank, i.e. per leg between depot/AFS refuels. The distance splitting bound. |
| `vrp.T_max_V` (aliased `Tmax`) | scalar double | Maximum total **time** allowed per vehicle/route (the duration bound). The time splitting bound. |
| `vrp.T_Customer` | scalar double | Service time spent at each customer. |
| `vrp.T_Afs` | scalar double | Refueling time at an AFS. **Read but NOT used** in either Split function (loaded into a local and ignored). Document it, but it has no effect here. |
| `vrp.T_Start` | scalar double | Start time. **Read but NOT used** in either Split function. |
| `vrp.last_F_location` | scalar/int | **Read but NOT used** in either Split function. |

> Fields `T_Afs`, `T_Start`, `last_F_location` are dead in these two routines — a faithful port may omit them from the Split code entirely. They matter elsewhere in METS, not here.

### 1.3 Inputs / outputs of the Split functions

**Input `tsp`:** a 1×N (or length-N) integer array — a permutation of the customer logical ids `0 .. nb_customer-1`. (In `Main_METS` it is a row `tsp_all(i,:)` or a crossover offspring permutation.)

**Input `par_hgs`:** the third argument (`~`) is **ignored** in both functions.

**Output `chromR`:** a cell array (→ in C++, a `vector<vector<int>>`), one entry per route. Each route is an integer array of **distance-table node ids** of the form:
```
[ 1, <interior nodes...>, 1 ]      // depot ... depot
```
where interior nodes are customers (`2..nb_customer+1`) and possibly AFS nodes (`> nb_customer+1`). Empty routes are stripped at the end (`chromR(cellfun(@isempty,...)) = []`).

### 1.4 Route representation note

Routes here are simple **explicit node-sequence arrays** `[1 ... 1]`. Elsewhere in METS solutions are also stored as predecessor/successor/routeID arrays, but the Split functions produce and consume only the explicit `chromR{i} = [1 ... 1]` array form. Your C++ Split should return `vector<vector<int>>` in this explicit form.

---

## 2. THE MISSING HELPER: `T_able`

`split_Tmax` calls `[isT_able] = T_able(distance_table, V_speed, T_Customer, a, Tmax)` but **`T_able` has no definition anywhere in the provided codebase** (confirmed by exhaustive filesystem search). It is an external path function. Its contract is fully determined by (a) its call signature, (b) how the boolean result is used, and (c) the identical time model defined in `chromR_detail_all.m` line 124. **You must implement it as follows.**

### `T_able(distance_table, V_speed, T_Customer, a, Tmax) -> bool`

- `a` is a prefix of `tsp`: `a = tsp(1:ii)`, i.e. the first `ii` customer logical ids under consideration for the current route.
- It answers: **"If we make a single route serving exactly the customers in `a` (in that order), depot→a(1)→…→a(end)→depot, is that route's total time ≤ Tmax?"**
- Returns `1` (feasible / "able") if time ≤ Tmax, `0` otherwise.

**Exact time model to use (must match `chromR_detail_all` line 124):**
```
route_nodes = [1, a[0]+1, a[1]+1, ..., a[k-1]+1, 1]      // depot, customers (logical+1), depot
total_distance = sum over consecutive pairs of distance_table(route_nodes[p], route_nodes[p+1])
num_customers  = k                                        // = numel(a)
time = total_distance / V_speed + num_customers * T_Customer
return (time <= Tmax) ? 1 : 0
```
Notes:
- Service time is charged once per customer (`num_customers = length(a)`), not for the depot. This is exactly `(numel(a)-2)*T_Customer` in `chromR_detail_all` where there the route already includes both depot endpoints (`numel(a)-2` interior = customers); here we count customers directly as `length(a)`.
- Use `<=` (inclusive) for feasibility. Tie (time exactly Tmax) counts as feasible/able.
- No AFS insertion, no refuel time, no start time in this check — it is a pure duration feasibility test of a depot-to-depot customer chain.

---

## 3. `split_Tmax` — FULL SPECIFICATION

### 3.1 Purpose
Greedily cut the customer permutation `tsp` into consecutive contiguous segments (routes), making each route as long as possible while staying within the **time** bound `Tmax`. First-fit / greedy longest-prefix packing. AFS nodes are **not** inserted by this routine.

### 3.2 Signature
`chromR = split_Tmax(vrp, tsp, ~)` → `vector<vector<int>> split_Tmax(const Vrp& vrp, const vector<int>& tsp)`

### 3.3 Local variables
- `split_tsp` : a working **copy** of `tsp` that gets consumed (elements removed from the front) as routes are cut off. IMPORTANT: keep `tsp` itself intact — the final route reconstruction (step 3.5) indexes the ORIGINAL `tsp`, not `split_tsp`.
- `c` : integer array, length initialized to `numel(tsp)`, all zeros. `c(i)` = number of customers assigned to route `i`. Later, zero entries are removed.
- `chromR` : output cell array, capacity `numel(tsp)` initially (over-allocated; trimmed at end).

### 3.4 Segmentation loop (produces the segment sizes `c`)

Outer loop runs `i = 1 .. numel(tsp)` (i.e. up to N potential routes; almost always breaks out earlier via the `while` consuming `split_tsp`).

```
for i = 1 .. length(tsp):            // route index
    ii = 0
    while length(split_tsp) != 0:
        ii = ii + 1
        a = split_tsp[1 .. ii]                       // first ii customers of the remaining permutation
        isT_able = T_able(distance_table, V_speed, T_Customer, a, Tmax)
        if isT_able == 0:
            // adding the ii-th customer broke the time bound → route i takes the first (ii-1) customers
            c[i] = ii - 1
            remove the first (ii-1) elements from split_tsp     // split_tsp(1:ii-1) = []
            ii = 0
            break                                    // break the while; move to next route i
        if length(split_tsp) == ii:
            // all remaining customers fit in one route → last route
            c[i] = ii
            split_tsp = []                           // empty it → subsequent i-iterations do nothing
            // (no break; while condition length==0 ends it)
    // when split_tsp is empty, inner while does nothing for remaining i
```

**Critical edge cases and tie-breaks in this loop:**
- The prefix is grown one customer at a time; the FIRST prefix that is NOT time-feasible determines the cut. Route `i` gets `ii-1` customers (the last feasible prefix).
- **Degenerate case — a single customer already exceeds Tmax:** if even `a = split_tsp(1:1)` is infeasible on the very first inner iteration (`ii==1`), then `c[i] = 0` and `split_tsp(1:0) = []` removes nothing. `split_tsp` is NOT shortened → the same customer is re-examined on the next outer `i`, again yielding `c=0`, and so on until `i` runs out. The resulting `c` entries are all zero for these and get deleted in step 3.5 (`c(c==0)=[]`), so **that customer is silently dropped from the output** (it never appears in any route). Reproduce this behavior exactly: a customer that cannot fit any route alone is omitted. (This is a latent bug in the source but must be matched for faithfulness; flag it to the caller if desired.)
- Removal `split_tsp(1:ii-1) = []` deletes the assigned customers from the front, so the next route starts at the next unassigned customer.
- When the remaining set fits entirely (`length(split_tsp)==ii` with feasible), that route absorbs all of them and `split_tsp` is emptied.

### 3.5 Building routes from segment sizes

```
c(c==0) = []                          // drop all zero-size segments
a = [1; cumsum(c)]                    // 1-based cut boundaries: a = [1, c1, c1+c2, c1+c2+c3, ...]
                                      // length = (#routes)+1
for i = 1 .. length(a)-1:             // one route per gap
    if i == 1:
        chromR{i} = [1, tsp(a(i)   : a(i+1)) + 1, 1]
    else:
        chromR{i} = [1, tsp(a(i)+1 : a(i+1)) + 1, 1]
chromR(cellfun(@isempty,chromR)) = []
```

**Exact index semantics (1-based, translate carefully to 0-based C++):**
- `a = [1, s1, s2, ...]` where `s_k = c1+...+ck` is the cumulative count.
- Route 1 uses ORIGINAL `tsp` positions `a(1) .. a(2)` = `1 .. c1` inclusive → the first `c1` customers.
- Route `i>1` uses `tsp` positions `a(i)+1 .. a(i+1)` inclusive = customers `(s_{i-1}+1) .. s_i`. The `+1` on the lower bound (only for `i>1`) prevents re-including the boundary customer that route 1 already implicitly consumed via the inclusive `a(1)..a(2)`. Net effect: routes partition `tsp` into consecutive blocks of sizes `c1, c2, c3, …` with NO overlap and NO gap.
- Each selected customer logical id `t` is converted to route node id `t+1` (0-based logical → distance-table id per §1.1). Wrap with depot `1` on both ends.

**Equivalent clean C++ formulation (recommended — avoids the off-by-one branch):** After computing the non-zero segment sizes `c = [c1, c2, ...]`, walk a running pointer `p = 0` over ORIGINAL `tsp`; for each segment size `s`, route = `[1] + [tsp[p+0]+1, ..., tsp[p+s-1]+1] + [1]`, then `p += s`. This produces identical output to the MATLAB `a`/`cumsum` scheme.

### 3.6 Output
`vector<vector<int>>` of routes, each `[1, ..., 1]`, customers only (no AFS). Empty routes already excluded.

---

## 4. `split_Dmax` — FULL SPECIFICATION

### 4.1 Purpose
Cut `tsp` into routes bounded by the **distance/range** bound `V_Dmax`, and — unlike `split_Tmax` — **conditionally insert exactly one AFS (refuel) node per route** when the range would otherwise be exceeded. This is the "Conditional AFS Insertion" of the segmentation phase. It builds routes incrementally, one customer at a time, tracking whether the current route has already had an AFS inserted.

### 4.2 Signature
`chromR = split_Dmax(vrp, tsp, ~)` → `vector<vector<int>> split_Dmax(const Vrp& vrp, vector<int> tsp)`
Note `tsp` is modified inside (elements removed); pass by value/copy.

### 4.3 State variables
- `chromR` : output cell array, over-allocated to `numel(tsp)`; trimmed at end.
- `insertnode` : current route's interior node list as **logical/route-value tokens** being accumulated (starts empty; grows as `[pre_insertnode, tsp(z)]`). NOTE these are the raw `tsp` customer values and the AFS token `nb_customer+j` (see below), NOT yet `+1`-shifted. The `+1` shift to distance-table ids happens when forming `route`/`pre_route`.
- `pre_insertnode` : the previous iteration's `insertnode` (the last accepted/committed interior list).
- `pre_route` : `[0, pre_insertnode, 0] + 1` → the committed route as distance-table ids (depot=1 at both ends). This is the "last known-good" route.
- `route` : `[0, insertnode, 0] + 1` → the candidate route (with the newly added node) as distance-table ids.
- `insert_f` : flag, "an AFS was just inserted on the previous iteration → this iteration must materialize the AFS token into `pre_insertnode` and re-process the same `tsp` element". 0/1.
- `has_f` : flag, "the current route already contains an AFS" (0 = no AFS yet, 1 = AFS present). Governs which distance test is applied.
- `j` : current route number (1-based). Used both as chromR index and to compute the private AFS id `nb_customer + j` / `nb_customer + j + 1`.
- `z` : index into `tsp` of the next customer to consider (1-based, advanced by `z=z+1` each iteration, but rolled back by `z=z-1` right after an AFS insertion so the same customer is retried).
- Loop counter `i` runs `1 .. nb_customer*2` (upper bound = twice the number of customers, giving room for AFS insertions; the loop always terminates earlier via `break` when `tsp` is exhausted).

### 4.4 AFS token numbering inside this function
- While building, the AFS token appended to the interior list is `nb_customer + j` (line 23: `pre_insertnode = [insertnode(1:end-1), nb_customer + j]`). After the global `+1` shift (`route = [0 insertnode 0]+1`) this becomes distance-table id `nb_customer + j + 1` = the private AFS for route `j`. Consistent with §1.1.
- The **generic AFS** used in the range-test detour is `nb_customer + 2` directly as a distance-table id (line 43), i.e. it's already a table id, no shift. For route 1 (`j=1`) this equals the private AFS id `nb_customer+j+1`.
- Terminal fix-up (lines 97–102): any route that ended up as exactly `[1, x, 1]` (length 3, a single node between depots) is rewritten to `[1, a, i+nb_customer+1, 1]` where `a = max(chromR{i})` (the single interior node) — i.e. a lone-customer route gets a trailing private AFS appended. See §4.9.

### 4.5 Per-iteration distance bookkeeping

Each iteration first (optionally) materializes a just-inserted AFS, then appends the next `tsp` customer to form the candidate `route`, then computes:

```
// 1) materialize pending AFS insertion, or carry forward
if insert_f == 1:
    pre_insertnode = [ insertnode[1 .. end-1], nb_customer + j ]   // replace the last appended (customer) with route-j AFS token
    z = z - 1                                                       // retry the SAME tsp customer next
    insert_f = 0
else:
    pre_insertnode = insertnode

// 2) append next customer (tsp(z)) as the trial node
insertnode = [ pre_insertnode, tsp(z) ]
pre_route  = [0, pre_insertnode, 0] + 1        // committed route (distance-table ids)
route      = [0, insertnode,     0] + 1        // candidate route (distance-table ids)

// 3) total candidate route distance D
D = sum_{p} distance_table(route[p], route[p+1])       // over all consecutive pairs

// 4) D_hou = distance of the SUFFIX from the max-id node to the end
[~, maxx] = max(route)                                  // position of the largest node id in route
D_hou = sum_{p = maxx .. end-1} distance_table(route[p], route[p+1])

// 5) D_f = distance if we detour to the generic AFS instead of returning to depot at the end
if has_f == 0:
    D_f = D - distance_table(route[end-1], route[end])          // remove last edge (…→depot)
              + distance_table(route[end-1], nb_customer+2)     // add edge (…→generic AFS)
```

**Interpretation of these three quantities:**
- `D` = total distance of the full candidate route depot→…→depot (with any AFS already in it).
- `D_f` (only meaningful when `has_f==0`, i.e. no AFS yet) = the distance of the route if, instead of the final customer→depot edge, the vehicle detoured from the last-before-depot node to the generic AFS. It's a look-ahead: "if this last leg is too long, would routing via an AFS keep us in range for that partial leg?"
- `D_hou` ("hou" = 后 = "after") = the distance of the route SUFFIX starting at the position of the highest-id node. Because AFS nodes have the highest ids (`> nb_customer+1`), `maxx` is the position of the AFS (when one is present). So when `has_f==1`, `D_hou` is the distance from the AFS to the end of the route = the range consumed on the post-refuel leg. (When no AFS is present, `maxx` is just the position of the highest-id customer — but `D_hou` is only consulted in the `has_f==1` branch, so that case is irrelevant.)

### 4.6 Decision logic — the core state machine

Two top-level branches on `has_f`. Preserve every branch, threshold, and `continue`/`break`/reset exactly.

#### Branch A: `has_f == 0` (route has no AFS yet — test the AFS-detour distance `D_f`)

```
if D_f <= V_Dmax:
    // The route (with a possible end-detour to AFS) is still within range.
    if (length(route) - 2) == length(tsp):
        // route interior count == number of remaining tsp customers → all remaining customers are on this route
        if D <= V_Dmax:
            chromR{j} = route            // full route fits WITHOUT needing an AFS
            break                        // DONE — entire tsp placed
        else:
            // The full return-to-depot distance D exceeds range, but D_f (via AFS) was OK →
            // append a private AFS then depot: replace final depot with [AFS_j+1, depot]
            chromR{j} = [ route[1 .. end-1], nb_customer + j + 1, 1 ]
            j = j + 1
            (if length(tsp)==0: print "splitD")     // debug print only; no logic effect
            tsp = []                                 // consume all remaining
            break                                    // DONE
    // else (not all customers placed yet): the added customer is fine → keep it and continue growing
    continue

elif D_f > V_Dmax:
    // Even detouring to an AFS at the end would exceed range → we must INSERT an AFS now
    insert_f = 1                          // next iteration materializes the AFS token into the route
    has_f   = 1                           // route now considered to have an AFS
    continue                              // note: z is NOT advanced here effectively (next iter does z=z-1)
```

Key points for Branch A:
- Acceptance of a newly added customer with no AFS needed = "`D_f <= V_Dmax` and not the final-customer case" → `continue` (customer stays, loop grows the route).
- The `(length(route)-2) == length(tsp)` test detects "this candidate route now contains every remaining customer" (interior length equals remaining count). `length(route)-2` strips the two depot endpoints.
- When range is exceeded even with the AFS detour, we don't drop the customer — we set flags so the NEXT iteration inserts the AFS BEFORE this customer and retries the customer after the refuel (`z=z-1`).

#### Branch B: `has_f == 1` (route already has an AFS — test the post-refuel suffix distance `D_hou`)

```
if D_hou <= V_Dmax:
    // post-AFS leg still in range
    if (length(route) - 3) == length(tsp):
        // interior count minus the AFS (hence -3 not -2) == remaining customers → all remaining placed
        chromR{j} = route
        j = j + 1
        (if length(tsp)==0: print "splitD")   // debug only
        tsp = []                               // consume all remaining
        break                                  // DONE
    continue                                   // keep growing the post-AFS leg

elif D_hou > V_Dmax:
    // post-AFS leg now exceeds range → the CURRENT (pre-add) route is complete; commit it
    chromR{j} = pre_route                       // commit the last good route (WITHOUT the customer that broke it)
    // remove the customers that were placed on this route from tsp
    for each node token in pre_insertnode:
        remove from tsp all elements equal to insertnode(that index)   // see note below
    j = j + 1
    has_f = 0
    z = 0                                        // restart customer pointer scan for the fresh route
    insertnode = []; pre_route = []; pre_insertnode = []
    // NO continue/break → falls through to next for-iteration, which starts a brand-new empty route
```

Key points for Branch B:
- `(length(route)-3) == length(tsp)`: interior of `route` = customers + 1 AFS. Subtract 2 depots and 1 AFS = number of customers on the route; compare to remaining `tsp` count.
- When the post-AFS leg overflows, we DO NOT put the just-added customer on this route. We commit `pre_route` (the route as it was before adding this customer), then delete the committed customers from `tsp`, bump the route number, reset `has_f`, and reset `z=0` so the next new route re-scans from the current front of `tsp`.

**The `tsp` removal loop (lines 84–86) — exact semantics and a hazard:**
```
for iii = 1 .. length(pre_insertnode):
    tsp( tsp == insertnode(iii) ) = []        // remove ALL elements equal to insertnode(iii)
```
- It iterates over `pre_insertnode` indices but removes values `insertnode(iii)` (value-based deletion, all matches). Since `insertnode = [pre_insertnode, tsp(z)]`, `insertnode(1..end-1) == pre_insertnode`, so `insertnode(iii)` for `iii<=length(pre_insertnode)` equals `pre_insertnode(iii)`. Effect: **remove from `tsp` every customer that was committed to this route.**
- One entry of `pre_insertnode` may be an AFS token `nb_customer+j`. Deleting by that value from `tsp` matches nothing (AFS ids aren't in `tsp`), harmless.
- **Hazard to reproduce faithfully:** deletion is by VALUE and removes ALL matches. `tsp` is a permutation (unique values), so this is safe here. In C++, implement as: build the set of committed customer values from `pre_insertnode` (excluding AFS tokens), then remove those values from the `tsp` working list.

### 4.7 The `insert_f` retry mechanism (how an AFS gets physically inserted)
This is the subtle part; trace it precisely:
1. In Branch A, when `D_f > V_Dmax`, set `insert_f=1`, `has_f=1`, `continue`. At this point `insertnode` still ends with the just-tried customer `tsp(z)`.
2. Next iteration top: `z=z+1` runs first (as always). Then because `insert_f==1`: `pre_insertnode = [insertnode(1:end-1), nb_customer+j]` — i.e. drop that last customer and append the route-`j` AFS token instead; and `z=z-1` (undo this iteration's increment AND the previous, netting the same `z` so the SAME customer is retried); set `insert_f=0`.
3. Then `insertnode = [pre_insertnode, tsp(z)]` re-adds the same customer, now AFTER the AFS. `route` now = depot → (customers so far) → AFS → (this customer) → depot.
4. Now `has_f==1`, so Branch B governs, testing `D_hou` (the AFS→…→depot suffix) against range.

Net effect: an AFS is placed in the route immediately before the customer whose addition would have blown the range, and the range accounting restarts from the AFS.

### 4.8 Loop termination
- Normal termination is via one of the `break` statements (all remaining customers placed) or by exhausting `tsp` and committing routes.
- Hard cap: the `for i = 1 .. nb_customer*2` bound guarantees termination even in pathological cases. In C++ keep an equivalent iteration cap of `2*nb_customer` as a safety bound.
- The `disp "splitD"` prints are diagnostics only — no functional effect; omit or log at debug level.

### 4.9 Post-processing (lines 96–102)
```
chromR(cellfun(@isempty,chromR)) = []        // drop empty route slots
for i = 1 .. length(chromR):
    if length(chromR{i}) == 3:               // route is exactly [1, x, 1] (single interior node)
        a = max(chromR{i})                   // = x, the single interior node id
        chromR{i} = [1, a, i + nb_customer + 1, 1]   // append a private AFS (id i+nb_customer+1) before final depot
```
- Any route that came out as a single-customer route `[1, x, 1]` is rewritten to `[1, x, AFS_i, 1]` where the AFS id is `i + nb_customer + 1` using the FINAL route index `i` (after empties were removed). This guarantees every singleton route has a trailing refuel. Reproduce exactly, using the post-compaction 1-based route index `i` to compute the AFS id (in 0-based C++: `AFS_id = (route_index_0based + 1) + nb_customer + 1`).

### 4.10 Output
`vector<vector<int>>` of routes `[1, ..., 1]`, each possibly containing exactly one AFS node (id `> nb_customer+1`) inserted at the point the range bound required a refuel, plus the singleton fix-up AFS. Empty routes excluded.

---

## 5. OBJECTIVE, PENALTIES, LOCAL SEARCH, CONTROL FLOW — SCOPE NOTE

The two files in scope (`split_Tmax`, `split_Dmax`) are **pure constructive segmentation routines**. They compute NO objective value, NO biased fitness, and apply NO penalty weights. Their only "cost/feasibility" logic is the inline range test (`D`, `D_f`, `D_hou` vs `V_Dmax`) and the time test (`T_able` vs `Tmax`) described above; those are hard construction bounds, not weighted penalties. They contain no local-search operators, no acceptance loop, no population/parent/survivor logic.

For completeness of the port, the surrounding context (NOT in these two files, do not implement from this section — noted only so the C++ author knows where Split sits):
- **Objective/penalties** (`chromR_detail_all.m`): after Split, the full solution is evaluated with weighted penalties `Penalty_all = [wT, wC, wD, wM]`:
  - `penalty_m = max(nb_V - V_nb, 0)`; contributes `wM * penalty_m` (vehicle-count penalty PM).
  - `penalty_D_v = max(leg_distance - V_Dmax, 0)` summed over legs; contributes `wD * sum(penalty_D_v)` (distance penalty PD).
  - `penalty_T_v = max(time_V - T_max_V, 0)` with `time_V = distance/V_speed + (#customers)*T_Customer` per route (time penalty PT).
  - AFS-capacity penalty PC is computed via `get_pc_now`/`AFSdelay_new` from AFS arrival times, `T_Afs`, `C_Afs` (overlap of refuel windows against station capacity).
- **Control flow** (`Main_METS.m`): each individual is built by choosing `split_Dmax` with probability `split_prob` (else `split_Tmax`) in the init phase, and 0.5/0.5 in the main loop; then `chromR_detail_all` → `ELS_mian` local search → `PopManagement` → best-sol/Last100 update → conditional `Repair_sol` (triggered when infeasible and `rand-0.5 > 1e-6`, i.e. ~50% probability). Termination by `maxIter`, `maxIterNonProd`, or `timeLimit`.

These belong to other files and are out of scope for this Split spec; implement them from their own sources.

---

## 6. IMPLEMENTATION CHECKLIST (both functions)

1. Node-id mapping: `tsp` value `t` → route id `t+1`; depot=1; AFS ids `> nb_customer+1`; generic AFS = `nb_customer+2`; route-`j` private AFS = `nb_customer+j+1`. Get this exactly right — it drives every `distance_table` lookup.
2. Time model (for `T_able`): `time = totalDistance/V_speed + numCustomers*T_Customer`; feasible iff `time <= Tmax` (inclusive).
3. Distance tests use raw `distance_table` sums vs `V_Dmax` (inclusive `<=`).
4. `split_Tmax`: greedy longest time-feasible prefix; partition original `tsp` by cumulative segment sizes; drop zero-size segments; **preserve the "single over-long customer gets silently dropped" behavior** (or explicitly document/deviate if the caller prefers — but the faithful port drops it).
5. `split_Dmax`: incremental build with the `has_f`/`insert_f` state machine; one AFS inserted per route when range demands; commit `pre_route` on overflow; remove committed customers from working `tsp`; reset `z=0` for each new route; singleton-route AFS fix-up at the end using the post-compaction index.
6. Strip empty routes at the end of both.
7. Ignore the 3rd argument, `T_Afs`, `T_Start`, `last_F_location` in both — they are dead here.

---

**Source files read in full:**
- `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Two novel segmentation/split_Tmax.m`
- `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Two novel segmentation/split_Dmax.m`

**Supporting files consulted for conventions (node ids, time/distance model, driver call sites):**
- `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Load/chromR_detail_all.m` (time model line 124, AFS-node selection line 75, penalty formulas)
- `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Load/get_vrp.m` (vrp struct source)
- `/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm/Main_METS.m` (call sites lines 86–89, 134–137; distance scaling lines 244–248)

**Important finding for the C++ author:** `T_able` (called by `split_Tmax`) has **no definition anywhere in the codebase** — it is an external MATLAB-path function. Section 2 above reconstructs its exact contract from the call signature, its boolean usage, and the confirmed time model in `chromR_detail_all.m`. Implement `T_able` as specified in §2; do not assume it exists elsewhere.