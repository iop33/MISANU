I have read the entire METS codebase. Below is the exhaustive, implementation-ready specification.

======================================================================
METS (GrVRP-PCAFS) — IMPLEMENTATION-READY C++ SPECIFICATION
======================================================================

Source read in full: Main_METS.m, Load/get_vrp.m, Load/chromR_detail_all.m, plus all supporting files needed to make the port self-contained: split_Dmax.m, split_Tmax.m, ELS_mian.m, phrase_chromR.m, get_chromR.m, get_pd_pt.m, get_pc_now.m, AFSdelay_new.m, AFSdelay_recursion.m, deleteAFS_delta.m, deleteAFS_node.m, m1..m9.m, Depot_m1/2/3/8/9.m, NewRoute_m1/2/3.m, Repair_sol.m, PopManagement.m, add2Pop.m, bestsolfind.m, infeasiblePop_updateBiasedFitnesses.m, update_feasiblePop.m, update_infeasiblePop.m, uti_addSol2Last100.m, uti_updateBestSol.m, selectparents.m, Crossover.m, GrVRP_PCAFS_MILP.py, and a dump of an instance .mat.

Absolute source dir: /Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm

IMPORTANT PORTING NOTE (per your instructions): the MATLAB operators maintain O(1) incremental deltas of distances/times/penalties. In C++ you will FULLY RE-EVALUATE each candidate route(s) after applying the move. So for each operator I give (a) the exact candidate-move set and ordering, (b) what the solution looks like after the move, (c) the acceptance rule, and (d) the AFS-insertion (CRI) rule. You do NOT need to reproduce the delta bookkeeping; re-run the evaluation of chromR_detail_all-style route metrics on the mutated solution and compare full cost.

----------------------------------------------------------------------
0. PROBLEM & NODE INDEXING (CRITICAL)
----------------------------------------------------------------------
Green VRP with Private Capacitated Alternative Fuel Stations. A single physical AFS (refueling station) exists. Each route may visit the AFS at most once (splitting the route into a "pre" segment d1 and a "post" segment d2). Fuel range is bounded by V_Dmax per segment; route duration by T_max_V; the AFS has capacity C_Afs (number of vehicles that may refuel simultaneously); optionally a vehicle-count penalty vs V_nb.

MATLAB is 1-indexed. Two coordinate systems are used and you MUST keep both:

(A) "distance_table space" (1-indexed node IDs used to index vrp.distance_table):
    - node 1              = depot D
    - nodes 2 .. nb+1     = customers C1..Cnb   (customer c -> row/col c+1)
    - nodes nb+2 .. 2*nb+2 = AFS clones. There are (nb+1) identical clone columns,
      all equal to the single physical AFS location. The AFS clone used by route j
      is distance-table index (nb+1 + j)  … but ALL clones are identical, so in
      practice `dAll(anything, nb+2)` (the first clone, index nbClients+2) is used
      everywhere as "distance to AFS". distance_table is symmetric, size
      (2*nb+2) x (2*nb+2). Depot->AFS and customer->AFS use the AFS column.

(B) "move space" / "chromR_move space" (0-indexed-ish node labels used inside ELS
    and populations): depot = 0, customer c = c (1..nb), AFS-of-route-j = nb + j.
    Conversion to distance-table index = label + 1  (so customer c is dAll(c+1,·),
    AFS label nb+j maps to dAll(nb+j+1,·); since all AFS clones are identical this
    equals the AFS column). Depot label 0 -> dAll(1,·).

Throughout, `nbClients = nb = vrp.nb_customer`. "AFS node value" in move space is any
label > nbClients. In chromR (route space, see §1) the AFS for route i is coded as
value `nb + 1 + i` (1-indexed, because depot=1 there).

----------------------------------------------------------------------
1. DATA STRUCTURES
----------------------------------------------------------------------

1.1 vrp (instance; loaded from .mat; treat as immutable inputs except the two
    ALL_brokenDIS caches which the algorithm writes back):
  Fixed scalar/vector fields (confirmed by .mat dump):
  - id[1+nb+1]            : cellstr labels 'D','C1'..'Cnb','S1'. (informational)
  - type[1+nb+1]          : 'd','c'..,'f'. (informational)
  - longitude[1+nb+1], latitude[1+nb+1] : coords (depot, customers, AFS). AFS coord
                            is the last entry. (Used only to build distance_table.)
  - nb_customer          : integer nb.
  - last_customer        : = nb+1  (distance-table index of last customer = nb+1).
  - last_F_location      : = nb+2  (1-index of first AFS clone in id list / AFS col).
  - last_D               : = 1 (depot index). (informational)
  - distance_table       : double, symmetric, (2*nb+2)x(2*nb+2). Euclidean distances.
                           Layout per §0(A). Self-distance 0 on diagonal.
  - V_fuel               : tank capacity Q (=32). (Not used by METS heuristic directly;
                           range is expressed via V_Dmax.)
  - V_fuel_rate          : r (=0.2) fuel per distance. (Not used directly by METS;
                           V_Dmax already encodes Q/r = 32/0.2 = 160.)
  - V_Dmax               : max distance per fuel segment (=160). PC/PD use this.
  - V_speed              : speed (=40). time = distance / V_speed.
  - V_nb                 : nominal fleet size (used for penalty_m; = nb for 15-set,
                           smaller for larger sets, e.g. 30 for 100).
  - T_max_V              : Tmax, max route duration (=7).
  - T_Start              : depot start time (=0). (Not used in cost.)
  - T_Afs                : AFS refuel service time (=0.5).
  - T_Customer           : per-customer service time (=0.5).
  - C_Afs                : AFS simultaneous-refuel capacity (=1,2,8,…).
  - correlatedVertices   : uint, nb x 5. Row c = the 5 granular neighbor CUSTOMERS of
                           customer c (move-space labels 1..nb), ordered nearest-first.
                           This is the granular candidate list for local search
                           (NOT length nbGranular; nbGranular=20 is used only as a
                           reshuffle modulus — see §5.3).
  Mutable caches (created/updated at runtime, sized popSizeLambda+popSizeMu+1):
  - ALL_brokenDIS          : infeasible-pop pairwise broken-pair distance matrix.
  - ALL_brokenDIS_feasible : feasible-pop pairwise broken-pair distance matrix.

1.2 Route representation. Three equivalent encodings coexist:

  (a) chromR (route/"1-indexed depot" space): cell array; chromR{i} is a row vector
      of distance-table node IDs for route i, always starting and ending at depot 1:
      chromR{i} = [1, ...customers/afs..., 1]. A customer c appears as c+1. The route's
      AFS appears as (nb+1+i). Example single-customer route with AFS:
      [1, c+1, nb+1+i, 1].

  (b) chromR_move (move space): chromR from (a) minus the leading/trailing depot and
      with every value decremented by 1 (so depot 0 removed). chromR_move{i} is the
      list of move-space labels visited by route i (customers 1..nb and AFS labels
      nb+i). Built in chromR_detail_all lines 333-341 and ELS_mian lines 6-11.

  (c) Linked-list arrays over move-space node index (the working representation of
      local search). Arrays are length = nbClients + (#AFS currently used). Index i
      = move-space node i (i in 1..nb are customers; i in nb+1.. are AFS nodes, the
      k-th AFS node corresponds to route whose routeID = k):
      - predecessors[i] : previous node in route (0 = depot / none).
      - successor[i]    : next node in route (0 = depot / none).
      - routeID[i]      : route number (1..R) the node belongs to; for AFS bookkeeping
                          the AFS node of route k sits at array position nb+k.
      - node_location[i]: relative position vs the route's AFS:
            -1 = node is BEFORE the AFS (or route has no AFS → whole route is "pre"),
             1 = node is AFTER the AFS,
           100 = the node IS the AFS.
      Built by phrase_chromR (pred/succ/routeID) + the node_location logic at the end
      of chromR_detail_all (lines 362-383). In ELS these four columns are packed as
      Node_related = [predecessors, successor, routeID, node_location] (N x 4).

  Conversion chromR <- linked list: get_chromR(pred,succ,routeID,nb): for each node
  with predecessors==0 (a route start), walk successors until 0, collect the path,
  place at chromR{routeID(start)}. Drop empty routes.

1.3 Route_related (per-route metrics; R x 7 matrix). Columns:
  - [1,2] = penalty_D_v = [pd_pre, pd_su] per route = [max(d1-Dmax,0), max(d2-Dmax,0)].
  - [3]   = penalty_T_v = max(time_v - T_max_V, 0) per route.
  - [4]   = afs_time    = arrival time at the AFS on that route (0 if no AFS). (In some
            patch code column 4 is also (ab)used as a "no-predecessor AFS" flag value 2;
            treat carefully — see §5.4.)
  - [5]   = time_V      = total route duration (incl. service, excl. AFS wait).
  - [6,7] = distance_pre_su = [d1, d2] = distance of pre-AFS and post-AFS segments.
            If no AFS, d1 = whole route distance, d2 = 0.

1.4 sol_table row (per individual). Fields materially used downstream:
  ID, tsp (giant-tour permutation), chromR (route space), chromR_move (move space),
  tsp_now, node_location, predecessors, successor, routeID,
  distance_window, time_window, afs_time, afs_time_delay, delay_duration, isdelay,
  time_V, time_V_shifting, time_V_max/min, time_Total,
  overtime_V, overtime_V_max/min, overtime_Total,
  distance_V, distance_V_max/min, distance_Total,
  distance_pre_su, penalty_D_v, penalty_T_v,
  nb_fueling, time_during, afs_V, afs_Total, nb_V,
  penalty_T, penalty_C, penalty_D, penalty_m  (these are the WEIGHTED penalties),
  cost_Total, IsFeasible, time.
  Population-only added fields: Fitness, fitRank, divRank, avgBrokenDist,
  brokenPairDistance/brokenDist.

1.5 Penalty_all: 2 x 4 matrix.
  Row 1 = weights [wT, wC, wD, wM]. wM is always 0 (see §3).
  Row 2 = current WEIGHTED penalty values [wT·penalty_T, wC·penalty_C, wD·penalty_D,
          wM·penalty_m] for the individual currently being processed.

1.6 Global parameters (Main_METS lines 6-18, 41-49):
  split_prob=0.5; PT=527; PC=195; PD=430; penaltyScaleFactor=1.2;
  penaltyDecreaseFactor=0.85; popSizeMu=154; popSizeLambda=68; targetFeasible=0.2;
  nbLast=20; maxIterNonProd=300; maxIter=2000; timeLimit=100000.
  par_hgs.el=0.5; eliteNum=floor(0.5*154)=77; nc=0.2; nClosest=floor(0.2*154)=30;
  nbGranular=20.
  Penalty_all initial = [[527,195,430,0],[0,0,0,0]].
  (Note a harmless display bug prints PD as PC; the value used is PD=430.)

----------------------------------------------------------------------
2. chromR_detail_all — FULL SOLUTION EVALUATION (reference evaluator)
----------------------------------------------------------------------
This is the canonical "evaluate a chromR (route space) from scratch" routine. Use its
logic as your full re-evaluation for both initial solutions and (adapted) for
re-checking local-search candidates. Inputs: vrp, chromR (route space), weights
Penalty_all row 1. Produces all metrics + node_location + Route_related.

Let wT,wC,wD,wM = Penalty_all(1,1..4).  R = numel(chromR).

STEP A. Vehicle-count penalty:
  nb_V = R.
  penalty_m = max(R - vrp.V_nb, 0).
  weighted: penalty_m_w = wM * penalty_m  (wM=0 → always 0).

STEP B. Per-arc time/distance windows. For each route i, take a = chromR{i} with the
  depot value 1 removed everywhere it equals 1 EXCEPT it iterates over the raw route.
  Actually: aaa = chromR{i}; remove entries == 1 (depot). Then for j=1..len(aaa)-1:
     dist_arc = dall(aaa(j), aaa(j+1))
     if aaa(j+1) != 1 (i.e. arrival node is not depot):
        time_arc = dist_arc / V_speed + T_Customer
     else:
        time_arc = dist_arc / V_speed
  Store into time_window(i, j+1), distance_window(i, j+1). (These per-step arrays are
  used only to compute afs_time cumulatively.) NOTE: because depot 1's are stripped
  from aaa first, the "service time added" logic effectively adds T_Customer for every
  arrival that is a customer or AFS. (The AFS also receives +T_Customer here; the AFS
  refuel time T_Afs is handled separately in the capacity model, not in route duration.)

STEP C. AFS location & arrival time per route. For route i:
  a = chromR{i} entries with value > (nb+1)  → the AFS nodes on the route.
  For each such AFS occurrence:
     afs_location(i) = position (1-based index) of that AFS within chromR{i}.
     afs_time(i)     = sum(time_window(i, 1 : afs_location(i))) - T_Customer.
     (Subtracting one T_Customer because the AFS itself is not a customer service.)
  (At most one AFS per route in valid solutions.)

STEP D. Segment distances & PD. For route i, let z = [1, afs_location(i), len(chromR{i})]
  with zeros removed. For each segment jj (there are (#AFS on route)+1 = 1 or 2):
     D = sum of dall(chromR{i}(k), chromR{i}(k+1)) for k in segment.
     distance_pre_su(i, jj) = D.
     penalty_D_v(i, jj) = max(D - V_Dmax, 0).
  penalty_D = sum(sum(penalty_D_v)); weighted = wD·penalty_D.

STEP E. Route duration & PT. For each route i:
  distance_V(i) = sum of all arc distances of chromR{i} (full route incl. return).
  time_V(i)     = distance_V(i)/V_speed + (len(chromR{i}) - 2)·T_Customer.
     (len-2 = number of intermediate nodes = customers + AFS. So AFS counts as a +T_Customer
      in duration too. This matches STEP B.)
  overtime_V(i) = max(time_V(i) - T_max_V, 0).
  overtime_Total = sum(overtime_V). penalty_T = max(overtime_Total,0); weighted = wT·penalty_T.
  Also penalty_T_v(i) = max(time_V(i) - T_max_V, 0) (per-route, stored in Route_related col3).
  time_V_shifting(i) = max(T_max_V - time_V(i), 0)  (slack available to delay AFS arrival).

STEP F. AFS capacity penalty PC (this is the PCAFS core). The AFS can serve C_Afs
  vehicles at once; refuel takes T_Afs. A vehicle arrives at afs_time(i) and occupies
  the pump on [afs_time(i), afs_time(i)+T_Afs]. If more than C_Afs vehicles overlap in
  any sub-interval, the overlap incurs penalty proportional to (excess count × duration).
  METS additionally allows delaying a vehicle's AFS arrival within its slack
  time_V_shifting to resolve conflicts (afs_time_delay). Two implementations exist and
  MUST match:

  F1. chromR_detail_all's inline delay pass (lines 136-221): pairwise conflict
      resolution loop over route pairs (i,ii). Definitions per pair: order so
      t1=earlier arrival, t2=later arrival; q1,q2 = the corresponding slacks; t3 flags
      which of i/ii is later. Cases:
        - if t1+T_Afs <= t2: no overlap, skip.
        - if (t1+T_Afs-t2 > q2) AND (t2+T_Afs-t1 > q1): unresolvable → conflict_table=1.
        - if only t1+T_Afs-t2 > q2: delay the LATER-finishing one to t2+T_Afs;
          set its afs_time_delay and recompute its slack = T_max_V - t1new - T_Afs;
          restart pair scan (i=1,ii=1).
        - symmetric case for the other side.
        - if both slacks suffice: delay the one with SMALLER slack; if q1>q2 delay t1
          side, else delay t2 side; restart scan.
      Result: afs_time_delay vector (delayed AFS arrival times) and isdelay flag.

  F2. The reusable evaluator used by local search — get_pc_now + AFSdelay_new/
      AFSdelay_recursion (this is the definitive PC computation you should port):

      get_pc_now(afs_time_delay, T_max, time_v, time_afs, T_Afs, C_Afs):
        remove zero entries from afs_time_delay.
        afs_time_end = afs_time_delay + T_Afs.
        c = sort([afs_time_delay ; afs_time_end]).
        time_during = diff(c)   (widths of consecutive time slices).
        For each slice k (left edge c(k)): d = afs_time_delay - c(k);
           overlap = count of vehicles with (eps < d + T_Afs) AND (d <= 0)
                     = # vehicles currently refueling in slice k.
        nb_fueling(k) = max(overlap(k) - C_Afs, 0)   (excess over capacity).
        pc_now = max( sum_k nb_fueling(k)·time_during(k), 0 ).
      RETURNS unweighted PC. (chromR_detail_all uses the same formula manually in
      lines 244-279, dropping the last slice: nb_fueling(end)=[].)

      AFSdelay_new(time_afs, time_v, T_max, T_Afs, C_Afs, afs_time_delay):
        if no AFS: pc_now=0.
        time_V_shifting = max(T_max - time_v, 0).
        afs_time_delay = time_afs (reset).
        call AFSdelay_recursion → returns (pc_now, afs_time_delay).

      AFSdelay_recursion: computes pc via get_pc_now, then runs the same pairwise
      backtracking-delay resolution as F1 (the greedy delay loop), attempting to reduce
      pc by delaying vehicles within slack. It restarts the scan whenever it applies a
      delay (delay_over flag), and recursively tries both delay-1 vs delay-2 options
      keeping whichever yields lower pc (the `backtrack` local function). Final pc_now
      = get_pc_now(final afs_time_delay).
      PORTING GUIDANCE: implement F2 exactly. This is the AFS-capacity penalty
      evaluator; the weighted value is wC·pc_now. It is invoked after every candidate
      move that changes any route's AFS timing.

STEP G. cost_Total = wT·penalty_T + wC·penalty_C + wD·penalty_D + wM·penalty_m
                     + distance_Total   (distance_Total = sum(distance_V)).
        IsFeasible = 1 iff (penalty_T + penalty_C + penalty_D + penalty_m == 0) using
        the UNWEIGHTED penalties; else 0.

STEP H. node_location construction (route space → move space linked list positions),
  lines 362-383: build tsp_now = concatenation of chromR_move routes. Sort A=tsp_now,
  get I. For each route: if route has no AFS (max label <= nb) → all its customers get
  node_location = -1. Else split at afs_loc (AFS position within route, 0-based-ish):
  nodes before AFS → -1, nodes from AFS onward → 1, the AFS itself → 100. (See §1.2d.)

----------------------------------------------------------------------
3. OBJECTIVE & PENALTY FORMULAS (EXACT)
----------------------------------------------------------------------
Weighted penalties (wM≡0):
  PT (time)     = wT · penalty_T,   penalty_T = Σ_i max(time_V(i) − T_max_V, 0).
  PD (distance) = wD · penalty_D,   penalty_D = Σ_i Σ_{seg∈{pre,su}} max(dseg(i) − V_Dmax, 0).
  PC (AFS cap)  = wC · penalty_C,   penalty_C = Σ_slice max(overlap−C_Afs,0)·width  (post-delay, §2.F2).
  PM (vehicles) = wM · penalty_m,   penalty_m = max(R − V_nb, 0).  wM=0 ⇒ PM=0 always.
Objective:
  cost_Total = distance_Total + PT + PC + PD + PM.
Feasibility:
  IsFeasible = (penalty_T==0 && penalty_C==0 && penalty_D==0 && penalty_m==0).
Time model:
  time_V(i) = distance_V(i)/V_speed + (#intermediate nodes)·T_Customer, where
  #intermediate nodes = (len(chromR{i})−2) = customers + AFS on the route.
  (AFS contributes T_Customer to duration; AFS refuel time T_Afs enters ONLY the
   capacity model PC, not route duration.)
Distance-segment model (per route with one optional AFS):
  d1 = depot→…→AFS distance (or whole route if no AFS); d2 = AFS→…→depot (or 0).
  time_afs(i) = d1/V_speed + (#nodes strictly before AFS incl. it appropriately)·T_Customer,
  computed via get_pd_pt (below) for candidate re-eval.

get_pd_pt(vrp, route_now(move-space, no depot), everTime=T_Customer, speed):
  route_now = ordered move-space labels of the route (customers + AFS label), NO depot.
  n = position of max label (the AFS, since AFS labels > customer labels).
  If route has an AFS (max label > nb):
     d1 = Σ arcs up to AFS + dAll(1, route_now(1)+1)   [depot→first + pre segment]
     d2 = Σ arcs after AFS + dAll(route_now(end)+1, 1) [post segment + last→depot]
     t1 = d1/speed + (n−1)·everTime      [pre duration; n−1 nodes before AFS incl?]
     t2 = d2/speed + (len − n + 1)·everTime
  Else (no AFS):
     d1 = Σ arcs + dAll(1,first+1) + dAll(last+1,1);  t1 = d1/speed + len·everTime;
     d2 = 0; t2 = 0.
  (Use +1 offset because move-space label k → distance-table index k+1; depot=1.)
  time_V = t1 + t2; PT_v = max(time_V − T_max_V, 0). PD_v = [max(d1−Dmax,0),max(d2−Dmax,0)].

Penalty adaptation (Main_METS lines 171-239): applied every nbLast(=20) iterations once
tspid ≥ 100, based on the last-nbLast solutions (Last100Sol). Compute fractions of the
last window that were feasible in each dimension separately:
  fractionFeasible_T = (#sols with penalty_T==0)/nbLast, similarly _C, _D.
Rules per dimension X∈{T,C,D}, using targetFeasible=0.2 and ±0.05 band:
  if fractionFeasible_X <= 0.15: wX = min(100000, wX·penaltyScaleFactor(1.2)).
  elif fractionFeasible_X >= 0.25: wX = max(0.1, wX·penaltyDecreaseFactor(0.85)).
  else: unchanged.
After changing wX, rescale every infeasible individual's stored weighted penalty:
  infeasiblePop.penalty_X(k) = (penalty_X(k)/origin_wX)·wX_new  (origin_wX = old weight).
Then infeasiblePop.cost_Total = penalty_D + penalty_C + penalty_T + distance_Total, and
re-sort the infeasible population by biased fitness (see §4.4) and rebuild vrp.ALL_brokenDIS.

----------------------------------------------------------------------
4. POPULATION MANAGEMENT
----------------------------------------------------------------------
Two subpopulations, feasiblePop and infeasiblePop (tables sorted by cost ascending).

4.1 add2Pop (insert one evaluated individual):
  Route to feasiblePop if IsFeasible==1 else infeasiblePop. If empty subpop: initialize
  (Fitness=1, fitRank=1, avgBrokenDist=0, divRank=0, ALL_brokenDIS row of zeros).
  Else: append the new individual at the end; if its cost_Total < worst current cost,
  insertion-sort it into ascending-cost position (shift others down), and correspondingly
  shift the pairwise broken-distance matrix rows/cols. Then:
    - fitRank(i) = (rank_by_cost(i) − 1)/(N−1)   (0 = best cost, 1 = worst).
    - compute broken-pair distance of the new individual vs every other (see 4.3).
    - avgBrokenDist(i) = −mean( nClosest smallest broken distances to others ).
    - divRank via sort of avgBrokenDist, normalized to [0,1] (0 = most diverse/isolated).
    - Fitness(i) = fitRank(i) + (1 − eliteNum/N)·divRank(i).   LOWER Fitness = better.

4.2 PopManagement (called after add2Pop): if subpop size > popSizeMu+popSizeLambda
  (=222), trim down to popSizeMu(=154) by repeatedly removing the "worst" individual:
    - A clone (its nearest broken-pair distance ≈ 0, i.e. sum of two smallest ≈ 0) is
      always considered worse than any non-clone. Among same clone-status, remove the
      one with the LARGEST Fitness (worst). Never remove index 1 (the best). After all
      removals, recompute fitRank, avgBrokenDist, divRank, Fitness as in 4.1.
  Update vrp.ALL_brokenDIS / ALL_brokenDIS_feasible accordingly.

4.3 Broken-pair distance between two individuals p,q (add2Pop / update_*Pop /
    infeasiblePop_updateBiasedFitnesses all use the same core):
    Take successors & predecessors arrays (move space) of both, truncate both to the
    common min length. Map every AFS label (>nb) to the sentinel (nb+1) so all AFS are
    identified. For each node jj in 1..len:
      - "extra" mismatch: if succ_p(jj) ≠ succ_q(jj) AND succ_p(jj) ≠ pred_q(jj) → +1.
      - "missing" mismatch: if pred_p(jj)==0 AND pred_q(jj)≠0 AND succ_q(jj)≠0 → +1.
    distance = (total mismatches) / vrp.last_customer  (= /(nb+1)); note update_feasiblePop
    and update_infeasiblePop divide by nbClients (nb) instead and iterate to
    min(#AFS_p,#AFS_q)+nbClients — a minor inconsistency; the runtime path used in the
    main loop is add2Pop (÷last_customer) and infeasiblePop_updateBiasedFitnesses
    (÷last_customer). Reproduce add2Pop's version for the live path.

4.4 infeasiblePop_updateBiasedFitnesses (called during penalty adaptation): re-sort
    infeasible pop by cost_Total (ties: larger ID first), recompute pairwise broken
    distances, avgBrokenDist (−mean of nClosest smallest), divRank = linspace(0,1,N),
    fitRank normalized, Fitness = fitRank + (1−eliteNum/N)·divRank, then final sort by
    Fitness (ties: larger ID first). Writes brokenDist per row (used to rebuild ALL_brokenDIS).

4.5 uti_updateBestSol: if sol is feasible and cost < bestSolRestart.cost − 1e-9, update
    bestSolRestart; if also < bestSolOverall.cost − 1e-9, update bestSolOverall; set
    isNewBest=true; else isNewBest=false. bestSolOverall initialized IsFeasible=0,
    cost_Total=999999.

4.6 uti_addSol2Last100: maintain a FIFO window of the most recent nbLast(=20)
    individuals (despite the name "100"). If full, drop the oldest, append newest.

----------------------------------------------------------------------
5. LOCAL SEARCH (ELS_mian) — DRIVER & NEIGHBORHOODS
----------------------------------------------------------------------
5.1 Entry. Inputs include chromR (route space unless isrepair), node_location,
    Route_related, Penalty_all. If not repairing, convert chromR→chromR_move. Build
    pred/succ/routeID via phrase_chromR. Assemble Node_related = [pred,succ,routeID,
    node_location] and Route_related as in §1.3. afs_time_delay initialized from
    Route_related col 4.

5.2 Book-keeping arrays:
  whenLastTestedRI[nbClients] = 0; whenLastModified[R] = 0; nbMoves=0; loopID=0;
  searchCompleted=false; isSuccess=0.

5.3 Outer loop `while ~searchCompleted`:
  - If loopID>0 set searchCompleted=true at loop top (guarantees ≥2 passes; the flag is
    reset to false by any successful move, so search continues until a full pass makes
    no move).
  - Reseed RNG: rng(SEED+tspid). For each customer i in 1..nbClients: with probability
    1/nbGranular (i.e. if mod(randi(1e9), nbGranular)==0) randomly permute that row of
    correlatedVertices (shuffle the 5 granular neighbors of i).
  - Early time cutoff (only when tspid≠1): if the wall-clock spent on THIS individual
    exceeds sol_table.time(end-1)·20/#individuals, break. (Adaptive per-individual time
    budget; for the C++ port you may keep an equivalent per-individual time cap or omit
    if using iteration caps — but to be faithful, implement it.)

  - For nodeU = 1..nbClients (each customer, in order):
      correlatedU = correlatedVertices(nodeU, :) (the 5 neighbors, possibly shuffled).
      lastTestRINodeU = whenLastTestedRI(nodeU); whenLastTestedRI(nodeU)=nbMoves.
      For posV = 1..5:  nodeV = correlatedU(posV).
        (5.4 "patch" resynchronization runs when whenLastModified changed since last —
         see below.)
        Gate: proceed with moves only if loopID==0 OR
          max(whenLastModified(routeID(nodeU)), whenLastModified(routeID(nodeV))) >
          lastTestRINodeU   (i.e. at least one involved route changed since U was last
          tested — the standard granular "don't retry unchanged pairs" rule).
        If gate passes, TRY operators m1..m9 IN THIS EXACT ORDER; each returns isSuccess.
        On the FIRST success, apply the move (Node_related/Route_related/Penalty_all/
        afs_time_delay are updated in place, nbMoves++, searchCompleted=false,
        whenLastModified updated) and `continue` to next (nodeU-inner) iteration:
            m1 → m2 → m3 → m4 → m5 → m6 → m7 → m8 → m9.
        THEN, only if nodeV == 0 (i.e. node_location(nodeV)==0 meaning… actually the
        code checks `if Node_related(nodeV,1)==0`, i.e. nodeV's predecessor is depot →
        nodeV is a route's first node): additionally try, in order:
            Depot_m1 → Depot_m2 → Depot_m3 → Depot_m8 → Depot_m9.
      After the posV loop, and only when loopID≠1 and nodeV was the LAST correlated
      neighbor (nodeV == correlatedU(end)): try New-route insertions, in order:
            NewRoute_m1 → NewRoute_m2 → NewRoute_m3
      (these create a brand-new route for nodeU / nodeU,x with a fresh AFS).
  - loopID++.

  First-improvement acceptance: the FIRST operator whose delta < −1e-6 is applied; the
  loop then moves on. There is no best-improvement scan.

5.4 The "patch" block (ELS_mian lines 60-144) runs at the start of each (U,V) pair when
  `ispatch != sum(whenLastModified)` (i.e. a modification happened). It re-normalizes
  Route_related/Node_related after route deletions and AFS bookkeeping drift:
    - Drop any Route_related row whose [d1,d2] are both ~0 (empty route) and decrement
      higher routeIDs.
    - For each route, if its AFS is redundant (removing the AFS keeps the merged segment
      ≤ V_Dmax) it deletes the AFS node, recomputes that route's PD/PT via get_pd_pt,
      and re-packs Route_related. It also reorders AFS node rows so the k-th AFS node
      sits at array position nb+k, and recomputes PC via AFSdelay_new when AFS timing
      changed. PORTING GUIDANCE: since you fully re-evaluate each candidate, you can
      replace this entire patch with a single canonicalize() step after every accepted
      move: (i) rebuild chromR via get_chromR, (ii) drop empty routes, (iii) for each
      route, remove the AFS if the whole route without AFS is ≤ V_Dmax (this is the
      "conditional AFS removal" mirror of CRI), (iv) recompute all Route_related metrics
      and Penalty_all row 2 from scratch, (v) rebuild node_location/pred/succ/routeID
      with AFS nodes reindexed to positions nb+1, nb+2, …. This yields identical
      semantics without the delta patch.

5.5 On exit: chromR = get_chromR(...). distance_Total = Σ Route_related[:,6:7].
  cost_Total = distance_Total + Σ Penalty_all(2,:). IsFeasible = (ΣPenalty_all(2,:)==0).
  Write all fields back to sol_table.

--- NEIGHBORHOOD OPERATOR CATALOG (semantics; re-evaluate fully in C++) ---
Common notation: preU=pred(U), X=succ(U), suX=succ(X), Y=succ(V), preV=pred(V),
suY=succ(Y). routeU=routeID(U), routeV=routeID(V). node_location values as §1.2d.
Acceptance for ALL operators: compute the full solution after the tentative move,
delta = cost_after − cost_before; accept iff delta < −1e-6 (MATLAB: `delta > −1e-6`
means reject). All weighted penalties recomputed (PT via per-route max, PD via segment
maxes, PC via AFSdelay_new, PM via route count).

CONDITIONAL AFS INSERTION (CRI) — shared rule for m1/m2/m3 (and analog in others):
When U (or U,X / X,U) is moved into routeV and routeV currently HAS NO AFS
(max node_location of routeV ≠ 100), the operator forms the tentative new routeV and
computes its metrics via get_pd_pt. It ALWAYS inserts a fresh AFS immediately AFTER the
inserted node(s) (route becomes [ …V, U, AFS, Y…] for m1; […V,U,X,AFS,Y…] for m2;
[…V,X,U,AFS,Y…] for m3). d1 = depot→…→AFS, d2 = AFS→…→depot. The new AFS is added as a
new move-space node at index numel(node_location)+1 with node_location=100, its own
routeID=routeV, and node_location of pre-AFS nodes set to −1, post-AFS to +1. The move
is accepted only if its full delta < −1e-6. (There is no separate "try with vs without
AFS" — insertion is unconditional when routeV had no AFS; the acceptance test decides.)
Conversely CONDITIONAL AFS REMOVAL happens implicitly when a route is reduced to only
[AFS] (isdelete): the route+AFS is deleted and delta gets an extra −2·dAll(1,nb+2)
(removing depot↔AFS round trip).

m1 — INSERT U AFTER V.  Move U from its route, place between V and Y.
  Reject if U==Y (already there).
  Same route (routeU==routeV): simple relocation; recompute route.
  Diff route: if routeV has AFS → plain insertion of U into routeV; if routeV has no
  AFS → apply CRI (insert U then AFS). If routeU becomes 2 nodes = [customer, AFS] only
  (isdelete) delete routeU (with −2·dAll(1,nb+2) credit).
  After: routeV = [ …V U (AFS) Y… ]; routeU loses U.

m2 — INSERT (U,X) AFTER V.  Move the pair U,X (X=succ(U)) after V.
  Reject if U==Y or V==X or X==0(depot). Reject if X is an AFS (X>nb) when cross-route.
  Same-route only if X is not the AFS (nodeX_loc≠100). Cross-route with CRI: route
  becomes [ …V U X (AFS) Y… ]. isdelete when routeU had exactly 3 nodes incl. AFS.

m3 — INSERT (X,U) AFTER V (reversed pair).  Same as m2 but the pair is placed reversed:
  route becomes [ …V X U (AFS) Y… ]. Reject if X==0, X==V, or U==Y. Same conditions on
  AFS/isdelete/CRI as m2.

m4 — SWAP U AND V.  Exchange positions of single nodes U and V (any routes).
  Reject if U==preV, U==Y, or U>V (canonical ordering to avoid double-count).
  After: U takes V's slot (pred preV, succ Y), V takes U's slot (pred preU, succ X);
  routeID/node_location swapped. Recompute both routes.

m5 — SWAP (U,X) AND V.  Exchange the pair (U,X) with single node V.
  Reject if U==preV, X==preV, U==Y, or X==0. Same-route requires X not AFS.
  Cross-route: if routeV has AFS → plain swap; if routeV has no AFS → CRI: after placing
  U,X into routeV, insert an AFS (route reconstructed via the routeV_now assembly and
  get_pd_pt). After: routeU gets V; routeV gets U,X (+AFS if inserted).

m6 — SWAP (U,X) AND (V,Y).  Exchange two consecutive pairs across/within routes.
  Reject if X==0,Y==0,Y==preU,U==Y,X==V,V==suX. Reject if BOTH X and Y are AFS, or if
  EXACTLY ONE of X/Y is AFS (only proceeds when neither X nor Y is the AFS, i.e.
  both node_location ≠100). After: swap the two pairs' positions & routeIDs.

m7 — 2-OPT (intra-route). (u,x)(v,y) → (u,v)(x,y): reverse the segment of routeU
  between X and V. Only when routeU==routeV. Reject if any predecessor chain from U back
  reaches V (would break the route), or X==nodeV. If the reversed segment contains the
  AFS (max node_location==100) recompute segments via get_pd_pt (the AFS may move
  between pre/post). After: segment X..V reversed; pred/succ rewired.

m8 — 2-OPT* type A. (u,x)(v,y) → (u,v)(x,y) ACROSS two routes: cut routeU after U and
  routeV before V, reconnect U→V and (tail of U)→(y). Only cross-route (routeU≠routeV).
  Handles the AFS carefully: xx = nodes after U in routeU; vv = nodes from V back to
  route start. Several sub-cases by whether xx/vv contain the AFS:
    - Both segments AFS-free or both with AFS in the "wrong" arrangement → reject.
    - The valid cases move vv into routeU (as the post part) and/or xx into routeV,
      swapping the AFS node ownership (xxafs/vvafs pointer swap) so each resulting route
      keeps exactly one AFS on the correct side. isdelete when a route collapses to only
      its AFS (extra −2·dAll(1,nb+2) or −dAll(1,nb+2) credit depending on branch).
  After: routeU = [depot..U, V..end], routeV = [depot.., X..end] with AFS reassigned.

m9 — 2-OPT* type B. (u,x)(v,y) → (u,y)(x,v): connect U→Y and X→V (the "other" 2-opt*).
  Cross-route only. Builds xx (after U), vv (V back to start), yy (after V). Many
  sub-cases by which of xx/yy contain the AFS; uses get_pd_pt to recompute the two
  rebuilt routes a=[uu yy], b=[fliplr(xx) vv] where uu = U back to route start. Swaps
  AFS node ownership (xxafs/vvafs) as needed. After: routeU = [depot..U, Y..end],
  routeV = [depot.., X.., V..end] with AFS reassigned.

Depot_m1/m2/m3 — same as m1/m2/m3 but INSERT U (or U,X / X,U) at the FRONT of routeV
  (immediately after the depot, before V which is routeV's first node). Triggered only
  when nodeV is a route-start (pred(V)==0). Depot_m2/m3 require routeV to already have
  an AFS (they do NOT do CRI new-AFS insertion; they reject if routeV has no AFS).
  Depot_m1 does support isdelete of an emptied routeU.

Depot_m8 — variant of m8 where V is a route start: relocate the tail segment xx of
  routeU to the front of routeV. Includes an `afs_swap` branch that, after moving,
  physically swaps routeU/routeV rows in Route_related and their routeIDs (so AFS
  ownership/order stays canonical).

Depot_m9 — variant of m9 with V a route start: exchanges the tail of routeU (xx) with
  the whole of routeV (yy), swapping AFS ownership; isdelete when routeV collapses to
  only its AFS.

NewRoute_m1/m2/m3 — create a BRAND-NEW route for U (m1), (U,X) (m2), or (X,U) (m3):
  routeAFS = max(routeID)+1; a new AFS node is appended (node_addAfs). The new route is
  [depot, U, AFS, depot] (m1) / [depot, U, X, AFS, depot] (m2) / [depot, X, U, AFS,
  depot] (m3). Its d1 = depot→U(→X)→AFS, d2 = AFS→depot. PM recomputed as
  max(0, (routeAFS − V_nb))·wM (=0 since wM=0). U (and X) removed from their old route;
  isdelete if the old route collapsed. m2/m3 reject (isSuccess=false) if the old route
  had exactly 3 nodes incl AFS (would create/destroy simultaneously). Accept iff full
  delta < −1e-6.

----------------------------------------------------------------------
6. INITIAL POPULATION CONSTRUCTION (giant tour → split)
----------------------------------------------------------------------
6.1 Generate popSizeMu·4 (=616) random customer permutations (giant tours):
    rng(SEED+1); tsp_all(i,:) = randperm(nbClients) for i=1..616.

6.2 For each i in 1..616 (respecting maxIter and timeLimit breaks):
    tsp = tsp_all(i,:).
    With prob split_prob(0.5): chromR = split_Dmax(vrp,tsp); else split_Tmax(vrp,tsp).
    (In the main loop after crossover the code passes `tsp` — note a bug: it splits the
     stale `tsp` variable, not the fresh offspring — see §7. In the INIT loop it
     correctly uses the current tsp.)
    Evaluate via chromR_detail_all → sol row + node_location + Route_related.
    Run ELS_mian (local search) on it.
    PopManagement (add2Pop + trim). uti_updateBestSol. uti_addSol2Last100.
    Repair (see §8) if infeasible and a coin flip passes.

6.3 split_Tmax(vrp, tsp): greedily cut the giant tour into routes so each route's
    duration ≤ Tmax. It walks prefix a = tsp(1:ii) and calls T_able(distance_table,
    V_speed,T_Customer,a,Tmax) to test whether the prefix route [depot, a, depot] fits
    within Tmax; when adding one more customer would exceed Tmax it closes the route at
    ii−1 and starts a new one. Produces chromR routes [1, tsp(seg)+1, 1] (NO AFS added
    at split time; AFS is introduced later by local search / by the Dmax splitter).
    NOTE: T_able is REFERENCED but its .m file is ABSENT from the repo. Its required
    semantics (infer & implement): given customer sublist a, build route depot→a→depot,
    compute duration = totaldist/V_speed + numel(a)·T_Customer, return
    isT_able = (duration ≤ Tmax) ? 1 : 0. Implement exactly this.

6.4 split_Dmax(vrp, tsp): greedily grow a route customer by customer; before the segment
    distance would exceed V_Dmax it inserts an AFS (has_f flag) splitting the route into
    pre/post around the AFS; when the post segment would also exceed V_Dmax it closes the
    route and starts a new one. AFS for route j is coded as node (nb+1+j) (route space).
    Special close-out: any route that ended up as just [1, c, 1] (3 elements) is rewritten
    as [1, c, (i+nb+1), 1] (single customer + its AFS). Produces chromR with AFS nodes
    already present. (Because the algorithm is delicate, port it literally following the
    has_f / insert_f state machine in split_Dmax.m lines 20-95, then the 3-element
    fixup lines 97-102.)

----------------------------------------------------------------------
7. MAIN GENETIC LOOP
----------------------------------------------------------------------
For tspid = popSizeMu·4+1 (=617) .. maxIter (=2000):
  Termination: break if nbIterNonProd > maxIterNonProd(300) OR toc > timeLimit.
  p1 = selectparents(...); p2 = selectparents(...).   (§7.1)
  offspring_tsp = Crossover(p1,p2,vrp).                (§7.2)
  rng(SEED+tspid); with prob 0.5 split_Dmax else split_Tmax.
    ⚠ FAITHFUL-PORT BUG: the split is called on `tsp` (a stale leftover variable from
      the init loop), NOT on offspring_tsp. So in the reference code the split ignores
      the crossover result. chromR_detail_all is then called with offspring_tsp as the
      recorded tsp, but chromR comes from splitting the stale tour. To reproduce results
      exactly, replicate this (split the last init-loop `tsp`); to fix the intended
      behavior, split offspring_tsp. DOCUMENT which you choose. (Recommend faithful port
      first, then a flagged fix.)
  Evaluate (chromR_detail_all) → ELS_mian → PopManagement → uti_updateBestSol →
  uti_addSol2Last100.
  Repair if infeasible & coin flip (§8).
  Update nbIterNonProd: if isNewBest → nbIterNonProd=1 (and record toBestTime) else ++.
  Penalty adaptation every nbLast iters once tspid≥100 (§3 adaptation).

7.1 selectparents (binary tournament over the COMBINED pool):
    a=|feasiblePop|, b=|infeasiblePop|. Draw p1=randi(a+b): index >a picks
    infeasiblePop(p1−a) else feasiblePop(p1). Draw p2 likewise. Return whichever has
    the SMALLER Fitness (lower Fitness = better). Called twice → p1,p2 (they may be
    equal). RNG note: no explicit reseed here, so uses the ambient RNG stream.

7.2 Crossover (OX-style order crossover on the customer-only giant tours):
    Extract each parent's chromR via get_chromR, flatten to x1,x2, DROP all AFS labels
    (>nb) → pure customer permutations of length nPoint=nb. c=randi(nb,2,1);
    point1=min(c), point2=max(c). y1[point1:point2]=x1[point1:point2] (inherit middle
    from p1). Fill the rest of y1 from x2 in wrap-around order starting after point2,
    skipping customers already in y1. Symmetrically y2 from p2's middle + p1 fill.
    Returns y1 (used as offspring_tsp) and y2 (unused).

----------------------------------------------------------------------
8. REPAIR
----------------------------------------------------------------------
Trigger (both loops): if sol.IsFeasible==0 AND rand−0.5 > 1e-6 (≈ 50% chance).
Repair_sol: set isrepair=1; multiply every penalty column that is currently >0 by
WP=10 (temporarily 10× the weights for the violated dimensions only). Re-run ELS_mian
with isrepair=1 (which reuses the current chromR_move & node_location rather than
re-splitting) to drive the solution toward feasibility under the amplified penalties.
If the repaired solution is now feasible, insert it via PopManagement and update best.
Reset isrepair=0. (Weights are restored implicitly because the ×10 is applied to a local
copy of Penalty_all passed into repair, not the persistent weights.)

----------------------------------------------------------------------
9. FINALIZATION / OUTPUT
----------------------------------------------------------------------
Round the distance_table to 2 decimals (×100, floor, /100) — this is a cosmetic
re-rounding applied only at the very end before recomputing the reported distance.
If bestSolOverall.IsFeasible==0 → Result = 99999 (no feasible solution found).
Else recompute Result = total distance of bestSolOverall.chromR_move: for each route,
a = [1, (route customers/AFS +1), 1]; Result += Σ dAll(a(k),a(k+1)). Resulttime =
bestSolOverall.time. Return (Result, Resulttime).

----------------------------------------------------------------------
10. INSTANCE LOADING (get_vrp)
----------------------------------------------------------------------
Maps INSTANCE index 1..60 to a .mat file under ./Instances (15_*,25_*,50_*,100_* and
jd200/400/600/800/1000 sets). Each .mat contains the single struct `vrp` with fields in
§1.1. For the C++ port, write a loader that reads these fields (or a converted
text/JSON export). The distance_table is precomputed Euclidean; you may instead rebuild
it from longitude/latitude (dist(i,j)=hypot(dx,dy)) and replicate the AFS clone columns
(nb+1 identical copies of the single AFS column) to reproduce the exact (2nb+2)×(2nb+2)
layout, but simpler is to keep one AFS column and treat any AFS reference as that column.

----------------------------------------------------------------------
11. CONSTANTS, TIE-BREAKS, EDGE CASES (do not lose these)
----------------------------------------------------------------------
- Acceptance threshold everywhere: accept move iff delta < −1e-6 (reject if ≥ −1e-6).
- Best-sol improvement threshold: strictly < current − 1e-9.
- Penalty adaptation band: targetFeasible ± 0.05 (0.15 / 0.25). Scale 1.2, decrease
  0.85, clamp [0.1, 100000].
- wM = 0 permanently → PM = 0; route-count never penalized in this configuration.
- IsFeasible uses UNWEIGHTED penalties (all three must be exactly 0).
- AFS refuel time T_Afs enters ONLY the capacity model PC (via afs_time windows), NOT
  route duration; route duration counts +T_Customer for the AFS visit itself.
- Granular neighbor list = correlatedVertices (5 nearest customers per customer).
  nbGranular(=20) is used ONLY as the reshuffle probability modulus, NOT the list width.
- ELS runs at least two full outer passes (loopID guard). First-improvement, operator
  order m1..m9 then (if V is route start) Depot_m1,2,3,8,9, then (if V is last neighbor
  & loopID≠1) NewRoute_m1,2,3.
- Depot_m2/Depot_m3 reject when target route lacks an AFS (no CRI there).
- NewRoute_m2/m3 reject when the source route has exactly 3 nodes incl AFS.
- Route deletion (isdelete): applies −2·dAll(1,nb+2) (or −dAll(1,nb+2) in one m8 branch)
  to the delta to credit the removed depot↔AFS travel, and removes the route + its AFS
  node, decrementing higher routeIDs and AFS array positions (deleteAFS_node/_delta).
- RNG: rng(SEED+tspid) is set (a) before init split of each individual, (b) inside
  ELS_mian before the granular reshuffle, (c) before the main-loop split. selectparents
  and Crossover use the ambient stream. For bit-reproducibility you'd need MATLAB's
  Mersenne-Twister; for a C++ port, match the STRUCTURE (per-iteration reseeding points)
  but exact numeric reproduction of random draws is not achievable without replicating
  MATLAB's RNG — treat SEED handling as "reseed at these points with a PRNG".

----------------------------------------------------------------------
12. KNOWN DISCREPANCIES / MISSING PIECES (flag in the port)
----------------------------------------------------------------------
(1) T_able is called by split_Tmax but has NO definition in the repo → implement per
    §6.3 (route-duration ≤ Tmax test).
(2) Main loop splits stale `tsp` instead of `offspring_tsp` (§7). Decide faithful vs
    fixed and document.
(3) The display line prints PD using PC's value (cosmetic only; PD=430 is used).
(4) update_feasiblePop/update_infeasiblePop divide broken distance by nbClients and use
    a slightly different node count than the live path (add2Pop uses ÷last_customer).
    The live main-loop path is add2Pop / infeasiblePop_updateBiasedFitnesses; port those.
(5) Route_related col 4 is dual-purposed (afs_time vs a flag "2") inside the ELS patch;
    if you adopt the recommended canonicalize() approach (§5.4 porting guidance) you
    avoid this entirely by recomputing afs_time from scratch each accepted move.

This specification, together with the file list above, is sufficient to implement the
METS C++ port without the MATLAB source. All key source files are under
/Users/matejacivkaroski/Documents/PyCharmProject/METS-Algorithm.