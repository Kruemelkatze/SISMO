# SISMO A* / intersection pipeline — performance analysis

Scope: `lean_code_mp_kw9.py` and `SISMO_V12_FS.nlogo` (+ `network-creation.nls`,
`a-star.nls`), investigating why every simulation step gets slower as the
number and length of pseudopodia grow.

## Root cause

**Critical finding** (`network-creation.nls:5-6`):

```netlogo
let pseudos [pseudopodia] of breeds
let coordinates-list [path-list] of item 0 pseudos
```

`pseudopodia` is not a turtle variable — it's the auto-generated **global
agentset** for that breed (like `turtles`, but scoped to the breed). So
`[pseudopodia] of breeds` evaluates that same global reporter once per member
of `breeds` and returns a list of copies of it; `item 0 pseudos` picks out
that whole agentset. The result: `coordinates-list` — and therefore
`raw_data` sent to Python — is the full `path-list` of **every currently-alive
pseudopodium in the simulation**, not just the one near the depleted food
source.

This is confirmed independently from the Python side: `process_intersections`
does `pseudos[reverse_index][2]`, which only type-checks if each element of
`raw_data` is itself a `[x, y, marker]` path-list — i.e. `raw_data` must be a
list of per-pseudopod paths, matching the NetLogo-side reading above.

Since `path-list` is append-only (`wiggle` only ever grows it, never trims),
this input **grows monotonically for the entire run**, on both axes the user
named: more pseudopodia → more entries; longer trails → more points per
entry.

**Compounding this**, `create-pseudopodia-network` → `process_intersections`
recomputes from scratch, via an O(n²) all-pairs loop (`process_segment_chunk`
— each chunk still scans against the *entire* `result_arr` regardless of
chunk boundaries), **every single time any food source is depleted** — not
incrementally.

## Empirical confirmation

`test_process_intersections.py` calls `process_intersections(raw_data)`
exactly as `network-creation.nls` does (no explicit `num_processes` /
`chunk_size` → defaults `num_processes=1`, `chunk_size=100` apply):

| pseudopodia | path length | segments | elapsed |
|---:|---:|---:|---:|
| 5 | 10 | 45 | 0.07s |
| 5 | 500 | 2,495 | 1.50s |
| 20 | 200 | 3,980 | 3.80s |
| 20 | 500 | 9,980 | 23.45s |
| 50 | 200 | 9,950 | 23.50s |
| 50 | 500 | 24,950 | 146.58s |
| 100 | 200 | 19,900 | 91.65s |
| 100 | 500 | 49,900 | **768.09s** |

Doubling segment count consistently ~quadruples the time (e.g.
9,980→19,900 segments: 23.4s→91.6s, a 3.9x jump for a 2.0x input increase) —
clean O(n²), confirmed. `build_time` (assembling `result_arr` from
`raw_data`) is negligible (microseconds) in every run, logged via
`output/py_intersection_timing.csv` — 100% of the cost is the pairwise
compute, ruling that phase out as a contributor.

There's also a fixed ~70-90ms tax per call regardless of size:
`Pool(num_processes)` spawns and tears down a fresh OS process on every
single invocation even though `num_processes=1` (the default NetLogo always
uses, since it calls `process_intersections(raw_data)` with no extra args)
means the pool provides zero real parallelism.

## Secondary bug found via the test cases (correctness, not performance)

`process_segment_chunk`'s pre-filter, before calling `calculate_intersection`:

```python
if seg1[0] != seg2[0] and seg1[1] != seg2[1] and seg1[2] != seg2[2] and \
   seg1[3] != seg2[3] and seg1[3] != seg2[2] and seg1[2] != seg2[0]:
```

`seg = [x1, y1, x2, y2]`. The clause `seg1[3] != seg2[2]` compares a
y-coordinate against an x-coordinate — mismatched axes, not a real
shared-endpoint check. It silently drops genuine intersections whenever an
unrelated y-coordinate of one segment numerically equals an x-coordinate of
the other (common near the NetLogo world's center, e.g. coordinate 0).
Documented in `test_process_intersections.py` as
`test_known_bug_axis_mismatched_prefilter_suppresses_intersection` — it
currently asserts the *buggy* behavior (zero hits for two segments that
geometrically cross at the origin) so it flags loudly if the prefilter is
ever "fixed" without the test being updated.

## Where A* (`a-star.nls`) fits

KW's instinct that this is "the A* implementation" isn't unreasonable — A*
runs on the same ever-growing graph, and its main loop does two O(graph size)
linear scans per iteration:

- `min-one-of (searchers with [active?]) [total-expected-cost]` — rescans all
  active searchers every iteration instead of using a priority queue.
- `searchers-in-loc` — a linear coordinate-comparison scan per neighbor
  expansion.

But A* operates on a graph *built from* the O(n²) intersection step above,
and per the measurements, the intersection computation looks like the
dominant cost by a wide margin. Rather than resolve this on paper, both
stages are now instrumented — the timing CSVs from a real sim run will show
the actual split.

## Instrumentation added

- **`lean_code_mp_kw9.py`**: `process_intersections` now logs every call to
  `output/py_intersection_timing.csv` (`build_time_s`, `compute_time_s`,
  `total_time_s`, segment/pseudopod counts, `num_processes`, `chunk_size`).
  New `log_context` param (default `"sim"`; the test harness passes `"test"`)
  tags rows by source without changing the NetLogo call signature.
- **`network-creation.nls`**: `create-pseudopodia-network` logs to
  `output/netlogo_timing.csv`, splitting the py roundtrip (marshal + Python
  compute) from the NetLogo build-loop time (the `one-of networkpoints with
  [...]` linear-scan section) via one intermediate `timer` read right after
  `py:runresult`.
- **`a-star.nls`**: `run-a-star` logs its own row (elapsed, graph size, path
  length found) to the same `network-timing-csv`.
- **`SISMO_V12_FS.nlogo`**: new `output/tick_timing.csv`, one row per
  simulation step (elapsed seconds, pseudopod/food counts, avg/max path
  length). Since `create-pseudopodia-network` only fires on food depletion
  (not every tick), this will show whether the "every step" slowdown is
  per-tick (something else growing) or concentrated at depletion events.

## Test cases added

`test_process_intersections.py` (repo root) calls `process_intersections(raw_data)`
exactly as `network-creation.nls` does, with `raw_data` shaped as the
confirmed `List[List[[x, y, marker]]]`. Includes:

- Crossing- and parallel-segment correctness checks.
- Marker-filter behavior (mirrors `replace-zeros` — a segment is only tested
  if both endpoints have `marker == 0`).
- The documented axis-mismatched prefilter bug (see above).
- A parameterized scaling benchmark (defaults trimmed to ~30s; the full
  100×500 sweep in the table above was run once separately and is reported
  here rather than left as the default, since it takes ~13 minutes).

## Proposed fixes (plan only — not implemented)

1. **Incremental network updates** (highest leverage): cache
   `result_arr`/network state across calls instead of recomputing all
   segments from scratch on every food depletion; only test new segments
   against the existing set. Turns O(calls × total²) into ~O(total²) once +
   O(Δ × total) per call.
2. **Real spatial algorithm**: replace the O(n²) all-pairs scan with a
   sweep-line (Bentley–Ottmann — there's already a dead
   `;;py:run "from bentley-ottman import ..."` line in `network-creation.nls`,
   suggesting this was considered before) or a spatial grid/bounding-box
   bucket to skip far-apart segment pairs.
3. **Fix the multiprocessing overhead**: either pass a real `num_processes`
   from NetLogo, or drop `Pool` entirely and call `process_segment_chunk`
   inline — `num_processes=1` currently pays process-spawn cost for zero
   parallel benefit.
4. **Fix the axis-mismatched prefilter bug** (correctness, not perf) — should
   likely compare endpoint pairs `(x, y)` together, not individual
   mismatched-axis components.
5. **Replace linear-scan `networkpoints` lookups** in `network-creation.nls`
   / `a-star.nls` with a `table:` extension keyed by `(xcor, ycor)`, removing
   the O(networkpoints) cost per lookup (several per intersection, plus
   `min-one-of` / `searchers-in-loc` per A* iteration).
6. **Trim/simplify `path-list`** (e.g. Douglas-Peucker or fixed-stride
   subsampling before sending to Python) to slow the growth rate of the input
   itself, independent of the algorithmic fixes above.

Recommended order: run a real sim with the new instrumentation first to get
the actual roundtrip/build/A* split, then do #1 + #3 (cheap, safe, big win)
before considering #2 (bigger rewrite).
