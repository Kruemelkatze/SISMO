"""
Test / benchmark harness for lean_code_mp_kw9.process_intersections.

process_intersections is the only piece of the pathfinding pipeline that is
actually callable from outside NetLogo (the A* search itself lives entirely
in a-star.nls and only runs inside the NetLogo world). This harness therefore
targets process_intersections, calling it exactly the way network-creation.nls
does it:

    py:set "raw_data" coordinates-list
    let pyresult py:runresult "process_intersections(raw_data)"

i.e. with a single positional argument and no explicit num_processes /
chunk_size, so the same defaults (num_processes=1, chunk_size=100) apply.

raw_data's shape, confirmed by reading both sides of the py boundary:
  - NetLogo's `coordinates-list` is built from `[pseudopodia] of breeds`,
    which (since `pseudopodia` is the auto-generated global agentset
    reporter for that breed, not a turtle variable) evaluates to the *same*
    global agentset for every member of `breeds`. `item 0 pseudos` therefore
    picks out that whole agentset, and `[path-list] of item 0 pseudos` is a
    list of the path-list of *every currently-alive pseudopodium* in the
    simulation, not just the one that triggered the event.
  - process_intersections then does `for pseudos in raw_data: ... len(pseudos)`
    and indexes `pseudos[reverse_index][2]`, which only works if each element
    of raw_data is itself a path-list of [x, y, marker] triples. This is
    independent confirmation of the point above: raw_data must be
    List[List[[x, y, marker]]], one path-list per pseudopodium.

So raw_data = [ path_list_of_pseudopod_0, path_list_of_pseudopod_1, ... ]
   path_list_of_pseudopod_i = [[x0, y0, m0], [x1, y1, m1], ...]

marker (index 2) is 0 for "real" trail points and non-zero for points
belonging to a hatched copy's duplicated prefix (see replace-zero /
replace-zeros in math_functions.nls). A segment is only considered for
intersection testing if *both* its endpoints have marker == 0.
"""

import math
import random
import time

from lean_code_mp_kw9 import process_intersections


# ---------------------------------------------------------------------------
# Synthetic data generation, mirroring what `wiggle` produces over time.
# ---------------------------------------------------------------------------

def make_pseudopod_path(num_points, origin=(0.0, 0.0), step=1.0, seed=None, marker=0):
    """A random-walk path-list of [x, y, marker] triples, like one
    pseudopodium's `path-list` after `num_points` calls to `wiggle`."""
    rng = random.Random(seed)
    x, y = origin
    heading = rng.uniform(0, 360)
    path = [[x, y, marker]]
    for _ in range(num_points - 1):
        heading += rng.uniform(-40, 40)
        x += math.cos(math.radians(heading)) * step
        y += math.sin(math.radians(heading)) * step
        path.append([x, y, marker])
    return path


def make_raw_data(num_pseudopodia, path_length, seed=None):
    """A list of `num_pseudopodia` random-walk path-lists, each with
    `path_length` points, all marker == 0 (the common case: no hatched
    copies in flight)."""
    return [
        make_pseudopod_path(path_length, origin=(i * 5.0, 0.0), seed=None if seed is None else seed + i)
        for i in range(num_pseudopodia)
    ]


# ---------------------------------------------------------------------------
# Correctness tests
# ---------------------------------------------------------------------------

def test_two_crossing_segments_report_intersection():
    # Pseudopod 0: horizontal segment (0,3) -> (4,3)
    # Pseudopod 1: vertical segment (2,1) -> (2,5)
    # These cross at (2, 3). Coordinates are chosen so no coordinate of one
    # segment numerically equals a coordinate of the other on a different
    # axis -- see test_known_bug_axis_mismatched_prefilter_suppresses_intersection
    # for why that matters.
    raw_data = [
        [[0.0, 3.0, 0], [4.0, 3.0, 0]],
        [[2.0, 1.0, 0], [2.0, 5.0, 0]],
    ]
    result = process_intersections(raw_data, log_context="test")
    hits = [r for r in result if r[0] is not None]
    assert len(hits) == 1, f"expected exactly one intersection, got {hits}"
    x, y = hits[0][0]
    assert abs(x - 2.0) < 1e-9 and abs(y - 3.0) < 1e-9, f"expected intersection at (2,3), got ({x}, {y})"
    print("PASS: test_two_crossing_segments_report_intersection")


def test_known_bug_axis_mismatched_prefilter_suppresses_intersection():
    # process_segment_chunk's pre-filter before calling calculate_intersection is:
    #   seg1[0]!=seg2[0] and seg1[1]!=seg2[1] and seg1[2]!=seg2[2] and
    #   seg1[3]!=seg2[3] and seg1[3]!=seg2[2] and seg1[2]!=seg2[0]
    # seg = [x1, y1, x2, y2]. The last two clauses compare seg1[3] (a
    # y-coordinate) against seg2[2] (an x-coordinate), and seg1[2] (an
    # x-coordinate) against seg2[0] (also an x-coordinate, so that one
    # clause is at least axis-consistent) -- but the seg1[3] vs seg2[2]
    # clause is comparing unlike axes. It appears intended as a cheap
    # "segments don't share an endpoint" pre-check but doesn't actually
    # compare endpoint pairs, so it can misfire whenever a y-coordinate of
    # one segment happens to numerically equal an x-coordinate of the other
    # -- e.g. two segments that both touch coordinate 0 near the NetLogo
    # world's center, which is common in this model. Two segments that
    # genuinely cross (at the origin, same geometry as the test above just
    # shifted) are silently dropped instead of reported.
    raw_data = [
        [[-1.0, 0.0, 0], [1.0, 0.0, 0]],   # horizontal, through the origin
        [[0.0, -1.0, 0], [0.0, 1.0, 0]],   # vertical, through the origin
    ]
    result = process_intersections(raw_data, log_context="test")
    hits = [r for r in result if r[0] is not None]
    # This asserts the CURRENT (buggy) behavior -- zero hits, even though the
    # segments geometrically cross at (0,0). If this assertion ever starts
    # failing, the prefilter bug has been fixed upstream and this test
    # (and its docstring note in the analysis) should be updated/removed.
    assert len(hits) == 0, (
        "expected the known prefilter bug to suppress this intersection; "
        f"got {hits} -- has process_segment_chunk's prefilter been fixed?"
    )
    print("PASS (documents known bug): test_known_bug_axis_mismatched_prefilter_suppresses_intersection")


def test_parallel_segments_report_no_intersection():
    raw_data = [
        [[0.0, 0.0, 0], [1.0, 0.0, 0]],
        [[0.0, 1.0, 0], [1.0, 1.0, 0]],
    ]
    result = process_intersections(raw_data, log_context="test")
    hits = [r for r in result if r[0] is not None]
    assert len(hits) == 0, f"parallel segments should not intersect, got {hits}"
    print("PASS: test_parallel_segments_report_no_intersection")


def test_marker_filters_out_hatched_copy_segments():
    # Same crossing geometry as the first test, but the vertical segment's
    # points are marked as belonging to a hatched copy (marker=1). Per the
    # filter in process_intersections (`pseudos[reverse_index][2] == 0`),
    # the segment should be skipped entirely, so no intersection is found.
    raw_data = [
        [[-1.0, 0.0, 0], [1.0, 0.0, 0]],
        [[0.0, -1.0, 1], [0.0, 1.0, 1]],
    ]
    result = process_intersections(raw_data, log_context="test")
    hits = [r for r in result if r[0] is not None]
    assert len(hits) == 0, f"marker-filtered segment should not be tested, got {hits}"
    print("PASS: test_marker_filters_out_hatched_copy_segments")


def test_single_pseudopod_self_path_no_error():
    # The common case: one pseudopod's own path-list, no other trails yet.
    # (network-creation.nls always passes ALL currently-alive pseudopodia,
    # but a single-pseudopod raw_data is the minimal valid input and is
    # useful as a smoke test.)
    raw_data = [make_pseudopod_path(25, seed=1)]
    result = process_intersections(raw_data, log_context="test")
    assert isinstance(result, list)
    print("PASS: test_single_pseudopod_self_path_no_error")


def run_correctness_tests():
    test_two_crossing_segments_report_intersection()
    test_known_bug_axis_mismatched_prefilter_suppresses_intersection()
    test_parallel_segments_report_no_intersection()
    test_marker_filters_out_hatched_copy_segments()
    test_single_pseudopod_self_path_no_error()


# ---------------------------------------------------------------------------
# Scaling benchmark: reproduces "every step takes longer as pseudopodia
# number/length grow" by sweeping both axes and timing the exact NetLogo
# call shape `process_intersections(raw_data)` (no explicit num_processes /
# chunk_size, so num_processes=1 / chunk_size=100 defaults apply, matching
# production). Each call also appends a row to output/py_intersection_timing.csv
# (log_context="test") via the instrumentation now built into
# lean_code_mp_kw9.process_intersections.
# ---------------------------------------------------------------------------

def run_scaling_benchmark(pseudopodia_counts=(5, 20, 50),
                           path_lengths=(10, 50, 200)):
    # Kept deliberately small by default (~30s total) since compute_time
    # scales roughly O((pseudopodia * path_length)^2) -- see the analysis
    # write-up for a full run out to 100 pseudopodia x 500 points (~13 min,
    # topping out at 768s for a single call at 49900 segments).
    print(f"{'num_pseudopodia':>16} {'path_length':>12} {'segments':>10} {'elapsed_s':>10}")
    for n in pseudopodia_counts:
        for length in path_lengths:
            raw_data = make_raw_data(n, length, seed=42)
            total_segments = n * (length - 1)
            start = time.perf_counter()
            process_intersections(raw_data, log_context="test")
            elapsed = time.perf_counter() - start
            print(f"{n:>16} {length:>12} {total_segments:>10} {elapsed:>10.4f}")


if __name__ == "__main__":
    print("Running correctness tests...")
    run_correctness_tests()
    print()
    print("Running scaling benchmark (mirrors NetLogo call: process_intersections(raw_data))...")
    print("Results are also logged to output/py_intersection_timing.csv (context=test)")
    run_scaling_benchmark()
