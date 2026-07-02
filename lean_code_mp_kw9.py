import time
import os
import csv
from multiprocessing import Pool

eps = 1e-10

TIMING_LOG_PATH = os.path.join("output", "py_intersection_timing.csv")
_TIMING_HEADER = ["timestamp", "context", "num_pseudopodia", "total_segments",
                   "build_time_s", "compute_time_s", "total_time_s",
                   "num_processes", "chunk_size"]


def _log_timing(context, num_pseudopodia, total_segments, build_time, compute_time,
                 total_time, num_processes, chunk_size):
    os.makedirs("output", exist_ok=True)
    write_header = not os.path.exists(TIMING_LOG_PATH)
    with open(TIMING_LOG_PATH, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(_TIMING_HEADER)
        writer.writerow([
            time.strftime("%Y-%m-%d %H:%M:%S"), context, num_pseudopodia, total_segments,
            f"{build_time:.6f}", f"{compute_time:.6f}", f"{total_time:.6f}",
            num_processes, chunk_size,
        ])

def calculate_intersection(A, B):
    x1, y1, x2, y2 = A
    x3, y3, x4, y4 = B

    denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(denominator) < eps:  # Use tolerance for floating point comparison
        return None  # Parallele Linien

    t_numerator = (x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)
    u_numerator = (x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)
    t = t_numerator / denominator
    u = u_numerator / denominator

    if 0 <= t <= 1 and 0 <= u <= 1:
        x = x1 + t * (x2 - x1)
        y = y1 + t * (y2 - y1)
        # Check if the intersection point is not on the endpoints at the same time
        # Use tolerance for floating point comparison
        is_on_endpoint_A_start = abs(x - x1) < eps and abs(y - y1) < eps
        is_on_endpoint_A_end = abs(x - x2) < eps and abs(y - y2) < eps
        is_on_endpoint_B_start = abs(x - x3) < eps and abs(y - y3) < eps
        is_on_endpoint_B_end = abs(x - x4) < eps and abs(y - y4) < eps
        
        # Nur filtern, wenn der Schnittpunkt genau auf beiden Segment-Endpunkten gleichzeitig liegt
        # (d.h. wenn beide Segmente den gleichen Endpunkt teilen)
        if not ((is_on_endpoint_A_start and is_on_endpoint_B_start) or
                (is_on_endpoint_A_start and is_on_endpoint_B_end) or
                (is_on_endpoint_A_end and is_on_endpoint_B_start) or
                (is_on_endpoint_A_end and is_on_endpoint_B_end)):
            return [x, y]
    return None


def process_segment_chunk(args):
    #print("Processing chunk %d" % count)
    start_i, end_i, result_arr = args
    #print("Index overview: start: {}, end: {}".format(start_i, end_i))
    intersections = []
    #  if (x-1 != x-1-compare) and (y-1 != y-2-compare) and (x-2 != x-2-compare) and (y-2 != y-2-compare) and (y-2 != y-1-compare) and (x-2 != x-1-compare)
    count_of_no_calcs = 0
    for i in range(start_i, end_i):
        added = False
        for j in range(i + 1, len(result_arr)):
            #print ("internal index overview: start: {}, end: {}".format(i,j))
            seg1 = result_arr[i]
            seg2 = result_arr[j]
            #if seg1[0] != seg2[0] and seg1[1] != seg2[1] and seg1[2] != seg2[2] and seg1[3] != seg2[3] and seg1[3] != seg2[2] and seg1[2] != seg2[0]:
            if seg1[0] != seg2[0] and seg1[1] != seg2[1] and seg1[2] != seg2[2] and seg1[3] != seg2[3] and seg1[3] != seg2[2] and seg1[2] != seg2[0]:
                temp = calculate_intersection(seg1, seg2)
                #print("temp: {}".format(temp))
                if temp is not None:
                    intersections.append([
                        temp,
                        [seg1[0], seg1[1]],
                        [seg1[2], seg1[3]],
                        [seg2[0], seg2[1]],
                        [seg2[2], seg2[3]]
                    ])
                    added = False #--100625
                elif added is False:
                    intersections.append([
                        temp,
                        [seg1[0], seg1[1]],
                        [seg1[2], seg1[3]]
                    ])
                    added = True #--100625
            else:
                count_of_no_calcs += 1
    #print("Count of no calculations:", count_of_no_calcs)
    return [end_i,intersections]

def remove_duplicate_intersections(intersections):
    eq_counter = 0
    l_begin = len(intersections)
    index = 0
    while index < l_begin:
        comp_index = index+1
        while comp_index<l_begin:

            if intersections[comp_index][0] is None and intersections[index][0] is None:

                #print(intersections[comp_index][1])

                if (intersections[comp_index][1][0] == intersections[index][1][0] and \
                        intersections[comp_index][1][1] == intersections[index][1][1] and \
                        intersections[comp_index][2][0] == intersections[index][2][0] and \
                        intersections[comp_index][2][1] == intersections[index][2][1]):

                    del intersections[comp_index]
                    l_begin-=1
                    eq_counter += 1
                    #print("Geht4")
                elif intersections[comp_index][0] != None and intersections[index][0] != None:
                    #print("Geht3")
                    if (intersections[comp_index][0][0] == intersections[index][0][0] and
                            intersections[comp_index][0][1] == intersections[index][0][1]):
                        del intersections[comp_index]
                        l_begin -= 1
            #print("copm index:", comp_index)
            comp_index+=1
    index+=1

   # print("equal shit {} ".format(eq_counter))

    return intersections

def process_intersections(raw_data, num_processes=1, chunk_size=100, log_context="sim"):
    # Erstellung der result_arr wie im ursprünglichen Code
    run_start = time.perf_counter()

    result_arr = []
    #print("{} raw_data length: {}".format(time.strftime("[%H:%M:%S]"), len(raw_data)))

    for pseudos in raw_data:
        coord_list_length = len(pseudos)
        # Process all n-1 segments: (0,1), (1,2), ..., (n-2, n-1)
        # Using reverse indexing with reverse_index and reverse_index-1
        for coord_list_index in range(0, coord_list_length - 1):
            reverse_index = coord_list_length - coord_list_index - 1
            if pseudos[reverse_index][2] == 0 and pseudos[reverse_index - 1][2] == 0:
                x1 = pseudos[reverse_index - 1][0]
                y1 = pseudos[reverse_index - 1][1]
                x2 = pseudos[reverse_index][0]
                y2 = pseudos[reverse_index][1]
                result_arr.append([x1, y1, x2, y2])

    #print("{} raw_data length: {}".format(time.strftime("[%H:%M:%S]"), len(result_arr)))
    build_time = time.perf_counter() - run_start

    # Aufteilung der Arbeit in Chunks
    compute_start = time.perf_counter()

    tasks = []
    total_segments = len(result_arr)
    for chunk_start in range(0, total_segments, chunk_size):
        chunk_end = min(chunk_start + chunk_size, total_segments)
        tasks.append((chunk_start, chunk_end, result_arr))

    # Parallele Verarbeitung
    with Pool(num_processes) as pool:
        results = pool.map(process_segment_chunk, tasks)

    # Sammeln der Ergebnisse
    intersections = []
    for res in results:
        #print("Result collector index : {}".format(res[0]))
        intersections.extend(res[1])

   # print("Len unfiltered: {}".format(len(intersections)))
    #filtered_list = remove_duplicate_intersections(intersections)
    #print("Len filtered: {}".format(len(filtered_list)))
    #sorted_list = sorted(filtered_list, key=lambda x: x[0] is not None)

    compute_time = time.perf_counter() - compute_start
    total_time = time.perf_counter() - run_start

    _log_timing(log_context, len(raw_data), total_segments, build_time, compute_time,
                total_time, num_processes, chunk_size)

    return intersections
    