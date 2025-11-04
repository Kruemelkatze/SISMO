import time
from multiprocessing import Pool

#start = time.perf_counter()
#my_function()
#end = time.perf_counter()
#print(f"Time taken: {end - start:.4f} seconds")

def calculate_intersection(A, B):
    x1, y1, x2, y2 = A
    x3, y3, x4, y4 = B

    denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if denominator == 0:
        return None  # Parallele Linien

    t_numerator = (x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)
    u_numerator = (x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)
    t = t_numerator / denominator
    u = u_numerator / denominator

    if 0 <= t <= 1 and 0 <= u <= 1:
        x = x1 + t * (x2 - x1)
        y = y1 + t * (y2 - y1)
        # Überprüfen, ob der Schnittpunkt nicht auf den Endpunkten liegt
        if (x != x1 and x != x2 and x != x3 and x != x4 and
                y != y1 and y != y2 and y != y3 and y != y4):
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
            print("copm index:", comp_index)
            comp_index+=1
    index+=1

    print("equal shit {} ".format(eq_counter))

    return intersections

def process_intersections(raw_data, num_processes=1, chunk_size=100):
    # Erstellung der result_arr wie im ursprünglichen Code
    start = time.perf_counter()

    result_arr = []
    #print("{} raw_data length: {}".format(time.strftime("[%H:%M:%S]"), len(raw_data)))

    for pseudos in raw_data:
        coord_list_length = len(pseudos)
        for coord_list_index in range(1, coord_list_length - 1):
            reverse_index = coord_list_length - coord_list_index
            if pseudos[reverse_index][2] == 0 and pseudos[reverse_index - 1][2] == 0:
                x1 = pseudos[reverse_index - 1][0]
                y1 = pseudos[reverse_index - 1][1]
                x2 = pseudos[reverse_index][0]
                y2 = pseudos[reverse_index][1]
                result_arr.append([x1, y1, x2, y2])

    print("{} raw_data length: {}".format(time.strftime("[%H:%M:%S]"), len(result_arr)))
    end = time.perf_counter()
    print(f"ITime taken: {end - start:.4f} seconds")

    # Aufteilung der Arbeit in Chunks
    scnd_start = time.perf_counter()

    tasks = []
    total_segments = len(result_arr)
    for start in range(0, total_segments, chunk_size):
        end = min(start + chunk_size, total_segments)
        tasks.append((start, end, result_arr))

    # Parallele Verarbeitung
    count = 0
    with Pool(num_processes) as pool:
        count+=1
        results = pool.map(process_segment_chunk, tasks)


    # Sammeln der Ergebnisse
    intersections = []
    for res in results:
        print("Result collector index : {}".format(res[0]))
        intersections.extend(res[1])
    
    print("Len unfiltered: {}".format(len(intersections)))
    #filtered_list = remove_duplicate_intersections(intersections)
    #print("Len filtered: {}".format(len(filtered_list)))
    #sorted_list = sorted(filtered_list, key=lambda x: x[0] is not None)

    scnd_end = time.perf_counter()

    print(f"IITime taken: {scnd_end - scnd_start:.4f} seconds")
    print(intersections)
    return intersections
    