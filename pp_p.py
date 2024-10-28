import time
from multiprocessing import Pool

def compute_sum_of_squares(start_end):
    start, end = start_end
    result = 0
    for i in range(start, end):
        result += i * i
    return [[result,0], [2*result, 1,3]]

if __name__ == "__main__":
    start_time = time.time()

    # Define the range of numbers
    start, end = 0, 10_000_000
    num_processes = 22  # Adjust this to the number of cores available

    # Create chunks of the range to distribute across processes
    chunk_size = (end - start) // num_processes
    ranges = [(i * chunk_size, (i + 1) * chunk_size) for i in range(num_processes)]
    print(ranges)
    # Use a Pool to parallelize the computation
    with Pool(processes=num_processes) as pool:
        results = pool.map(compute_sum_of_squares, ranges)
    print(results)
    print(results[1])
    print(results[1][1])
    # Sum up the results from all processes
    #total_result = sum(results)

    #print(f"Result: {total_result}")
    print(f"Time taken with parallel processing: {time.time() - start_time:.2f} seconds")
