import time

def compute_sum_of_squares(start, end):
    result = 0
    for i in range(start, end):
        result += i * i
    return result

if __name__ == "__main__":
    start_time = time.time()

    # Define the range of numbers
    start, end = 0, 10_000_000

    # Compute the sum of squares without parallel processing
    result = compute_sum_of_squares(start, end)

    print(f"Result: {result}")
    print(f"Time taken without parallel processing: {time.time() - start_time:.2f} seconds")