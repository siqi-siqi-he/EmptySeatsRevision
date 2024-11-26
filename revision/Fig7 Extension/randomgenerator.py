import random
import os

def generate_random_files():
    # Loop to create files t1.txt to t16.txt
    folder_name = "revision/Fig7 Extension/RandomNumbers"
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    for i in range(100):
        file_name = os.path.join(folder_name, f"randomnumbers_sim{i}.txt")  # Create file name
        with open(file_name, "w") as file:
            # Generate 100 random numbers between 0 and 1
            random_numbers = [random.uniform(0, 1) for _ in range(20)]
            # Write each number to the file on a new line
            file.write("\n".join(map(str, random_numbers)))
            print(f"Created {file_name} with 100 random numbers.")

# Call the function
generate_random_files()