import os
import shutil

def move_files_with_prefix(folder_a, folder_b, prefix):
    # Ensure folder B exists
    if not os.path.exists(folder_b):
        os.makedirs(folder_b)

    # Iterate through all files in folder A
    for file_name in os.listdir(folder_a):
        file_path = os.path.join(folder_a, file_name)

        # Check if it's a file and starts with the specified prefix
        if os.path.isfile(file_path) and file_name.startswith(prefix):
            destination_path = os.path.join(folder_b, file_name)
            shutil.move(file_path, destination_path)
            print(f"Moved: {file_name} to {folder_b}")

# Define the folders and the prefix
folder_a = "DBD"  # Replace with the path to folder A
folder_b = "DBD_ke"  # Replace with the path to folder B
prefix = "DBD_NL_ke"

# Call the function
move_files_with_prefix(folder_a, folder_b, prefix)
