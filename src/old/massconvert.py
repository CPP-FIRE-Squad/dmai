import os
import sys

# Check if the correct number of command line arguments are provided
if len(sys.argv) != 3:
    print("Usage: python massconvert.py <input_path> <output_path>")
    sys.exit(1)

# Get the input and output paths from the command line arguments
input_path = sys.argv[1]
output_path = sys.argv[2]



# Get a list of files in the input directory
files = os.listdir(input_path)
# Create the output directory if it doesn't exist
os.makedirs(output_path, exist_ok=True)
# Loop through each file and execute gh52gb3.py
for file in files:
    # Check if the file name starts with "snapshot_"
    if file.startswith("snapshot_"):
        # Create the command to call gh52gb3.py with the input and output paths
        command = f"python gh52gb3.py {input_path}/{file} {output_path}/{os.path.splitext(file)[0]}"
        # Execute the command
        print(f"Executing: '{command}'")
        os.system(command)
