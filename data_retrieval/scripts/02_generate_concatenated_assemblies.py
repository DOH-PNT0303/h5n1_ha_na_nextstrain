import os
import pandas as pd

# Define paths relative to the script's location
script_dir = os.path.dirname(os.path.abspath(__file__))  # Get the directory of the script
results_dir = os.path.join(script_dir, "..", "results")  # Navigate up one level to access the results folder

# Ensure results directory exists
os.makedirs(results_dir, exist_ok=True)

# Input and output file paths
metadata_file = os.path.join(results_dir, "merged_metadata_sequences_mod.csv")
output_fasta_file = os.path.join(results_dir, "concatenated_assemblies.fasta")

# Load metadata file
metadata_df = pd.read_csv(metadata_file)

# Ensure necessary columns are present
required_columns = {"Isolate_Name", "Segment", "sequence"}
if not required_columns.issubset(metadata_df.columns):
    raise ValueError(f"Missing required columns: {required_columns - set(metadata_df.columns)}")

# Step 1: Group by Isolate_Name and filter for those with all 8 segments
isolate_df = metadata_df.groupby("Isolate_Name").filter(lambda x: len(x["Segment"].unique()) == 8)

# Step 2: Concatenate sequences for Isolates with all 8 segments
concatenated_sequences = {
    isolate: "".join(group.sort_values("Segment")["sequence"])
    for isolate, group in isolate_df.groupby("Isolate_Name")
}

# Step 3: Write concatenated sequences to FASTA file
with open(output_fasta_file, "w") as fasta_out:
    for isolate_name, seq in concatenated_sequences.items():
        fasta_out.write(f">{isolate_name}\n{seq}\n")

# Confirmation message
print(f"Concatenated sequences saved to: {output_fasta_file}")

