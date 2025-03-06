import os
import re
from Bio import SeqIO
import pandas as pd

# Define paths relative to the script's location
script_dir = os.path.dirname(os.path.abspath(__file__))  # Get the directory of the script
data_dir = os.path.join(script_dir, "..", "data")  # Navigate up one level to access the data folder
results_dir = os.path.join(script_dir, "..", "results")  # Navigate up one level to access the results folder

# Ensure results directory exists
os.makedirs(results_dir, exist_ok=True)

# Function to extract accession from FASTA header
def extract_accession(header):
    match = re.match(r"([^|]+)", header)
    return match.group(1).strip() if match else None

# Function to extract Isolate_Name from GenBank_Title
def extract_isolate_name(genbank_title):
    match = re.search(r"\((.*?)\((H5N1)\)\)", genbank_title)
    return match.group(1).strip() if match else None

# Step 1: Read the outgroup FASTA file and create a dictionary of sequences
outgroup_fasta_file = os.path.join(data_dir, "outgroup_sequences.fasta")  # Modify the filename if needed
fasta_dict = {extract_accession(record.id): str(record.seq)
              for record in SeqIO.parse(outgroup_fasta_file, "fasta")
              if extract_accession(record.id)}

# Step 2: Read the metadata CSV file for outgroup
outgroup_metadata_file = os.path.join(data_dir, "outgroup_metadata.csv")  # Modify the filename if needed
outgroup_metadata_df = pd.read_csv(outgroup_metadata_file)

# Step 3: Merge sequences with metadata for outgroup
outgroup_metadata_df["sequence"] = outgroup_metadata_df["Accession"].map(fasta_dict)

# Step 4: Extract 'Isolate_Name' from 'GenBank_Title'
outgroup_metadata_df["Isolate_Name"] = outgroup_metadata_df["GenBank_Title"].apply(extract_isolate_name)

# Step 5: Save the merged outgroup data
merged_outgroup_metadata_file = os.path.join(results_dir, "merged_outgroup_metadata_sequences.csv")
outgroup_metadata_df.to_csv(merged_outgroup_metadata_file, index=False)

# Confirmation message
print(f"Outgroup results saved to: {merged_outgroup_metadata_file}")

# -------- Concatenate Outgroup with Merged Metadata -------- #

# Read the previous merged metadata (without outgroup)
output_merged_csv = os.path.join(script_dir, "../results/cleaned_all_metadata_sequences_no_outgroup.csv")  # Adjust the file path as needed

# Check if the file exists
if os.path.exists(output_merged_csv):
    merged_metadata_df = pd.read_csv(output_merged_csv)
else:
    print(f"File not found: {output_merged_csv}")
    exit()  # Exit if the file doesn't exist

# Read the merged outgroup data
outgroup_df = pd.read_csv(merged_outgroup_metadata_file)

# Concatenate the merged metadata with the outgroup data
final_merged_df = pd.concat([merged_metadata_df, outgroup_df], ignore_index=True)

# Step 6: Save the final merged data
final_output_csv = os.path.join(results_dir, "cleaned_all_metadata_sequences.csv")
final_merged_df.to_csv(final_output_csv, index=False)

# Confirmation message
print(f"Final merged metadata (with outgroup) saved to: {final_output_csv}")
