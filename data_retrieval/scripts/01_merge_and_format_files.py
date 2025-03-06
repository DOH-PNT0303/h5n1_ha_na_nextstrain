import os
from Bio import SeqIO
import pandas as pd
import re

# Define paths relative to the script's location
script_dir = os.path.dirname(os.path.abspath(__file__))  # Get the directory of the script
data_dir = os.path.join(script_dir, "..", "data")  # Navigate up one level to access the data folder
results_dir = os.path.join(script_dir, "..", "results")  # Navigate up one level to access the results folder

# Ensure results directory exists
os.makedirs(results_dir, exist_ok=True)

# Step 1: Read the FASTA file and create a dictionary of sequences
fasta_file = os.path.join(data_dir, "2024_09_01_2024_12_31_WA_OR_NCBI_sequences.fasta")

# Function to extract accession from FASTA header
def extract_accession(header):
    match = re.match(r"([^|]+)", header)
    return match.group(1).strip() if match else None

# Function to extract Isolate_Name from GenBank_Title
def extract_isolate_name(genbank_title):
    match = re.search(r"\((.*?)\((H5N1)\)\)", genbank_title)
    return match.group(1).strip() if match else None

# Create a dictionary of sequences
fasta_dict = {extract_accession(record.id): str(record.seq)
              for record in SeqIO.parse(fasta_file, "fasta")
              if extract_accession(record.id)}

# Step 2: Read the metadata CSV file
metadata_file = os.path.join(data_dir, "2024_09_01_2024_12_31_WA_OR_NCBI_metadata.csv")
metadata_df = pd.read_csv(metadata_file)

# Step 3: Merge sequences with metadata
metadata_df["sequence"] = metadata_df["Accession"].map(fasta_dict)

# Step 4: Extract 'Isolate_Name' from 'GenBank_Title'
metadata_df["Isolate_Name"] = metadata_df["GenBank_Title"].apply(extract_isolate_name)

# Step 5: Save the merged data
merged_metadata_file = os.path.join(results_dir, "merged_metadata_sequences.csv")
metadata_df.to_csv(merged_metadata_file, index=False)

# Confirmation message
print(f"Results saved to: {merged_metadata_file}")
