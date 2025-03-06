import pandas as pd
import subprocess
import os
import csv

# Define paths relative to the script's directory
script_dir = os.path.dirname(os.path.abspath(__file__))  # Get script's directory
data_dir = os.path.join(script_dir, "..", "data")  # Path to data/
results_dir = os.path.join(script_dir, "..", "results")  # Path to results/

# Ensure results directory exists
os.makedirs(results_dir, exist_ok=True)

# Input and output file paths
tsv_file = os.path.join(data_dir, "nextstrain_fetch_genome.ucsc.edu_trash_hgPhyloPlace_subtreeAuspice1_genome_30f0cb_92e850.json_metadata.tsv")
output_csv = os.path.join(results_dir, "samn_assembly_accessions.csv")

# Load TSV and extract biosample_accession column
df = pd.read_csv(tsv_file, sep="\t")  # Use sep="\t" to read TSV
biosamples = df["biosample_accession"].dropna().unique()

# Function to fetch assembly accessions for a biosample using Entrez Direct
def get_assembly_accessions(biosample):
    """Fetch assembly accessions from NCBI using Entrez Direct."""
    cmd = f'esearch -db assembly -query "{biosample}" | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout.strip().split("\n") if result.stdout else []

# Create an empty list to hold the results
results = []

# Process each biosample and fetch assembly accessions
for biosample in biosamples:
    print(f"Fetching assembly accessions for Biosample: {biosample}")
    accessions = get_assembly_accessions(biosample)
    if accessions:
        for accession in accessions:
            results.append([biosample, accession])
    else:
        results.append([biosample, "No assembly found"])

# Write the results to a CSV file
with open(output_csv, mode='w', newline='') as file:
    writer = csv.writer(file)  # Default delimiter is comma for CSV
    writer.writerow(["Biosample Accessions", "Assembly Accessions"])  # Write header
    writer.writerows(results)  # Write data

print(f"Results have been written to {output_csv}")
