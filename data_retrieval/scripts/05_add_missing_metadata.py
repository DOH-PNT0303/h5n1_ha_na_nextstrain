import pandas as pd
import subprocess
import os
import csv

# Configuration
script_dir = os.path.dirname(os.path.abspath(__file__))  # Get script directory
tsv_file = os.path.join(script_dir, "../data/nextstrain_fetch_genome.ucsc.edu_trash_hgPhyloPlace_subtreeAuspice1_genome_30f0cb_92e850.json_metadata.tsv")
output_csv = os.path.join(script_dir, "../results/samn_accessions_with_metadata.csv")

# Load TSV and extract relevant columns
df = pd.read_csv(tsv_file, sep="\t")  # Use sep="\t" to read TSV

# Columns you need to extract
metadata_columns = ["biosample_accession", "date", "host", "authors"]

# Ensure only the necessary columns are available
df = df[metadata_columns].dropna(subset=["biosample_accession"])

# Create an empty list to store results
results = []

# Function to fetch nucleotide accessions using Entrez Direct
def get_nucleotide_accessions(biosample):
    """Fetch nucleotide accessions from NCBI using Entrez Direct."""
    cmd = f'esearch -db nucleotide -query "{biosample}" | efetch -format docsum | xtract -pattern DocumentSummary -element Caption'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout.strip().split("\n") if result.stdout else []

# Process each biosample and fetch nucleotide accessions
for _, row in df.iterrows():
    biosample = row["biosample_accession"]
    date = row["date"]
    host = row["host"]
    authors = row["authors"]

    # Fetch nucleotide accessions
    print(f"Fetching accession for Biosample: {biosample}")
    accessions = get_nucleotide_accessions(biosample)

    if accessions:
        for accession in accessions:
            # Add ".1" to each accession
            accession_with_version = f"{accession}.1"
            results.append([biosample, accession_with_version, date, host, authors])
    else:
        results.append([biosample, "No sequence found", date, host, authors])

# Ensure the results directory exists
os.makedirs(os.path.join(script_dir, "../results/"), exist_ok=True)

# Filter out rows where the 'Accession' is "No sequence found"
df_results = pd.DataFrame(results, columns=["BioSample", "Accession", "Collection_Date", "Host", "Submitters"])
df_results = df_results[df_results['Accession'] != "No sequence found"]

# Convert 'Collection_Date' to datetime, ensuring it is in the correct format (YYYY-MM-DD)
df_results['Collection_Date'] = pd.to_datetime(df_results['Collection_Date'], errors='coerce').dt.strftime('%Y-%m-%d')

# Drop rows with invalid dates (which will be NaT after conversion)
df_results = df_results.dropna(subset=['Collection_Date'])

# Write results to CSV file
df_results.to_csv(output_csv, index=False)

print(f"Results have been written to {output_csv}")

# -------- Merging and Filling Missing Values Logic -------- #

# Absolute path to the existing metadata
existing_metadata_csv = os.path.join(script_dir, "../results/all_merged_metadata_with_sequences.csv")

# Check if the file exists before attempting to read it
if os.path.exists(existing_metadata_csv):
    existing_metadata = pd.read_csv(existing_metadata_csv)
else:
    print(f"File not found: {existing_metadata_csv}")
    # Handle the case where the file is missing (e.g., skip processing or exit the script)
    exit()

# Read the new metadata table (samn_accessions_with_metadata.csv)
new_metadata = pd.read_csv(output_csv)  # The file we've just created

# Strip any leading/trailing whitespace from column names
existing_metadata.columns = existing_metadata.columns.str.strip()
new_metadata.columns = new_metadata.columns.str.strip()

# Ensure 'Collection_Date' in existing metadata is in datetime format (YYYY-MM-DD) before merging
existing_metadata['Collection_Date'] = pd.to_datetime(existing_metadata['Collection_Date'], errors='coerce').dt.strftime('%Y-%m-%d')

# Drop rows with invalid dates in the existing metadata (which will be NaT after conversion)
existing_metadata = existing_metadata.dropna(subset=['Collection_Date'])

# Perform a left join to merge the tables on 'Accession'
merged_df = pd.merge(existing_metadata, new_metadata, how="left", on="Accession")

# Dynamically handle filling missing values for Host, Collection_Date, and Authors
columns_to_merge = ['Host', 'Collection_Date', 'Submitters', 'BioSample']

for column in columns_to_merge:
    # Check for '_x' and '_y' columns in the merged dataframe
    if f'{column}_x' in merged_df.columns and f'{column}_y' in merged_df.columns:
        merged_df[column] = merged_df[f'{column}_x'].fillna(merged_df[f'{column}_y'])
        # Drop the '_x' and '_y' columns after merging
        merged_df = merged_df.drop(columns=[f'{column}_x', f'{column}_y'])

# Save the merged data to a new CSV file
output_merged_csv = os.path.join(script_dir, "../results/cleaned_all_metadata_sequences_no_outgroup.csv")  # You can adjust the output path here
merged_df.to_csv(output_merged_csv, index=False)

print(f"Merged metadata has been written to {output_merged_csv}")
