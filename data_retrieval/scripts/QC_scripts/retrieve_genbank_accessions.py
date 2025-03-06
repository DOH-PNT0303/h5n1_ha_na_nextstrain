import pandas as pd
import subprocess
import os
import glob
from concurrent.futures import ThreadPoolExecutor

# Define paths relative to the script's directory
script_dir = os.path.dirname(os.path.abspath(__file__))  # Get script's directory
data_dir = os.path.join(script_dir, "..", "data")  # Path to data/
results_dir = os.path.join(script_dir, "..", "results")  # Path to results/
nucleotides_dir = os.path.join(results_dir, "nucleotides_from_NCBI")  # Path for nucleotide FASTA files
final_fasta = os.path.join(results_dir, "all_sequences.fna")  # Output final FASTA
metadata_output_file = os.path.join(results_dir, "metadata_table.csv")  # Output metadata table

# Ensure output directories exist
os.makedirs(results_dir, exist_ok=True)
os.makedirs(nucleotides_dir, exist_ok=True)

# Input TSV file path
tsv_file = os.path.join(data_dir, "nextstrain_fetch_genome.ucsc.edu_trash_hgPhyloPlace_subtreeAuspice1_genome_30f0cb_92e850.json_metadata.tsv")

# Load TSV and extract biosample_accession column
df = pd.read_csv(tsv_file, delimiter='\t')  # Use \t for tab-separated values
biosamples = df["biosample_accession"].dropna().unique()

# Function to fetch nucleotide accessions using Entrez Direct
def get_nucleotide_accessions(biosample):
    """Fetch nucleotide accessions from NCBI using Entrez Direct."""
    cmd = f'esearch -db nucleotide -query "{biosample}" | efetch -format docsum | xtract -pattern DocumentSummary -element Caption'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout.strip().split("\n") if result.stdout else []

# Set to hold nucleotide accessions to avoid duplicates
nucleotide_accessions = set()

# Fetch nucleotide accessions for each biosample
for biosample in biosamples:
    print(f"Fetching nucleotide accessions for Biosample: {biosample}")
    accessions = get_nucleotide_accessions(biosample)
    if accessions:
        nucleotide_accessions.update(accessions)
    else:
        print(f"No nucleotide accession found for {biosample}")

# Function to fetch metadata for each biosample (date, host, author)
def get_metadata_for_biosample(biosample):
    """Fetch metadata (date, host, author) for a given biosample from the TSV file."""
    metadata_row = df[df["biosample_accession"] == biosample]
    if not metadata_row.empty:
        date = metadata_row["date"].values[0]
        host = metadata_row["host"].values[0]
        author = metadata_row["author"].values[0]
        return biosample, date, host, author
    return biosample, None, None, None  # Return None for missing values

# Function to download and extract FASTA files for a given nucleotide accession using Entrez Direct
def download_and_extract_nucleotide(nucleotide):
    fasta_path = os.path.join(nucleotides_dir, f"{nucleotide}.fasta")

    # Skip if the file already exists
    if os.path.exists(fasta_path):
        print(f"Skipping download: {nucleotide} already exists.")
    else:
        print(f"Downloading: {nucleotide}")
        # Fetch the nucleotide sequence using efetch
        cmd = f'esearch -db nucleotide -query "{nucleotide}" | efetch -format fasta > {fasta_path}'
        subprocess.run(cmd, shell=True)

    return fasta_path

# Parallelize downloading and extracting using ThreadPoolExecutor
with ThreadPoolExecutor() as executor:
    fasta_files = list(executor.map(download_and_extract_nucleotide, nucleotide_accessions))

# Remove None values from the list (in case any FASTA files weren't found)
fasta_files = [file for file in fasta_files if file]

# Merge all FASTA files into one
if fasta_files:
    print("Merging all FASTA files into one...")
    with open(final_fasta, "w") as outfile:
        for fname in fasta_files:
            with open(fname) as infile:
                outfile.write(infile.read())

    print(f"Final FASTA file saved as: {final_fasta}")
else:
    print("No FASTA files were downloaded.")

# Fetch metadata for all biosample_accession and create metadata table
metadata_list = []
for biosample in biosamples:
    biosample_metadata = get_metadata_for_biosample(biosample)
    metadata_list.append(biosample_metadata)

# Create the metadata DataFrame
metadata_df = pd.DataFrame(metadata_list, columns=["biosample_accession", "date", "host", "author"])

# Merge metadata with the nucleotide accessions
metadata_df["Accession"] = metadata_df["biosample_accession"].apply(lambda biosample: get_nucleotide_accessions(biosample)[0] if get_nucleotide_accessions(biosample) else None)

# Save the metadata table to a CSV file
metadata_df.to_csv(metadata_output_file, index=False)

print(f"Metadata table saved to: {metadata_output_file}")
print("Process complete!")
