import pandas as pd
import subprocess
import os
import glob
from concurrent.futures import ThreadPoolExecutor

# Define paths relative to the script's directory
script_dir = os.path.dirname(os.path.abspath(__file__))  # Get script's directory
data_dir = os.path.join(script_dir, "..", "data")  # Path to data/
results_dir = os.path.join(script_dir, "..", "results")  # Path to results/
assemblies_dir = os.path.join(results_dir, "assemblies_from_NCBI")  # Path for assemblies
final_fasta = os.path.join(results_dir, "all_sequences.fna")  # Output final FASTA

# Ensure output directories exist
os.makedirs(results_dir, exist_ok=True)
os.makedirs(assemblies_dir, exist_ok=True)

# Input TSV file path
tsv_file = os.path.join(data_dir, "nextstrain_fetch_genome.ucsc.edu_trash_hgPhyloPlace_subtreeAuspice1_genome_30f0cb_92e850.json_metadata.tsv")

# Load TSV and extract biosample_accession column
df = pd.read_csv(tsv_file, delimiter='\t')  # Use \t for tab-separated values
biosamples = df["biosample_accession"].dropna().unique()

def get_assembly_accessions(biosample):
    """Fetch assembly accessions from NCBI using Entrez Direct."""
    cmd = f'esearch -db assembly -query "{biosample}" | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout.strip().split("\n") if result.stdout else []

assembly_accessions = set()  # Use a set to avoid duplicates

# Fetch assembly accessions for each biosample
for biosample in biosamples:
    print(f"Fetching assembly accessions for Biosample: {biosample}")
    accessions = get_assembly_accessions(biosample)
    if accessions:
        assembly_accessions.update(accessions)
    else:
        print(f"No assembly found for {biosample}")

# Function to download and extract FASTA files for a given assembly
def download_and_extract_assembly(assembly):
    zip_file = os.path.join(assemblies_dir, f"{assembly}.zip")
    fasta_path = os.path.join(assemblies_dir, f"{assembly}")

    # Skip if the file already exists
    if os.path.exists(zip_file):
        print(f"Skipping download: {assembly} already exists.")
    else:
        print(f"Downloading: {assembly}")
        subprocess.run(["datasets", "download", "genome", "accession", assembly, "--include", "genome", "--filename", zip_file])

    # Extract FASTA file
    subprocess.run(["unzip", "-o", zip_file, "-d", fasta_path])

    # Locate the extracted FASTA file
    fasta_pattern = os.path.join(fasta_path, "ncbi_dataset", "data", assembly, "*.fna")
    fasta_files_found = glob.glob(fasta_pattern)

    if fasta_files_found:
        return fasta_files_found[0]  # Return the first found FASTA file
    else:
        print(f"FASTA file not found for {assembly}. Check extracted contents manually.")
        return None

# Parallelize downloading and extracting using ThreadPoolExecutor
with ThreadPoolExecutor() as executor:
    fasta_files = list(executor.map(download_and_extract_assembly, assembly_accessions))

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

print("Process complete!")
