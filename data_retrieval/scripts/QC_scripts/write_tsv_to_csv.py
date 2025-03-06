import pandas as pd
import os

# Define paths
script_dir = os.path.dirname(os.path.abspath(__file__))  # Get script directory
tsv_file = os.path.join(script_dir, "../data/nextstrain_fetch_genome.ucsc.edu_trash_hgPhyloPlace_subtreeAuspice1_genome_192df2_918c00.json_metadata.tsv")
csv_file = os.path.join(script_dir, "../results/converted_metadata.csv")

# Load TSV file
df = pd.read_csv(tsv_file, sep="\t")

# Ensure results directory exists
os.makedirs(os.path.join(script_dir, "../results/"), exist_ok=True)

# Write to CSV
df.to_csv(csv_file, index=False)

print(f"TSV file successfully converted to CSV: {csv_file}")

