import os
import pandas as pd
from Bio import SeqIO
import re

# Define file paths relative to the root directory
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))  # Get scripts/ directory
ROOT_DIR = os.path.dirname(ROOT_DIR)  # Move up to project root

FASTA_FILE = os.path.join(ROOT_DIR, "results/all_sequences.fna")
CSV_FILE = os.path.join(ROOT_DIR, "results/merged_metadata_sequences_mod.csv")
PARSED_FASTA_FILE = os.path.join(ROOT_DIR, "results/parsed_fasta.csv")
OUTPUT_FILE = os.path.join(ROOT_DIR, "results/all_merged_metadata_with_sequences.csv")

def clean_whitespace_and_linebreaks(df, column_name):
    """Remove any line breaks, tabs, and trim leading/trailing whitespace in the specified column."""
    if column_name in df.columns:
        # Remove line breaks (both \n and \r) and other whitespace characters
        df[column_name] = df[column_name].str.replace(r'[\n\r\t]+', '', regex=True)  # Remove newlines, tabs
        df[column_name] = df[column_name].str.strip()  # Trim leading/trailing whitespace
    return df

# Function to extract Isolate_Name from GenBank_Title
def extract_isolate_name(genbank_title):
    match = re.search(r"\((.*?)\((H5N1)\)\)", genbank_title)
    return match.group(1).strip() if match else None

def extract_segment_number(genbank_title):
    """Extract segment number from GenBank_Title."""
    segment_match = re.search(r'segment (\d+)', genbank_title)
    if segment_match:
        return segment_match.group(1)
    return ""

def extract_state(isolate_name):
    """Extract state from Isolate_Name."""
    parts = isolate_name.split('/')
    if len(parts) >= 4:  # Assuming the state is at position 3 (index 2)
        return parts[2]  # The third part is the state
    return None  # Return None if state cannot be extracted

def parse_fasta(fasta_path, output_path):
    """Parse FASTA file and extract Accession, GenBank_Title, and sequence, then save to CSV."""
    fasta_data = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        header_parts = record.description.split(" ", 1)  # Split into Accession and GenBank_Title
        accession = header_parts[0]  # First part is Accession
        genbank_title = header_parts[1] if len(header_parts) > 1 else ""  # Everything else is GenBank_Title
        sequence = str(record.seq)  # Get sequence

        isolate_name = extract_isolate_name(genbank_title)
        segment_number = extract_segment_number(genbank_title)
        state = extract_state(isolate_name)

        # Prepare data for the FASTA dataframe
        fasta_data.append({
            "Accession": accession,
            "GenBank_Title": genbank_title,
            "sequence": sequence,
            "Isolate_Name": isolate_name,
            "Segment": segment_number if segment_number else "",
            "USA": state,
            "Country": "USA"  # Default country value is "USA"
        })

    fasta_df = pd.DataFrame(fasta_data)

    # Clean whitespace and line breaks from GenBank_Title and sequence columns
    fasta_df = clean_whitespace_and_linebreaks(fasta_df, 'GenBank_Title')
    fasta_df = clean_whitespace_and_linebreaks(fasta_df, 'sequence')
    fasta_df = clean_whitespace_and_linebreaks(fasta_df, 'Isolate_Name')

    # Fill missing values in the 'USA' column with the extracted state from Isolate_Name
    fasta_df['USA'] = fasta_df['USA'].fillna(fasta_df['Isolate_Name'].apply(extract_state))

    # Fill missing values in the 'Country' column with "USA"
    fasta_df['Country'] = fasta_df['Country'].fillna("USA")

    fasta_df.to_csv(output_path, index=False)
    print(f"Parsed and cleaned FASTA data saved to {output_path}")
    return fasta_df

def merge_metadata(fasta_df, metadata_path, output_path):
    """Merge parsed FASTA data with metadata CSV, keeping entries in metadata if already present, and merge the rest."""
    if not os.path.exists(metadata_path):
        print(f"Error: Metadata CSV file not found at {metadata_path}")
        return

    metadata_df = pd.read_csv(metadata_path)

    # Clean whitespace and line breaks from GenBank_Title and sequence columns in the metadata dataframe
    metadata_df = clean_whitespace_and_linebreaks(metadata_df, 'GenBank_Title')
    metadata_df = clean_whitespace_and_linebreaks(metadata_df, 'sequence')

    # Find the Accessions that are already in the metadata
    existing_accessions = metadata_df['Accession'].unique()

    # Filter out those rows from parsed FASTA data that are already in metadata
    new_fasta_df = fasta_df[~fasta_df['Accession'].isin(existing_accessions)]

    # Merge the metadata with the filtered parsed FASTA data (excluding existing accessions)
    merged_df = pd.concat([metadata_df, new_fasta_df], ignore_index=True)

    # Save merged data
    merged_df.to_csv(output_path, index=False)
    print(f"Merging complete. Output saved to {output_path}")

def main():
    # Check if FASTA file exists
    if not os.path.exists(FASTA_FILE):
        print(f"Error: FASTA file not found at {FASTA_FILE}")
        return

    # Parse FASTA and save to CSV
    fasta_df = parse_fasta(FASTA_FILE, PARSED_FASTA_FILE)

    # Merge parsed FASTA with metadata and save the final output
    merge_metadata(fasta_df, CSV_FILE, OUTPUT_FILE)

if __name__ == "__main__":
    main()
