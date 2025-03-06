import pandas as pd
import os

# Get the absolute path to the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Path to the cleaned metadata file (relative to the script's directory)
cleaned_metadata_csv = os.path.join(script_dir, "../results/cleaned_all_metadata_sequences.csv")

# Load the cleaned metadata file
df = pd.read_csv(cleaned_metadata_csv)

# Ensure the necessary columns are available
columns_needed = ["Isolate_Name", "Segment", "Collection_Date", "Country", "Host", "USA", "sequence"]
df = df[columns_needed]

# --- Filtering out Isolate_Names that end with -egg or -cell ---
df = df[~df['Isolate_Name'].str.endswith(('-egg', '-cell'))]

# Dictionary to map state abbreviations to full state names
state_abbreviations = {
    "AL": "Alabama", "AK": "Alaska", "AZ": "Arizona", "AR": "Arkansas", "CA": "California",
    "CO": "Colorado", "CT": "Connecticut", "DE": "Delaware", "FL": "Florida", "GA": "Georgia",
    "HI": "Hawaii", "ID": "Idaho", "IL": "Illinois", "IN": "Indiana", "IA": "Iowa",
    "KS": "Kansas", "KY": "Kentucky", "LA": "Louisiana", "ME": "Maine", "MD": "Maryland",
    "MA": "Massachusetts", "MI": "Michigan", "MN": "Minnesota", "MS": "Mississippi", "MO": "Missouri",
    "MT": "Montana", "NE": "Nebraska", "NV": "Nevada", "NH": "New Hampshire", "NJ": "New Jersey",
    "NM": "New Mexico", "NY": "New York", "NC": "North Carolina", "ND": "North Dakota", "OH": "Ohio",
    "OK": "Oklahoma", "OR": "Oregon", "PA": "Pennsylvania", "RI": "Rhode Island", "SC": "South Carolina",
    "SD": "South Dakota", "TN": "Tennessee", "TX": "Texas", "UT": "Utah", "VT": "Vermont",
    "VA": "Virginia", "WA": "Washington", "WV": "West Virginia", "WI": "Wisconsin", "WY": "Wyoming"
}

# Group by Isolate_Name and filter to include only those that have both Segment 4 and 6
valid_isolates = df.groupby('Isolate_Name')['Segment'].unique()
valid_isolates = valid_isolates[valid_isolates.apply(lambda x: set([4, 6]).issubset(set(x)))].index

# Filter the dataframe to include only rows for Isolate_Names that have both Segment 4 and 6
filtered_df = df[df['Isolate_Name'].isin(valid_isolates)]

# Create the new columns based on the provided logic
filtered_df['strain'] = filtered_df['Isolate_Name']
filtered_df['virus'] = 'avian_flu'  # All rows will have this value
filtered_df['isolate_id'] = filtered_df['Isolate_Name']  # Assuming 'Isolate_Name' is the isolate_id
filtered_df['date'] = filtered_df['Collection_Date']

# --- Assign region and country ---
filtered_df['region'] = filtered_df['Country'].apply(lambda x: 'North America' if x == 'USA' else ('Asia' if x == 'Japan' else 'Other'))

# Set the country column as the actual Country from filtered_df
filtered_df['country'] = filtered_df['Country']

# Replace state abbreviation with full state name in the 'division' column
filtered_df['division'] = filtered_df['USA'].map(state_abbreviations)

# If there are any state abbreviations that don't exist in the mapping (NaN or unknown abbreviations), handle them
filtered_df['division'].fillna('Iwate', inplace=True)

filtered_df['host'] = filtered_df['Host']

# Select the columns in the order needed for the final metadata file
final_columns = ['strain', 'virus', 'isolate_id', 'date', 'region', 'country', 'division', 'host']
final_metadata_df = filtered_df[final_columns]

# Remove redundant rows based on all columns
final_metadata_df = final_metadata_df.drop_duplicates()

# Define the output path for the results (relative to the script's directory)
output_metadata_csv = os.path.join(script_dir, "../results/nextstrain_files/metadata_segment_HA_NA.csv")

# Ensure the nextstrain_files/ directory exists
os.makedirs(os.path.dirname(output_metadata_csv), exist_ok=True)

# Save the filtered and transformed metadata to the output CSV file
final_metadata_df.to_csv(output_metadata_csv, index=False)

print(f"Metadata for Isolates with both Segment 4 and 6 (without redundant rows) has been written to {output_metadata_csv}")

# --- Generating the FASTA files for HA and NA sequences ---

# Filter the dataframe for HA (Segment 4) and NA (Segment 6) sequences
ha_df = filtered_df[filtered_df['Segment'] == 4]
na_df = filtered_df[filtered_df['Segment'] == 6]

# Define the output directory for the FASTA files (relative to the script's directory)
output_fasta_dir = os.path.join(script_dir, "../results/nextstrain_files/")

# Ensure the nextstrain_files/ directory exists
os.makedirs(output_fasta_dir, exist_ok=True)

# Write the HA sequences to a FASTA file
ha_fasta_file = os.path.join(output_fasta_dir, "ha_sequences.fasta")
with open(ha_fasta_file, "w") as ha_file:
    for _, row in ha_df.iterrows():
        header = row["Isolate_Name"]
        sequence = row["sequence"]
        ha_file.write(f">{header}\n{sequence}\n")

# Write the NA sequences to a FASTA file
na_fasta_file = os.path.join(output_fasta_dir, "na_sequences.fasta")
with open(na_fasta_file, "w") as na_file:
    for _, row in na_df.iterrows():
        header = row["Isolate_Name"]
        sequence = row["sequence"]
        na_file.write(f">{header}\n{sequence}\n")

print(f"HA sequences written to {ha_fasta_file}")
print(f"NA sequences written to {na_fasta_file}")
