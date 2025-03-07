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

# Treat 'Collection_Date' as a string and assign it to 'date' column
filtered_df['date'] = filtered_df['Collection_Date'].astype(str)

# Reformat the date column to a consistent format (YYYY-MM-DD) before exporting to CSV
filtered_df['date'] = pd.to_datetime(filtered_df['date'], errors='coerce')
filtered_df['date'] = filtered_df['date'].dt.strftime('%Y-%m-%d')

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

# Save the filtered and transformed metadata to the output CSV file (comma-separated)
output_metadata_csv = os.path.join(script_dir, "../results/nextstrain_files/metadata_segment_HA_NA.csv")
final_metadata_df.to_csv(output_metadata_csv, index=False)

# Print confirmation for CSV save
print(f"Metadata for Isolates with both Segment 4 and 6 (without redundant rows) has been written to {output_metadata_csv}")

# Define the output path for the TSV file (tab-separated)
output_metadata_tsv = os.path.join(script_dir, "../results/nextstrain_files/metadata_segment_HA_NA.tsv")

# Save the filtered and transformed metadata to the output TSV file (tab-separated)
final_metadata_df.to_csv(output_metadata_tsv, sep='\t', index=False)

# Print confirmation for TSV save
print(f"Metadata for Isolates with both Segment 4 and 6 (without redundant rows) has been written to {output_metadata_tsv}")
