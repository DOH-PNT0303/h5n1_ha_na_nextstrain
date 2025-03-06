# H5N* Oregon and Washington Phylogenetic Tree Construction

## Goal
The objective is to create a concatenated phylogenetic tree for the HA (hemagglutinin) and NA (neuraminidase) segments of recent H5N* sequences from Oregon and Washington. This will focus on sequences available only from NCBI due to data-sharing policies from GISAID.

## How to run the scripts
1. Clone the Repo and change directories to the `data_retrieval` subfolder
```
git clone https://github.com/DOH-PNT0303/h5n1_ha_na_nextstrain.git
cd data_retrieval
```
2. You can execute each script using Python individually to see each step outlined below.
```
python scripts/01_merge_and_format_files.py
python scripts/02_generate_concatenated_assemblies.py
python scripts/03_fetch_influenza_seqs.py
python scripts/04_parse_and_merge_fasta_with_metadata.py
python scripts/05_add_missing_metadata.py
python scripts/06_add_outgroup.py
python scripts/07_create_fasta_metadata_files.py
```

## Explanation of Steps

### 1. Pull H5N* Sequences from NCBI

We start by pulling the H5N* sequences from NCBI and curating them down to the sequences of interest from Oregon and Washington that we have already identified through our HA tree.

#### 1.1 Query NCBI Virus

- Use NCBI Virus to pull all H5N1 sequences and segments with the following filters:
  - **Geographic Region**: Oregon and Washington
  - **Genotype**: H5N* or H5
  - **Collection Date**: From 09/1/2024 to 12/31/2024

[Link to NCBI Virus query](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Alphainfluenzavirus,%20taxid:197911&Serotype_s=H5N*%20H5&USAState_s=OR&USAState_s=WA&CollectionDate_dr=2024-09-01T00:00:00.00Z%20TO%202024-12-31T23:59:59.00Z)

- **Downloaded Data**:
  - Sequences: `data/2024_09_01_2024_12_31_WA_OR_NCBI_sequences.fasta`
  - Metadata: `data/2024_09_01_2024_12_31_WA_OR_NCBI_metadata.csv`

#### 1.2 Merge and Format Sequences

- Use the script `scripts/01_merge_and_format_files.py` to merge and format the downloaded sequences.

#### 1.3 Subset and Clean Data

- Manually subset the merged metadata file `results/merged_metadata_sequences_mod.csv` to keep only the Oregon and Washington sequences of interest based on the subclade of interest in our HA tree. The sequences of interest are listed in `data/washington_oregon_sequences_of_interest.tsv`.

### 2. Pull Contextual Sequences from the D1.1 Clade

To get relevant contextual data we need to ensure that we are using contextual sequences that are the same genotype (D1.1) to minimize the impact of reassortment on phylogenetic signals. To pull contextual sequences from the D1.1 clade we perform phylogenetic placement of the Washington and Oregon sequences of interest that have whole genomes onto the UShER concatenated whole genome tree (N=50 sequences displayed) and pull all sequences from the D1.1 subclade that our Washington sequences cluster within. We do that through the following steps:

#### 2.1 Generate Concatenated Assemblies

- Use the script `scripts/02_generate_concatenated_assemblies.py` to generate the concatenated assemblies in `results/concatenated_assemblies.fasta`.

#### 2.2 Perform Phylogenetic Placement

- Upload the whole genome concatenated sequences (from `results/concatenated_assemblies.fasta`) into [UShER](https://genome.ucsc.edu/cgi-bin/hgPhyloPlace) for placement on the H5N1 D1.1 2024 outbreak concatenated segments tree.

#### 2.3 Export Subtree Metadata

- Export the metadata from the subtree (subtree 1) that contains both the Washington (homo sapiens) cluster and the Oregon chicken sequences. This metadata contains NCBI sequences within the D1.1 subclade, which will provide context for our sequences of interest.

- **Exported Metadata**: `data/nextstrain_fetch_genome.ucsc.edu_trash_hgPhyloPlace_subtreeAuspice1_genome_30f0cb_92e850.json_metadata.txt`

#### 2.4 Fetch Contextual Sequences from NCBI

- Use the script `scripts/03_fetch_influenza_seqs.py` to programmatically pull all contextual sequences from NCBI using the biosample accession from the metadata.

  **Note**: There are 32 sequences with a biosample accession, and of these, 16 have assemblies available in NCBI. The remaining 12 sequences are based on SRA reads and are excluded from the analysis as we do not have access to their assemblies.

- **Downloaded Assemblies**: `results/assemblies_from_NCBI/`
- **Merged Fasta File**: `results/all_sequences.fna`

### 3. Merge and Clean Data

- Merge the contextual sequences with the sequences of interest using the script `scripts/04_parse_and_merge_fasta_with_metadata.py`.
- Deduplicate sequences and fill in any missing metadata using the script `scripts/05_add_missing_metadata.py`.
- Add in the outgroup sequence `A/jungle crow/Iwate/0304I001/2022` using the script `scripts/06_add_outgroup.py`.
- The final cleaned dataset will be saved as `results/cleaned_all_metadata_sequences.csv`.

### 4. (Optional) Prepare Sequences for Nextstrain Build

- Take the cleaned dataset (`results/cleaned_all_metadata_sequences.csv`) and extract the sequences that have both HA and NA segments.
- Write the HA and NA sequences into separate FASTA files, along with their metadata, for use in the Nextstrain phylogenetic tree building pipeline using the script `scripts/07_create_fasta_metadata_files.py`

## Conclusion

Following these steps, we curate the data, generate the required contextual sequences, and prepare the dataset for constructing a concatenated HA-NA phylogenetic tree for recent H5N* sequences from Oregon and Washington.
