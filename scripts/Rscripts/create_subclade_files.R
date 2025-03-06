library(tidyverse)
library(readxl)

raw <- read_excel("data/metadata.xlsx")
metadata <- read_excel("data/metadata.xlsx") %>% 
  rename(strain = Isolate_Name,
         isolate_id = Isolate_Id, 
         date = Collection_Date, 
         host = Host) %>% 
  mutate(Location = trimws(Location)) %>% # Trim whitespace around the Location variable
  separate(Location, into = c("region", "country", "division", "location"), sep = " / ", extra = "merge")
  
subclade <- read.table(file = 'data/nextstrain_groups_wadoh.tsv', sep = '\t', header = TRUE)

subclade_meta <- metadata  %>% 
  right_join(subclade) %>% 
  left_join(metadata) %>% 
  bind_rows(
    metadata %>% 
      filter(strain %in% c("A/Washington/254/2024-egg", 
                                 "A/chicken/Iwate/22A5T/2023")) %>% 
      distinct(strain, .keep_all = TRUE)
  )

full_genomes <- subclade_meta %>% 
  filter(
    !is.na(`PB2 Segment_Id`) &
      !is.na(`PB1 Segment_Id`) &
      !is.na(`PA Segment_Id`) &
      !is.na(`HA Segment_Id`) &
      !is.na(`NP Segment_Id`) &
      !is.na(`NA Segment_Id`) &
      !is.na(`MP Segment_Id`) &
      !is.na(`NS Segment_Id`)
  ) %>% 
  mutate(virus = "avian_flu") %>% 
  select(strain, virus, isolate_id, date, region, country, division, host)

fragmented_genomes <- subclade_meta %>% 
  filter(!(strain %in% full_genomes$strain))

inferred_genomes <- subclade_meta %>% 
  mutate(virus = "avian_flu") %>% 
  select(strain, virus, isolate_id, date, region, country, division, host) %>% 
  filter(!strain %in% c("A/Washington/240/2024", 
                              "A/Washington/254/2024", 
                              "A/Wisconsin/179/2024", 
                              "A/Iowa/124/2024"))

# Function to read FASTA and extract relevant information
read_fasta_to_df <- function(fasta_file) {
  # Read the FASTA file
  fasta_data <- readLines(fasta_file)
  
  # Identify header lines
  header_indices <- grep("^>", fasta_data)
  
  # Extract headers
  headers <- fasta_data[header_indices]
  
  # Extract sequences by grouping lines between headers
  sequences <- sapply(seq_along(header_indices), function(i) {
    start <- header_indices[i] + 1
    end <- if (i < length(header_indices)) header_indices[i + 1] - 1 else length(fasta_data)
    paste(fasta_data[start:end], collapse = "")  # Combine sequence lines
  })
  
  # Create a data frame with headers and sequences
  parsed_data <- tibble(
    Header = sub("^>", "", headers),  # Remove the leading ">" from headers
    Sequence = sequences  # Assign combined sequences
  ) %>%
    separate(
      Header, 
      into = c("Isolate_Id", "Isolate_Name", "Segment"), 
      sep = "\\|",  # Split using "|" delimiter
      remove = TRUE  # Remove the original Header column
    )
  
  return(parsed_data)
}

############
# Read in genomes that have all segments 
############
fasta_file <- "data/Franklin_County.fasta"
fasta_df <- read_fasta_to_df(fasta_file) %>%
  select(Isolate_Name, Sequence, Segment) %>%
  mutate(
    Isolate_Name = if_else(
      Isolate_Name == "A/Washington/254/2024", 
      "A/Washington/254/2024-egg", 
      Isolate_Name)) %>% 
  mutate(
    Sequence = toupper(Sequence),
    Segment = tolower(Segment)
  )

######
# Read in fragmented genomes 
######
fasta_file2 <- "data/Franklin_County_all.fasta"
# Read and preprocess the fasta file
fasta_df2 <- read_fasta_to_df(fasta_file2) %>%
  mutate(
    Isolate_Name = if_else(
      Isolate_Name == "A/Washington/254/2024" & Isolate_Id == "EPI_ISL_19666173", 
      "A/Washington/254/2024-egg", 
      Isolate_Name),
    Sequence = toupper(Sequence),
    Segment = tolower(Segment)
  ) %>% 
  select(Isolate_Name, Sequence, Segment)

# Define expected segments and their corresponding sequence lengths
segment_lengths <- c(pb1 = 2341, pb2 = 2341, pa = 2233, ha = 1760, 
                     np = 1565, na = 1458, mp = 1027, ns = 865)
segment_lengths <- c(pb1 = 100, pb2 = 100, pa = 100, ha = 100, 
                     np = 100, na = 100, mp = 100, ns = 100)

# Ensure every Isolate_Name has all segments
fasta_df_complete <- fasta_df2 %>%
  complete(Isolate_Name, Segment = names(segment_lengths)) %>%
  mutate(
    Sequence = if_else(
      is.na(Sequence),
      strrep("N", segment_lengths[Segment]),
      Sequence
    )
  ) %>% 
  filter(!Isolate_Name %in% c("A/Washington/240/2024", "A/Washington/254/2024"))

write_fasta_by_segment <- function(data, output_dir = ".") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  data %>%
    group_by(Segment) %>%
    group_walk(~ {
      # Create FASTA content
      fasta_content <- paste0(">", .x$Isolate_Name, "\n", .x$Sequence, collapse = "\n")
      
      # Generate the file name following the pattern
      file_name <- file.path(output_dir, paste0("sequences_h5n1_", .y$Segment, ".fasta"))
      
      # Write the content to the file
      writeLines(fasta_content, con = file_name)
      
      message("Written: ", file_name)
    })
}

# 1. Read the `.fas` file
fas_file <- "outbreak_files_franklin/Franklin03.fas"  # Replace with your file path

# Read the FASTA file
sequences <- readLines(fas_file)

# Extract the headers (lines that start with '>')
headers <- sequences[str_detect(sequences, "^>")]

# 2. Extract strain names from the headers
strain_names <- str_remove(headers, "^>")  # Remove '>' from the headers to get strain names

# 3. Define the metadata to add for all sequences
metadata <- data.frame(
  strain = strain_names,
  virus = rep("avian_flu", length(strain_names)),
  region = rep("North America", length(strain_names)),
  country = rep("United States", length(strain_names)),
  division = rep("Washington", length(strain_names)),
  host = rep("Chicken", length(strain_names)),
  date = NA  # We will set the date based on the strain name
)

# 4. Assign dates based on strain names
metadata$date <- ifelse(str_detect(strain_names, "_2024$"), "2024-10-11", 
                        ifelse(str_detect(strain_names, "_2025$"), "2025-XX-XX", NA))

# 5. Merge the new metadata with your existing full_genomes dataframe
# Assuming 'full_genomes' already has a 'strain' column for matching

full_genomes <- full_genomes %>%
  full_join(metadata)


# Call the function
write_fasta_by_segment(fasta_df, output_dir = "outbreak_files_franklin")
write.table(full_genomes, file = "outbreak_files_franklin/sequences_h5n1_franklin_county.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(full_genomes$strain, file = "outbreak_files_franklin/include_strains_outbreak.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

## add in sequences from Franklin03.fas 



write_fasta_by_segment(fasta_df_complete, output_dir = "output_fasta")
write.table(inferred_genomes, file = "output_fasta/sequences_h5n1_franklin_county.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(inferred_genomes$strain, file = "output_fasta/include_strains_outbreak.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 
# 
# fasta_df2 <- read_fasta_to_df(fasta_file2) %>%
#   mutate(
#     Isolate_Name = if_else(
#       Isolate_Name == "A/Washington/254/2024" & Isolate_Id == "EPI_ISL_19666173", 
#       "A/Washington/254/2024-egg", 
#       Isolate_Name),
#     Sequence = toupper(Sequence),
#     Segment = tolower(Segment)
#   ) %>% 
#   select(Isolate_Name, Sequence, Segment) %>% 
#   filter(!Isolate_Name %in% c("A/Washington/240/2024"))
# write_fasta_by_segment(fasta_df2, output_dir = "output_fasta")
# write.table(inferred_genomes, file = "output_fasta/sequences_h5n1_franklin_county.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
# write.table(inferred_genomes$strain, file = "output_fasta/include_strains_outbreak.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
