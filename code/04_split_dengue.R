# List of required packages
required_packages <- c("optparse", "dplyr", "lubridate", "tidyr", "readr", "ape", "seqinr", "readr", "purrr", "stringr")

# Function to check and install packages
# and load necessary utility functions
source("code/_headers.R")

# Define and parse command-line options
opt_parser <- OptionParser(
  option_list = list(
    make_option(c("-m", "--metadata"), type = "character", help = "Input csv file containing cleaned metadata, with sequence identifiers matching those in the input fasta file."),
    make_option(c("-f", "--fasta"), type = "character", help = "Input fasta file, with sequence identifiers matching those in the metadata file."),
    make_option(c("-s", "--serotype"), type = "character", help = "Input serotype name"),
    make_option(c("-o", "--outfile"), type = "character", default = "subsampled", help = "Base name for output files. Files will be named as '<outfile>_'serotype'.fasta', '<outfile>__'serotype'_infoTbl.tsv', and '<outfile>__'serotype'_infoTbl.csv'")
  )
)
opt <- parse_args(opt_parser)

########################################################################
## main
########################################################################

## read in input cleaned metadata file
metadata_df <- safe_read_file_param(opt[["metadata"]], read_csv, show_col_types = FALSE, required = TRUE)

# merge different serotype naming schemes
metadata_df <- metadata_df %>%
  mutate(
    Serotype = case_when(
      grepl("dengue_virus_2|Dengue_virus_2|dengue_virus_type_2|Dengue_virus_type_2", Virus_name) ~ "Dengue_2",
      grepl("dengue_virus_1|Dengue_virus_1|dengue_virus_type_1|dengue_virus_type_I|Dengue_virus", Virus_name) ~ "Dengue_1",
      .default = Virus_name
    )
  )

info_msg("Serotype table:")
metadata_df %>%
  group_by(Serotype) %>%
  tally()

## read in input fasta file
seqs <- safe_read_file_param(opt[["fasta"]], read.fasta, required = TRUE)

# Filter serotype to match, rename, and write to separate files
chosen_sero <- opt[["serotype"]]

## Subset the metadata for the current serotype
serotype_metadata <- metadata_df %>%
  filter(Serotype == chosen_sero) %>%
  # reformat Sequence_name, dropping `Virus_name` and use `Serotype` instead
  mutate(
    old_seq_name = Sequence_name,
    Sequence_name = paste(GenBank_ID, Country, State, City, Serotype, Date, Decimal_Date, sep = "|")
  )

# Filter fasta seqs for the current serotype and rename sequence names
fasta_acc_ids <- str_split_i(names(seqs), "\\|", 1)
seqs_in_metadat <- seqs[fasta_acc_ids %in% serotype_metadata$GenBank_ID]

## convert sequences to a tibble format
seqs_as_tibble <- map(seqs_in_metadat, ~ {
  tibble::tibble_row(
    nus = paste(unlist(.x)) %>% list(),
    name = attr(.x, "name"),
    Annot = attr(.x, "Annot")
  )
}) %>% bind_rows()

info_msg("Sequences converted to `tidyverse` style")

## left join with metadata
seqs_joined_metadat <- seqs_as_tibble %>%
  left_join(
    serotype_metadata %>% select(old_seq_name, Sequence_name),
    by = join_by("name" == "old_seq_name")
  )

## rename the sequences based on the metadata
## and convert back to "SeqFastadna" format
renamed_seqs <- seqs_joined_metadat %>%
  pmap(\(nus, name, Annot, Sequence_name) {
    seq_obj <- nus

    attributes(seq_obj) <- list(
      name = Sequence_name,
      Annot = Annot,
      class = "SeqFastadna"
    )

    seq_obj
  }) %>%
  set_names(seqs_joined_metadat$Sequence_name)

info_msg("Sequences renamed and converted back to `seqinr`'s format")


seq_name <- as.data.frame(as.matrix(attributes(renamed_seqs)$names))
remove <- filter(seq_name, grepl("NA|NA|NA", V1, fixed = TRUE))
species.to.remove <- remove$V1
vec.names <- unlist(lapply(strsplit(names(renamed_seqs), ";"), function(x) x[length(x)]))
vec.tokeep <- which(!vec.names %in% species.to.remove)

# Write the sequences to a new FASTA file
write.fasta(
  sequences = renamed_seqs[vec.tokeep], names = names(renamed_seqs[vec.tokeep]),
  file.out = paste0(opt[["outfile"]], chosen_sero, ".fasta")
)

# Prepare and write the metadata information table
out_metadata <- serotype_metadata %>% select(
  GenBank_ID, Country, State, City, Serotype,
  Date, Decimal_Date, Source, Sequence_name
)

write_tsv(
  out_metadata,
  file = paste0(opt[["outfile"]], chosen_sero, "_infoTbl.txt"),
  col_names = TRUE,
  quote = "none"
)

write_csv(
  out_metadata,
  file = paste0(opt[["outfile"]], chosen_sero, "_infoTbl.csv"),
  quote = "none"
)

info_msg("Processing completed for ", chosen_sero)
