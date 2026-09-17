# List of required packages
required_packages <- c(
  "optparse", "dplyr", "lubridate", "tidyr", "anytime", "purrr",
  "readr", "ape", "seqinr", "stringr", "magrittr"
)

# Function to check and install packages
# and load necessary utility functions
source("code/_headers.R")

# Define and parse command-line options
opt_parser <- OptionParser(
  option_list = list(
    make_option(c("-m", "--metadata"), type = "character", help = "Input tsv file containing metadata from GenBank, with sequence identifiers matching those in the input fasta file."),
    make_option(c("-f", "--fasta"), type = "character", help = "Input fasta file, with sequence identifiers matching those in the metadata file."),
    make_option(c("-x", "--extra_metadata"), type = "character", help = "Option to add a csv file containing metadata from sequencing but not on Genbank"),
    make_option(c("-z", "--extra_fasta"), type = "character", help = "Option to add a fasta file containing metadata from sequencing but not on Genbank"),
    make_option(c("-j", "--outfile_tsv"), type = "character", help = "Outfile tsv"),
    make_option(c("-k", "--outfile_fasta"), type = "character", help = "Outfile fasta"),
    make_option(c("-l", "--outfile_csv"), type = "character", help = "Outfile csv"),
    make_option(c("-s", "--start_date"), type = "character", default = "2000-01-01", help = "Start date for filtering (format YYYY-MM-DD)"),
    make_option(c("-e", "--end_date"), type = "character", default = "2024-05-31", help = "End date for filtering (format YYYY-MM-DD)"),
    make_option(c("-H", "--host"), type = "character", default = "Homo sapiens", help = "Host type sample for sequencing was taken from")
  )
)

opt <- parse_args(opt_parser)

########################################################################
## main
########################################################################

## read in input metadata file from GenBank
metadata_df <- safe_read_file_param(opt$metadata, read_tsv, show_col_types = FALSE, required = TRUE) %>%
  select(-`Geographic State`) %>%
  mutate(source = "GenBank") %>%
  filter(Accession != "OP909734.1")

## read in and combine extra metadata from sequencing
## Nb will need to be in a specific format to match GenBank metadata
## See Github for specific formatting requirements
extra_metadata <- safe_read_file_param(opt$extra_metadata, read_tsv, show_col_types = FALSE)
problems(extra_metadata)

extra_metadata <- extra_metadata %>% mutate(source = "Extra")

if (nrow(extra_metadata) == 0 || ncol(extra_metadata) == 0) {
  warn_msg("Extra metadata NOT included")
} else {
  info_msg("`metadata_df` columns: ")
  print(colnames(metadata_df))

  info_msg("`extra_metadata` columns: ")
  print(colnames(extra_metadata))
}

# Process dates
metadata_df <- metadata_df %>%
  mutate(
    Date = case_when(
      str_detect(`Isolate Collection date`, "^\\d{4}$") ~ sprintf("%s-06-15", `Isolate Collection date`),
      str_detect(`Isolate Collection date`, "^\\d{4}-\\d{2}$") ~ sprintf("%s-15", `Isolate Collection date`),
      .default = `Isolate Collection date`
    ),
    Date = as.Date(anytime(`Isolate Collection date`))
  )

extra_metadata <- extra_metadata %>% mutate(
  Date = as.Date(`Isolate Collection date`, format = "%d-%m-%Y")
)

metadata_df <- bind_rows(metadata_df, extra_metadata)
info_msg("Extra metadata included")

# Filter host and date
metadata_df <- metadata_df %>%
  filter(`Host Name` == opt[["host"]]) %>%
  filter(between(Date, as.Date(opt[["start_date"]]), as.Date(opt[["end_date"]])))

# Extract state level information
metadata_df <- metadata_df %>%
  separate_wider_delim(
    cols = `Geographic Location`,
    delim = ":",
    names = c("Country", "State"),
    too_many = "error",
    too_few = "align_start"
  ) %>%
  separate_wider_delim(
    cols = State,
    delim = ",",
    names = c("State", "City"),
    too_many = "merge",
    too_few = "align_start"
  ) %>%
  mutate(
    City = str_split_i(City, ",", 1)
  ) %>%
  mutate(across(
    .cols = c(Country, State, City),
    .fns = ~ trimws(.x, which = "both")
  ))


# Select desired columns and remove any with NA in Date and Country
metadata_df <- metadata_df %>%
  select(Accession, `Virus Name`, Date, Country, State, City, source) %>%
  filter(!is.na(Date) & !is.na(Country))

# Match metadata to fasta file and rename sequences
## read in input fasta file
seqs <- safe_read_file_param(opt[["fasta"]], read.fasta, required = TRUE)

extra_seqs <- safe_read_file_param(opt[["extra_fasta"]], read.fasta)
if (is.null(extra_seqs)) {
  warn_msg("Extra fasta sequences are NOT included")
} else {
  full_seqs <- c(seqs, extra_seqs)
  info_msg("Extra fasta sequences are included")
}

seq_names <- names(full_seqs)

## subsetting metadata to sequences in the fasta
## and build new sequence names
metadata_df <- metadata_df %>%
  filter(Accession %in% seq_names) %>%
  select(Accession, Date, `Virus Name`, Country, State, City, source) %>%
  mutate(
    Date = as.character(Date),
    decDate = map(Date, ~ calcDecimalDate_fromTxt(.x, dayFirst = FALSE, namedMonths = FALSE, sep = "-")) %>% unlist()
  ) %>%
  mutate(
    across(
      c(Country, State, City, `Virus Name`),
      ~ str_replace_all(.x, " ", "_") %>% str_replace_all("[()]", "")
    )
  ) %>%
  mutate(Sequence_name = paste(Accession, Country, State, City, `Virus Name`, Date, decDate, sep = "|"))

info_msg("Metadata df subsetted and new sequence names built")

## subsetting the fasta sequences to the metadata
seqs_have_metadat <- full_seqs[seq_names %in% metadata_df$Accession]

info_msg("Fasta sequences subsetted")

## convert sequences to a tibble format
seqs_as_tibble <- map(seqs_have_metadat, ~ {
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
    metadata_df %>% select(Accession, Sequence_name),
    by = join_by("name" == "Accession")
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
  file.out = opt[["outfile_fasta"]]
)

# Prepare and write the metadata information table
out_metadata <- metadata_df %>% rename(
  GenBank_ID = Accession,
  Country = Country,
  State = State,
  City = City,
  Virus_name = `Virus Name`,
  Date = Date,
  Decimal_Date = decDate,
  Source = source,
  Sequence_name = Sequence_name
)

write_tsv(
  out_metadata,
  file = opt[["outfile_tsv"]],
  col_names = TRUE,
  quote = "none"
)

write_csv(
  out_metadata,
  file = opt[["outfile_csv"]],
  quote = "none"
)

info_msg("Processing metadata and sequences completed")
