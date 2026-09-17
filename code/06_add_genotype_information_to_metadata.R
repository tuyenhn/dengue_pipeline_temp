# List of required packages
required_packages <- c(
  "optparse", "dplyr", "lubridate", "tidyr",
  "readr", "ape", "seqinr", "stringr", "readr"
)

# Function to check and install packages
# and load necessary utility functions
source("code/_headers.R")

# Define and parse command-line options
opt_parser <- OptionParser(
  option_list = list(
    make_option(c("-m", "--metadata"), type = "character", help = "Input csv file containing cleaned metadata from GenBank, with sequence identifiers matching those in the input fasta file."),
    make_option(c("-b", "--genotype"), type = "character", help = "Input csv file containing typing data, with sequence identifiers matching those in the input metadata file."),
    make_option(c("-l", "--outfile_csv"), type = "character", help = "Outfile csv")
  )
)

opt <- parse_args(opt_parser)

########################################################################
## main
########################################################################


## read in typing file
typing_df <- safe_read_file_param(
  opt[["genotype"]], read_delim,
  delim = ";", show_col_types = FALSE, required = TRUE
)
.read_probs <- problems(typing_df)
if (nrow(.read_probs) > 0) {
  warn_msg("Warnings when reading file")
  .read_probs
}

typing_df <- typing_df %>%
  select(c("seqName", "clade")) %>% # select genotype column
  rename(genotype = clade) # rename to genotype

# read in metadata file
metadata_df <- safe_read_file_param(opt[["metadata"]], read_csv, show_col_types = FALSE, required = TRUE)

metadata_df <- left_join(metadata_df, typing_df, by = c("Sequence_name" = "seqName"))

# give any empty clades a value of unassigned
metadata_df <- metadata_df %>%
  replace_na(list(genotype = "unassigned"))

# write out metadata file
write_csv(metadata_df, file = opt[["outfile_csv"]])
info_msg("Adding serotype information completed")
