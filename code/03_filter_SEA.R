# List of required packages
required_packages <- c("dplyr", "lubridate", "tidyr", "optparse", "ape", "seqinr", "readr")

# Function to check and install packages
# and load necessary utility functions
source("code/_headers.R")

# Define and parse command-line options
opt_parser <- OptionParser(
  option_list = list(
    make_option(c("-e", "--metadata"), type = "character", help = "Input tsv file containing metadata from GenBank, with sequence identifiers matching those in the input fasta file."),
    make_option(c("-m", "--metadata_vietnam"), type = "character", help = "Input tsv file containing metadata from GenBank, with sequence identifiers matching those in the input fasta file."),
    make_option(c("-g", "--metadata_china"), type = "character", help = "Input tsv file containing metadata from GenBank, with sequence identifiers matching those in the input fasta file."),
    make_option(c("-x", "--metadata_extra"), type = "character", help = "Input tsv file containing extra metadata, with sequence identifiers matching those in the input fasta file."),
    make_option(c("-f", "--fasta"), type = "character", help = "Input fasta file, with sequence identifiers matching those in the metadata file."),
    make_option(c("-o", "--outfile"), type = "character", default = "subsampled", help = "Base name for output files. Files will be named as '<outfile>.fasta'.")
  )
)

opt <- parse_args(opt_parser)

########################################################################
## main
########################################################################

## read in input metadata file
metadata_df <- safe_read_file_param(opt[["metadata"]], read_csv, show_col_types = FALSE, required = TRUE)

## read in input Vietnam metadata file
metadata_vietnam <- safe_read_file_param(opt[["metadata_vietnam"]], read_csv, show_col_types = FALSE, required = TRUE)

## read in input China metadata file
metadata_china <- safe_read_file_param(opt[["metadata_china"]], read_csv, show_col_types = FALSE, required = TRUE)

## read in input extra metadata file
metadata_extra <- safe_read_file_param(opt[["metadata_extra"]], read_tsv, show_col_types = FALSE, required = TRUE)

## read in input fasta file
seqs <- safe_read_file_param(opt[["fasta"]], read.fasta, required = TRUE)

########################################################################

metadata_vietnam <- select(metadata_vietnam, c("Accession", "Region"))
metadata_china <- select(metadata_china, c("Accession", "Region"))

metadata_vietnam <- metadata_vietnam %>%
  mutate(Region = case_when(
    grepl("Central|Central|Central (maybe north)", Region, ignore.case = TRUE) ~ "Central_Vietnam",
    grepl("North|North (maybe central)", Region, ignore.case = TRUE) ~ "Northern_Vietnam",
    grepl("South", Region, ignore.case = TRUE) ~ "Southern_Vietnam"
  ))

# Remove NA
metadata_vietnam <- metadata_vietnam %>% filter(!is.na(Region))

metadata_china <- metadata_china %>%
  mutate(Region = case_when(
    grepl("Fujian", Region, ignore.case = TRUE) ~ "Fujian",
    grepl("Guangdong", Region, ignore.case = TRUE) ~ "Guangdong",
    grepl("Hainan", Region, ignore.case = TRUE) ~ "Hainan",
    grepl("Henan", Region, ignore.case = TRUE) ~ "Henan",
    grepl("Hunan", Region, ignore.case = TRUE) ~ "Hunan",
    grepl("Jiangxi", Region, ignore.case = TRUE) ~ "Jiangxi",
    grepl("Shandong", Region, ignore.case = TRUE) ~ "Shandong",
    grepl("Shanghai", Region, ignore.case = TRUE) ~ "Shanghai",
    grepl("Sichuan", Region, ignore.case = TRUE) ~ "Sichuan",
    grepl("Yunnan", Region, ignore.case = TRUE) ~ "Yunnan",
    grepl("Zhejiang", Region, ignore.case = TRUE) ~ "Zhejiang"
  )) %>%
  filter(!is.na(Region))

# always include sequences from extra metadata
joint_metadata <- bind_rows(metadata_vietnam, metadata_china, metadata_extra)

# Filter metadata

SEA <- c(
  "Brunei", "Burma", "Myanmar", "Cambodia", "Indonesia", "Laos", "Malaysia",
  "Philippines", "Singapore", "Thailand", "India", "Taiwan"
)

metadata_vietnam_china <- filter(metadata_df, GenBank_ID %in% joint_metadata$Accession)
metadata_sea <- filter(metadata_df, Country %in% SEA)

metadata_final <- bind_rows(metadata_vietnam_china, metadata_sea)
info_msg("Metadata filtering completed")

# filter sequences

seq_name <- as.data.frame(as.matrix(attributes(seqs)$names))

keep <- seq_name %>%
  filter(V1 %in% metadata_final$Sequence_name)

taxa_split <- data.frame(do.call("rbind", strsplit(as.character(keep$V1), "|", fixed = TRUE)))
taxa_split$name <- keep$V1
species.to.keep <- taxa_split$name

vec.names <- unlist(lapply(strsplit(names(seqs), ";"), function(x) x[length(x)]))

vec.tokeep <- which(vec.names %in% species.to.keep)

# write files

write.fasta(
  sequences = seqs[vec.tokeep], names = names(seqs)[vec.tokeep],
  file.out = paste0(opt$outfile, ".fasta")
)

write.csv(metadata_final,
  file = paste0(opt$outfile, "_infoTbl.tsv"),
  row.names = FALSE, quote = FALSE
)

write.csv(metadata_final,
  file = paste0(opt$outfile, "_infoTbl.csv"),
  row.names = FALSE, quote = FALSE
)

info_msg("Sequence filtering completed")
