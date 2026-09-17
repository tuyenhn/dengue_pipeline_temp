## load packages
required_packages <- c(
  "dplyr", "tidyverse", "tidyr", "dplyr",
  "seqinr", "ape", "optparse", "readr"
)

# Function to check and install packages
# and load necessary utility functions
source("code/_headers.R")

option_list <- list(
  make_option(c("-f", "--fasta"), type = "character", help = "Input fasta file."),
  make_option(c("-t", "--tree"), type = "character", help = "Input tree file."),
  make_option(c("-m", "--metadata_china"), type = "character", help = "Input metadata file."),
  make_option(c("-v", "--metadata_vietnam"), type = "character", help = "Input metadata file."),
  make_option(c("-o", "--outfile"), type = "character", default = "subsampled", help = "Base name for output files."),
  make_option(c("-l", "--DTA"), type = "character", default = "location", help = "Name of the column to treat as location.")
)

opt <- parse_args(OptionParser(option_list = option_list))

##########################################################
## main
##########################################################

# Load FASTA
fasta <- safe_read_file_param(opt[["fasta"]], read.fasta, required = TRUE)

seq_name <- data.frame(name = sapply(fasta, function(x) attr(x, "name")))

# Load Tree
tree <- safe_read_file_param(opt[["tree"]], read.nexus, required = TRUE)

tree_name <- data.frame(name = tree$tip.label)

# Load granular Vietnam locations
denv_vietnam <- safe_read_file_param(opt[["metadata_vietnam"]], read_csv, show_col_types = FALSE, required = TRUE)

process_date <- function(x) {
  date <- ifelse(nchar(x$"Collection_Date") == 4, paste(x$"Collection_Date", "06-15", sep = "-"),
    ifelse(nchar(x$"Collection_Date") == 7, paste(x$"Collection_Date", 15, sep = "-"), x$"Collection_Date")
  )
  x$Date <- as.Date(parse_date_time(date, orders = c("mdy", "dmy", "myd", "y", "my", "m", "ymd", "ym")))

  x <- dplyr::select(x, -c("Collection_Date"))
}

metadata_df_vietnam <- process_date(denv_vietnam)

# select desired columns and remove any with NA in Date and Country

metadata_df_vietnam <- metadata_df_vietnam %>%
  dplyr::select(Accession, Organism_Name, Region, Date)

# merge different serotype naming schemes
metadata_df_vietnam <- metadata_df_vietnam %>%
  mutate(serotype = case_when(
    grepl("dengue virus 2|Dengue virus 2|dengue virus type 2", Organism_Name, ignore.case = TRUE) ~ "Dengue_2",
    grepl("dengue virus 3|Dengue virus 3|dengue virus type 3", Organism_Name, ignore.case = TRUE) ~ "Dengue_3",
    grepl("dengue virus 4|Dengue virus 4|dengue virus type 4", Organism_Name, ignore.case = TRUE) ~ "Dengue_4",
    grepl("dengue virus 1|Dengue virus 1|dengue virus type 1|dengue virus type I|Dengue virus", Organism_Name, ignore.case = TRUE) ~ "Dengue_1"
  ))

metadata_df_complete_vietnam <- metadata_df_vietnam %>%
  filter(!Region %in% c("Imported case", "Unknown"))

metadata_df_complete_vietnam <- metadata_df_complete_vietnam %>% mutate(Region = case_when(
  grepl("Central|Central (maybe north)", Region, ignore.case = TRUE) ~ "Central_Vietnam",
  grepl("North (maybe central)|North", Region, ignore.case = TRUE) ~ "North_Vietnam",
  grepl("South", Region, ignore.case = TRUE) ~ "Southern_Vietnam"
))

# load granular china locations
denv_china <- safe_read_file_param(opt[["metadata_china"]], read_csv, show_col_types = FALSE, required = TRUE)

metadata_df_china <- process_date(denv_china)

# select desired columns and remove any with NA in Date and Country
metadata_df_china <- metadata_df_china %>%
  dplyr::select(Accession, Organism_Name, Region, Date)

# merge different serotype naming schemes
metadata_df_china <- metadata_df_china %>%
  mutate(
    serotype = case_when(
      grepl("dengue virus 2|Dengue virus 2|dengue virus type 2", Organism_Name, ignore.case = TRUE) ~ "Dengue_2",
      grepl("dengue virus 3|Dengue virus 3|dengue virus type 3", Organism_Name, ignore.case = TRUE) ~ "Dengue_3",
      grepl("dengue virus 4|Dengue virus 4|dengue virus type 4", Organism_Name, ignore.case = TRUE) ~ "Dengue_4",
      grepl("dengue virus 1|Dengue virus 1|dengue virus type 1|dengue virus type I|Dengue virus", Organism_Name, ignore.case = TRUE) ~ "Dengue_1"
    )
  )

table(metadata_df_china$Region)

metadata_df_complete_china <- metadata_df_china %>%
  filter(!Region %in% c("Imported case", "Unknown"))

metadata_df_complete_china$Region <- iconv(metadata_df_complete_china$Region, to = "UTF-8")

# Replace "\xcaHubei" with "Hubei" in the Region column
metadata_df_complete_china$Region <- gsub("\\\\xcaHubei", "Hubei", metadata_df_complete_china$Region)

# Replace "Yunnan\\xca" with "Yunnan" in the Region column
metadata_df_complete_china$Region <- gsub("Yunnan\\\\xca", "Yunnan", metadata_df_complete_china$Region)

# merge the two dataframes

metadata_df_complete <- rbind(metadata_df_complete_vietnam, metadata_df_complete_china)

metadata_df_complete <- dplyr::select(metadata_df_complete, c("Accession", "Region"))

# Filter FASTA based on tree tips
to_keep <- seq_name$name %in% tree_name$name
filtered_fasta <- fasta[to_keep]

# Metadata for DTA

fasta_split <- data.frame(do.call("rbind", strsplit(as.character(names(filtered_fasta)), "|", fixed = TRUE)))

fasta_split <- left_join(fasta_split, metadata_df_complete, by = c("X1" = "Accession"))

# if NA in region, place with value from X2

fasta_split$Region[is.na(fasta_split$Region)] <- fasta_split$X2[is.na(fasta_split$Region)]

# new dataframe for DTA metadata

dta_metadata <- data.frame(traits = names(filtered_fasta), location = fasta_split$Region)

# write table
write.table(dta_metadata, file = paste0(opt[["outfile"]], "_metadata.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
# Write filtered FASTA
write.fasta(sequences = filtered_fasta, names = seq_name$name[to_keep], file.out = paste0(opt[["outfile"]], "_filtered.fasta"))
