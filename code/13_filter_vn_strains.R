## load packages
required_packages <- c("tidyverse", "optparse")

# Function to check and install packages
# and load necessary utility functions
source("code/_headers.R")

option_list <- list(
  make_option(c("-c", "--csv"), type = "character", help = "Input csv file."),
  make_option(c("-o", "--outfile"), type = "character", help = "Base name for output files.")
)

opt <- parse_args(OptionParser(option_list = option_list))

##########################################################
## main
##########################################################

# Load infoTbl csv
info_tbl <- safe_read_file_param(opt[["csv"]], read_csv, required = TRUE)

# filter VN strains
vn_strains <- info_tbl %>%
  filter(Country == "Vietnam") %>%
  pull(Sequence_name)

# save filtered results
write_lines(vn_strains, opt[["outfile"]])
