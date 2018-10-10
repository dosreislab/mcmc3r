# Load raw data from inst/extdata
vulpes21x29.raw <- read.table( file = "inst/extdata/vulpes21x29.tsv", header = T, dec = ".", sep = "\t" )

# Save raw data in rda format
devtools::use_data( vulpes21x29.raw )
