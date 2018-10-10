# Load raw data from inst/extdata
carnivores19x29.raw <- read.table ( file = "inst/extdata/carnivores19x29.tsv", header = F, skip = 1, dec = ".", sep = "\t" )

# Save raw data in rda format
devtools::use_data( carnivores19x29.raw )
