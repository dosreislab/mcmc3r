# Load raw data from inst/extdata

canids19x29.raw <- read.table(file="inst/extdata/canids19x29.csv", sep=",")

# Save raw data in rda format

devtools::use_data(canids19x29.raw)
