# Load raw data from inst/extdata

vulpes21x29.raw <- read.table(file="inst/extdata/vulpes21x29.csv", sep=",")

# Save raw data in rda format

devtools::use_data(vulpes21x29.raw)
