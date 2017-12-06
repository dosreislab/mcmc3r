# Load raw data from inst/extdata

vulpes21x29.raw <- read.table(file="inst/extdata/vulpes21x29.csv", sep=",")

# Get rownames and delete first column of raw data
# to get a matrix object

rownames(vulpes21x29.raw) <- vulpes21x29.raw[,1]
vulpes21x29.raw.matrix    <- vulpes21x29.raw[,2:dim(vulpes21x29.raw)[2]]
vulpes21x29.raw.matrix    <- as.matrix(vulpes21x29.raw.matrix)

# Run morpho::matrix2array to generate an object of
# class array readable by geomorph

vulpes21x29 <- morpho::matrix2array(proc=vulpes21x29.raw.matrix, coords=3)

# Run geomorph::gpagen and save coordinates in an
# array (rda format)

proc              <- geomorph::gpagen(vulpes21x29)
vulpes21x29.array <- proc$coords

devtools::use_data(vulpes21x29.array)


