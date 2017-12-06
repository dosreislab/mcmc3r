# Load raw data from inst/extdata

canids19x29.raw <- read.table(file="inst/extdata/canids19x29.csv", sep=",")

# Get rownames and delete first column of raw data
# to get a matrix object

rownames(canids19x29.raw) <- canids19x29.raw[,1]
canids19x29.raw.matrix    <- canids19x29.raw[,2:dim(canids19x29.raw)[2]]
canids19x29.raw.matrix    <- as.matrix(canids19x29.raw.matrix)

# Run morpho::matrix2array to generate an object of
# class array readable by geomorph

canids19x29 <- morpho::matrix2array(proc=canids19x29.raw.matrix, coords=3)

# Run geomorph::gpagen and save coordinates

proc              <- geomorph::gpagen(canids19x29)
canids19x29.array <- proc$coords

# Save matrix as rda format

canids19x29.matrix <- morpho::array2matrix(proc=canids19x29.array, coords=3)

devtools::use_data(canids19x29.matrix)

