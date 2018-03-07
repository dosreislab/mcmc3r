# Load raw data from inst/extdata

canids19x29.raw <- read.table( file = "inst/extdata/canids19x29.csv", sep = "," )

# Get rownames and delete first and second columns of raw data
# to get a matrix object

rownames( canids19x29.raw ) <- canids19x29.raw[ ,1 ]
canids19x29.raw.matrix      <- canids19x29.raw[ , 3:dim( canids19x29.raw )[2] ]
canids19x29.raw.matrix      <- as.matrix( canids19x29.raw.matrix )

# Run morpho::matrix2array to generate an object of
# class array readable by geomorph

canids19x29 <- morpho::matrix2array( X = canids19x29.raw.matrix, coords = 3 )

# Run geomorph::gpagen and save coordinates

C.PS <- geomorph::gpagen( canids19x29 )

# Remove vulpes coords from proc.C
index.V   <- which( dimnames( C.PS$coords )[[3]] == "Vul_vul" )
C.PS.filt <- C.PS$coords[ , , -c( index.V ) ]

# Calculate mean shape of all Vulpes sp. specimens
# We use the object previously created with all the vulpes
# aligned in a previous PS analysis (see data-raw/vulpes21.29.matrix.R)
mean.shape.V <- mshape( vulpes21x29.array )

# Create an array to have the filtered carnivores
# data set (no Vulpes) and add the mean shape of all Vulpes vulpes
# to work as the mean shape of the population of Vulpes vulpes
C.PS.filt.arr              <- array( dim = c( 29, 3, 19 ) )
dimnames( C.PS.filt.arr )  <- list( paste( "lmk", seq( 1:29 ), sep="" ),
                                    c( "x", "y", "z" ),
                                    c( dimnames( C.PS.filt )[[3]],
                                       "Vul_vul" ) )
C.PS.filt.arr[ , , 1:18 ]  <- C.PS.filt
C.PS.filt.arr[ , , 19 ]    <- mean.shape.V

# Perform last PS analysis to get morphological alignment
M.PS              <- geomorph::gpagen( C.PS.filt.arr )
canids19x29.array <- M.PS$coords

# Save array as rda format

devtools::use_data( canids19x29.array )



