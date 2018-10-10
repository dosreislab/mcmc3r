# Process object with raw data previously saved in rda format (see canids19x29.raw.R)
# for the carnivores specimens and perform Procrustes alignment with
# Morpho.

# 1. Convert raw data into class matrix and add row names
C.mat              <- as.matrix( carnivores19x29.raw[ , 3:dim( carnivores19x29.raw )[2] ] )
row.names( C.mat ) <- carnivores19x29.raw[ ,1 ]

# 2. Save this object that will be used in some examples in the package
C.mat.unal <- C.mat
devtools::use_data( C.mat.unal )

# 3.Run mcmc3r::matrix2array to generate an object of
# class array readable by Morpho
# Save the object as it will be used in some examples in the package
C.arr      <- matrix2array( X = C.mat, coords = 3 )
C.arr.unal <- C.arr
devtools::use_data( C.arr.unal )

# 4. Use Morpho::procSym for Procrustes alignment
# 4.1. Get right and left vectors with corresponding symmetric
# landmarks
right <- c( 11, 22, 13, 19, 15, 20, 24, 5, 7, 2, 9, 26, 29 )
left  <- c( 10, 21, 12, 17, 14, 18, 23, 4, 6, 1, 8, 25, 28 )
pairedLM <- cbind( left, right )
# 4.2. Run Morpho:procSym
C.PS   <- Morpho::procSym( dataarray = C.arr,
                           pairedLM = pairedLM )

# 5. Save morphological alignment with 19 carnivores in an rda
#    object (class matrix)
C <- array2matrix( X = C.PS$rotated, coords = 3 )
devtools::use_data( C )

# 6. Save procSym object output by Morpho, which mean shape is going to
# be used to align the Vulpes vulpes specimens later in the analysis
devtools::use_data( C.PS )


