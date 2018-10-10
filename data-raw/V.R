# Process object with raw data previously generated and saved in rda format (see vulpes21x29.raw.R)
# for the foxes ("Vulpes vulpes") specimens and perform Procrustes alignments with
# Morpho.

# 1. Convert to class matrix
V.mat <- as.matrix( vulpes21x29.raw[ , 2:dim( vulpes21x29.raw )[2] ] )

# 2. Add row names and save object to be used in examples in the package
row.names( V.mat ) <- vulpes21x29.raw[ , 1 ]
V.mat.unal         <- V.mat
devtools::use_data( V.mat.unal )

# 3. Run matrix2array to generate an object of
#    class array readable by Morpho and then save it as rda object
V.arr      <- matrix2array( X = V.mat, coords = 3 )
V.arr.unal <- V.arr
devtools::use_data( V.arr.unal )

# 4. We have the 21 foxes specimens in object "V.arr". 
#    So now we remove the common fox "Vulpes_1" in C.PS from V.arr and 
#    align it to the mean shape of C.PS so the "Vulpes vulpes" have
#    the same orientation than the rest of carnivores
#    Save the procSym object
V.arr.nov1            <- V.arr[,,-1]
V.PS.nov1             <- Morpho::align2procSym( x = C.PS, newdata = V.arr.nov1 )
dimnames( V.PS.nov1 ) <- dimnames( V.arr.nov1 )
devtools::use_data( V.PS.nov1 )

# 5. Get the vulpes (Vulpes_1) from C.PS and add it to V.PS array
V.PS              <- array( dim = c( 29, 3, 21 ) )
dimnames( V.PS )  <- list( paste( "lmk", seq( 1:29 ), sep="" ),
                                  c( "x", "y", "z" ),
                                  c( "Vulpes_1",
                                     dimnames( V.PS.nov1 )[[3]] ) )
V.PS[,,1]    <- C.PS$rotated[,,13] # Vulpes_1
V.PS[,,2:21] <- V.PS.nov1

# 6. Save alignment with 21 foxes (Vulpes vulpes) in rda format (class matrix)
V <- array2matrix( X = V.PS, coords = 3 )
devtools::use_data( V )

