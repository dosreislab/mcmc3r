#' Generate a phylip file for MCMCTree
#'
#' @description
#' Generate an alignment file in phylip format for MCMCTree.
#' The option "seqfile" in the control file used by MCMCTree
#' should read the path to the file output by this function.
#'
#' @param filename Character, name for the output file.
#'
#' @param proc 3D array (p x k x n) or matrix (n x p), landmark coordinates
#' after a Procrustes analysis (see details).
#'
#' @param coords Integer, only needed if \code{class(proc) == "array"}
#' (see details).
#'
#' @param names (Optional) list, species name included in the morphological
#' alignment (see examples B and C).
#'
#' @param ages (Optional) list, ages of the species included in the
#' morpholical alignment (see example C).
#'
#' @param popvar Vector of integers, population variance (see details).
#'
#' @param R Matrix, correlation matrix.
#'
#' @param method Character, method to decompose the inverse of the
#' correlation matrix R, \code{eigen} or \code{chol} (see details).
#'
#' @details
#'
#' If \code{class(proc) == "array"}, the array has format p x k x n, where
#' 'p' is the number of landmarks, 'k' the number of coordinates,
#' and 'n' the number of specimens. See \code{canids19x29.array} in the
#' \code{data} directory for an example.
#' You also need to provide \code{coords = 2} (2D landmark points) or
#' \code{coords = 3} (3D landmark points) if the \code{proc} object is an array.
#'
#' If \code{class(proc) == "matrix"}, the matrix has 'n' rows, one for specimen,
#' and 'p' columns.Each landmark can be given in 2D or 3D. For instance,
#' if the landmarks are 3D, the first 3 columns will be the
#' coordinates x, y, and z for the first landmark, the next 3
#' columns for the second landmark, and so on. See \code{canids19x29.matrix} in the
#' \code{data} directory for an example.
#'
#' If \code{names} is not provided, the name for each species will be
#' "Species_1", "Species_2", and so on.
#'
#' The parameter \code{popvar} accounts for the population variance. First of all,
#' a PA regarding the landmark points collected from different specimens
#' belonging to different species should be carried out, which will result into
#' a matrix of superimposed coordinates, say \code{M}. If there are many specimens of
#' one of the species included in \code{M}, then another PA is carried out with
#' the landmark points collected from these specimens belonging to the same species.
#' The resulting matrix, say \code{V}, is used to obtain the population variance.
#' The variance-covariance matrix of \code{V}, say \code{VCV}, is first calculated
#' and, afterwards, the diagonal of \code{VCV}, the vector of variances, is used
#' to scale the matrix \code{M} according to the standard deviation.
#' Therefore, \code{popvar = diag(VCV)}. If this parameter is not provided, then
#' the population variance is assumed to be 0.
#'
#' If \code{method = "eigen"}, the inverse of the correlation matrix provided,
#'  \code{R}, will be decomposed following the spectral eigendecomposition
#' (see \code{eigen}).
#' If \code{method = "chol"}, the inverse of \code{R} will be decomposed
#' following the Cholesky decomposition (see \code{chol}).
#'
#' @author Sandra Alvarez-Carretero
#'
#' @examples
#' # A.1) Providing only the morphological alignment (proc) after the
#' #      Procrustes analysis (PA) in an object of class \"array\".
#'
#'        write.morpho( filename = "seqfile.aln", proc = canids19x29.array,
#'                      coords = 3 )
#'
#' # A.2) Providing only the morphological alignment (proc) after the
#' #      Procrustes analysis (PA) in an object of class \"matrix\".
#'
#'        write.morpho( filename = "seqfile.aln", proc = canids19x29.matrix,
#'                      coords = 3 )
#'
#' # B) Providing the morphological alignment (proc) after the
#' #    PA in an object of class \"array\" and a list with the
#' #    names of the species
#'
#'      names <- list( sp1  = "Ael_sp.", sp2  = "Can_dir", sp3  = "Epi_hay", sp4  = "Hes_sp.",
#'                     sp5  = "Mes_cor", sp6  = "Tom_sp.", sp7  = "Enh_pah", sp8  = "Cuo_alp",
#'                     sp9  = "Spe_ven", sp10 = "Can_lup", sp11 = "Cer_tho", sp12 = "Oto_meg",
#'                     sp13 = "Vul_vul", sp14 = "Urs_ame", sp15 = "Ail_ful", sp16 = "Nan_bio",
#'                     sp17 = "Par_her", sp18 = "Tha_won", sp19 = "Smi_fat"
#'                   )
#'
#'      write.morpho( filename = "seqfile.aln", proc = canids19x29.array,
#'                    coords = 3, names = names )
#'
#' # C) Providing the morphological alignment (proc) after the
#' #    PA in an object of class \"array\", a list with the names
#' #    of the species, and a list with the ages of the species
#'
#'      names <- list( sp1  = "Ael_sp.", sp2  = "Can_dir", sp3  = "Epi_hay", sp4  = "Hes_sp.",
#'                     sp5  = "Mes_cor", sp6  = "Tom_sp.", sp7  = "Enh_pah", sp8  = "Cuo_alp",
#'                     sp9  = "Spe_ven", sp10 = "Can_lup", sp11 = "Cer_tho", sp12 = "Oto_meg",
#'                     sp13 = "Vul_vul", sp14 = "Urs_ame", sp15 = "Ail_ful", sp16 = "Nan_bio",
#'                     sp17 = "Par_her", sp18 = "Tha_won", sp19 = "Smi_fat"
#'                   )
#'
#'      ages <- list( sp1  = 15.97, sp2  =  1.80, sp3  = 13.6,  sp4  = 39.74,
#'                    sp5  = 30.80, sp6  = 15.97, sp7  = 30.80, sp8  =  0,
#'                    sp9  =  0,    sp10 =  0,    sp11 =  0,    sp12 =  0,
#'                    sp13 =  0,    sp14 =  0,    sp15 =  0,    sp16 =  0,
#'                    sp17 =  0,    sp18 =  0,    sp19 =  0
#'                   )
#'
#'      write.morpho( filename = "seqfile.aln", proc = canids19x29.array,
#'                    coords = 3, names = names, ages = ages )
#'
#' # D) Providing an object of class \"array\" after having
#' #    carried out a PA. As an example, we use the function
#' #    geomorph::gpagen, but you can use your preferred
#' #    function meanwhile the format of the \"proc\" object
#' #    (class array) for write.morpho is
#' #    "p landmarks x k coordinates x n species"
#'
#'      df <- canids19x29.raw
#'
#'      # Process raw data
#'
#'      mm <- df[,2:dim( df )[2]] # 19 sp x 87 coords (87/3=29 lmks)
#'      rownames( mm ) <- df[,1]
#'      mm <- as.matrix( mm )
#'
#'      # Convert matrix into 3D array
#'
#'      ma <- matrix2array( proc = mm, coords = 3 )
#'
#'      # Get procrustes analysis done
#'
#'      ma.paln <- geomorph::gpagen( ma )
#'
#'      # Check object is class \"array\"
#'
#'      class( ma.paln$coords )
#'
#'      # Run write.morpho
#'
#'      write.morpho( filename = "seqfile.aln", proc = ma.paln$coords, coords = 3 )
#'
#' # E) Providing the morphological alignment (proc) after the
#' #    PA in an object of class \"matrix\" and population variance
#'
#'      write.morpho( filename = "seqfile.aln", proc = canids19x29.array,
#'                    coords = 3, popvar = var.fx )
#'
#' # F) Providing the morphological alignment (proc) after the
#' #    PA in an object of class \"matrix\", population variance,
#' #    the correlation matrix, and choosing to use the Cholesky
#' #    decomposition
#'
#'      write.morpho( filename = "seqfile.aln", proc = canids19x29.array,
#'                    coords = 3, popvar = var.fx, R = R, method = "chol" )
#'
#'
#'
#'
#' @export
write.morpho <- function( filename, proc, coords = c( 2, 3 ), names = NULL, ages = NULL,
                          popvar = NULL, R = NULL, method = c( "eigen", "chol" ) ) {

  #\\ Check a name for the output file has been given

  if ( missing( filename ) ){
    stop( "Please add a name for the output file" )
  }

  #\\ Check object with lmks is provided

  if ( missing( proc ) ){
         stop( "Please use an object of class \"matrix\" or \"array\" with the result of a Procrustes analysis" )
  }

  if ( class( proc ) == "array" | class( proc ) == "matrix" ){

    #\\ Check class of object with lmks and
    #\\ get variables

    if( class( proc ) == "array" ){

      # The coords info is needed

      if ( is.null( coords ) | length( coords ) != 1 | class( coords ) != "numeric" ){
        stop( "Please enter coords = 2 or coords = 3 depending on using 2D or 3D landmarks, respectively" )
      }

      chars  <- dim( proc )[1] * dim( proc )[2]
      num.sp <- dim( proc )[3]

    }

    else if( class( proc ) == "matrix" ){

      num.sp <- dim( proc )[1]
      chars  <- dim( proc )[2]

    }


    #\\ If names are provided...

    if ( !is.null( names ) ){

      # Check class for names

      if ( class( names ) != "list" ){

        stop( "You need to provide an object of class list with the species included in the alignment file.
               E.g. species <- list( sp1 = \"sp1\", sp2 = \"sp2\"" )
      }

      # Get species names

      names   <- unlist( names )

      #-- If ages are not provided ...

      if ( is.null( ages ) ){

        # Get spaces to justify text

        num.spaces <- max( nchar( names ) ) - nchar( names )
        x          <- sprintf( paste( "% ", num.spaces, "s", sep = ""), c( "" ) )

      }

      #-- If ages are provided ...

      else{

        # Check class for ages

        if ( class( ages ) != "list" ){

          stop( "You need to provide an object of class list with the ages of the species included in the alignment file.
                 E.g. ages <- list( ag1=30, ag2=15" )
        }

        # Put together names and transformed ages and get spaces

        ages       <- unlist( ages )
        names      <- paste( names, - ages + max( ages ) + 0.01, sep = "^" ) # 0.01 = ct for MCMCTree
        num.spaces <- max( nchar( names ) ) - nchar( names )
        x          <- sprintf( paste( "% ", num.spaces, "s", sep = "" ), c( "" ) )

      }

    }


    #\\ If names are not provided...

    else{

      # Name species as "Species_1", "Species_2", etc.

      names <- paste( "Species_", seq( 1:num.sp ), sep = "" )

      #-- If ages are not provided ...

      if ( is.null(ages) ){

        # Get spaces to justify text

        num.spaces <- max( nchar( names ) ) - nchar( names )
        x          <- sprintf( paste( "% ", num.spaces, "s", sep = "" ), c( "" ) )

      }

      #-- If ages are provided ...

      else{

        # Check class for ages

        if ( class(ages) != "list" ){

          stop( "You need to provide an object of class list with the ages of the species included in the alignment file.
              E.g. ages <- list(ag1=30, ag2=15" )
        }

        # Put together names and transformed ages and get spaces

        ages       <- unlist( ages )
        names      <- paste( names, - ages + max( ages ) + 0.01, sep = "^" ) # 0.01 = ct for MCMCTree
        num.spaces <- max( nchar(names) ) - nchar( names )
        x          <- sprintf( paste( "% ", num.spaces, "s", sep = "" ), c( "" ) )

      }

    }

    #\\ Get names with data

    names <- paste( names, x, sep = "    " )

    if ( class( proc ) != "array" ){

      rownames( proc ) <- names

    }

    else{

      faln             <- array2matrix( proc = proc, coords = coords )
      rownames( faln ) <- names
      proc             <- faln

    }

    #\\ Get popvar

    if ( !is.null( popvar ) ){

      # If the variances are different...
      if ( length( unique( popvar ) ) != 1 ){

        proc <- proc %*% diag( 1 / sqrt( popvar ) )
        pvar <- 1 # The matrix has been scaled, so now

      }

      # If the variances are all equal...
      if ( length( unique( popvar ) ) == 1 ){

        pvar <- 1 # No need for scale

      }


    }

    else{

      pvar <- 0 # No variance at all

    }

    #\\ Decompose R and get lnd

    if ( !is.null( R ) ){

      # Check method is not missing

      if ( missing(method) ){
        stop( "Please select a method to decompose the correlation matrix")
      }

      # If R is not class "matrix", stop

      if ( class( R ) != "matrix" ){

        stop( "Please provide an object of class \"matrix\" for the correlation matrix" )
      }

      # Match argument

      method <- match.arg(method)

      if ( method == "chol" ){

        # chol returns upper triangular matrix

        U <- chol( R )
        L <- t( U )

        # R = L %*% U = L %*% t( L )
        # solve( R ) = t( solve(L) ) %*% solve( L )
        # solve( R ) = t( A )        %*% A
        #
        # all.equal( t( solve( L ) ) %*% solve( L ), solve( R ) )

        Linv <- solve( L )

        # Transform data "M", object proc, so
        # Z = M %*% t( A ) = proc %*% t( Linv )

        Z <- proc %*% t( Linv )

      }

      else if ( method == "eigen" ){

        # eigen returns a list with the eigenvectors
        # and the eigenvalues
        #
        # solve( R ) = Rinv = t( A ) %*% A
        #
        # t( A ) = V %*% D
        # t( A ) = eigen( Rinv )$vectors %*% diag( sqrt( eigen( Rinv )$values ) )
        #
        # all.equal( tA %*% t( tA ), solve( R ) )

        tA <- eigen( solve(R) )$vectors %*% diag( sqrt( eigen( solve(R) )$values ) )

        # Transform data "M", object proc, so
        # Z = M %*% t( A ) = proc %*% t( A )

        Z <- proc %*% tA

      }

      # Get logarithm of the determinant of the correlation matrix
      # This variable is used to correct the likelihood value when the data have
      # been transformed to account for the correlation.
      # This is done to allow the correct calculation of bayes factors.

      lnd <- determinant( R )$modulus # lnd = lod(det(R))

    }

    else{

      # There is no correlation at all

      lnd <- 0
      Z   <- proc

    }

    #\\ Generate file

    string <- paste( length( names ), chars, "M", pvar, lnd, sep = "  " )
    cat( "\n", file = filename )
    write( paste( "  ", string, sep = "" ), file = filename, append = T )
    #cat(paste("  ", string, sep=""))
    cat( "\n", file = filename, append = T )
    #cat("\n")
    write.table( Z, file = filename, append = T,
                 sep = " ", row.names = T, col.names = F, quote = FALSE )
    cat( "\n", file = filename, append = T )
    #Z

  }

  #\\ If object proc is not a matrix nor an array...

  else{

    stop( "You need to use an object of class \"matrix\" or \"array\" with the result of a Procrustes analysis" )

  }

}






#' Convert a matrix into a 3D array
#'
#' @description
#' Convert a matrix with landmark points into a 3D array object.
#'
#' @param proc Matrix (n x p), 'p' landmark points for 'n' specimens
#' (see details).
#'
#' @param coords Integer, 2 or 3 for 2D or 3D landmarks, respectively.
#'
#' @details
#'
#' The matrix has format n x p, with 'n' rows, one for specimen, and 'p' columns.
#' Each landmark can be given in 2D or 3D. For instance,
#' if the landmarks are 3D, the first 3 columns will be the
#' coordinates x, y, and z for the first landmark, the next 3
#' columns for the second landmark, and so on. See \code{canids19x29.matrix} in the
#' \code{data} directory for an example.
#'
#' @return
#'
#' An object of class array with format p x k x n, where 'p' is the number of
#' landmarks, 'k' the number of coordinates, and 'n' the number of specimens.
#'
#' Note that if the matrix provided does not have rownames, the specimens in the
#' returned array (dimension 'n') will be labelled as \"1\", \"2\", and so on.
#'
#' @author Sandra Alvarez-Carretero
#'
#' @export

matrix2array <- function( proc, coords = c( 2, 3 ) ){

  #\\ Check object with lmks is provided

  if ( missing( proc ) ){
    stop( "Please provide an object of class \"matrix\"" )
  }

  #\\ Only proceed it the object is a matrix

  if ( class( proc ) == "matrix" ){

    if ( length( dim( proc ) ) != 2 ){
      stop( "Please provide a matrix \"n x p\", with 'n' specimens and 'p' characters" )
    }

    if ( is.null( coords ) | length( coords ) != 1 | class( coords ) != "numeric" ){
      stop( "Please enter coords = 2 or coords = 3 depending on using 2D or 3D landmarks, respectively" )
    }

    chars       <- dim( proc )[2]
    num.sp      <- dim( proc )[1]
    list.coords <- c("x", "y", "z")

    # Create an empty array to store the coordinates in the format
    # p x k x n (chars x coord x num.sp)

    ma             <- array( dim = c( chars / coords, coords, num.sp ) )
    dimnames( ma ) <- list( paste( "lmk", seq( 1:( chars/coords ) ), sep="" ),
                            list.coords[1:coords],
                            rownames( proc ) )

    # Select x, y, z positions

    if ( coords == 3 ){

      xi <- seq( from = 1, to = chars, by = 3 )
      yi <- seq( from = 2, to = chars, by = 3 )
      zi <- seq( from = 3, to = chars, by = 3 )

      # Fill in array p x k x n

      for (i in 1:num.sp) {
        ma[,1,i] <- unlist( proc[i,xi] )
        ma[,2,i] <- unlist( proc[i,yi] )
        ma[,3,i] <- unlist( proc[i,zi] )
      }

    }

    else if ( coords == 2 ){

      xi <- seq( from = 1, to = chars, by = 2 )
      yi <- seq( from = 2, to = chars, by = 2 )

      # Fill in array p x k x n

      for (i in 1:num.sp) {
        ma[,1,i] <- unlist( proc[i,xi] )
        ma[,2,i] <- unlist( proc[i,yi] )
      }

    }


  #\\ Return 3D array

    return(ma)

  }

  else{
    stop( "Please use an object of class \"matrix\"" )
  }

}







#' Convert a 3D array into a matrix
#'
#' @description
#' Convert a 3D array with landmark points into an object of class matrix.
#'
#' @param proc 3D array, 'p' landmark points, 'k' coordinates, and
#' 'n' specimens.
#'
#' @param coords Integer, 2 or 3, for 2D or 3D landmarks, respectively.
#'
#' @details
#'
#' See \code{canids19x29.array} in the \code{data} directory for an example
#' of the format of a 3D array.
#'
#' @return
#'
#' An object of class matrix, with 'n' rows, one for specimen, and 'p' columns.
#' Each landmark can be given in 2D or 3D. For instance,
#' if the landmarks are 3D, the first 3 columns will be the
#' coordinates x, y, and z for the first landmark, the next 3
#' columns for the second landmark, and so on.
#'
#' Note that if the 'n' dimension (specimens) of the array provided does
#' not have names, #' the specimens in the returned matrix will be
#' labelled as \"1\", \"2\", and so on.
#'
#' @author Sandra Alvarez-Carretero
#'
#' @export

array2matrix <- function( proc, coords = c( 2, 3 ) ){

  #\\ Check object with lmks is provided

  if ( missing( proc ) ){
    stop( "Please provide an object of class \"array\"" )
  }

  #\\ Only proceed if the object is an array

  if ( class( proc ) == "array" ){

    if ( length( dim( proc ) ) != 3 ){
      stop( "Please provide a 3D array \"p x k x n\", with 'p' landmarks, 'k' coordinates, and 'n' specimens" )
    }

    if ( is.null( coords ) | length( coords ) != 1 | class( coords ) != "numeric" ){
      stop( "Please enter coords = 2 or coords = 3 depending on using 2D or 3D landmarks, respectively" )
    }

    chars  <- dim( proc )[1]
    coords <- coords
    num.sp <- dim( proc )[3]

    faln <- matrix( 0, nrow = num.sp, ncol = chars * coords )

    if ( coords == 2 ){

      xi <- seq( from = 1, to = chars, by = 2 ); yi <- seq( from = 2, to = chars, by = 2 )

      for ( i in 1:num.sp ) {
        faln[i,xi] <- proc[,1,i]
        faln[i,yi] <- proc[,2,i]
      }

      rownames(faln) <- dimnames( proc )[[ 3 ]]

    }

    else if ( coords == 3 ){
      xi <- seq( from = 1, to = chars * coords, by = 3 )
      yi <- seq( from = 2, to = chars * coords, by = 3 )
      zi <- seq( from = 3, to = chars * coords, by = 3 )

      for ( i in 1:num.sp ) {
        faln[i,xi] <- proc[,1,i]
        faln[i,yi] <- proc[,2,i]
        faln[i,zi] <- proc[,3,i]

      }

      rownames(faln) <- dimnames( proc )[[ 3 ]]

    }

  #\\ Return matrix

    return(faln)

  }

  else{
    stop( "Please use an object of class \"array\"" )
  }

}












