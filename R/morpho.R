#' Generate a phylip file for MCMCTree
#'
#' @description
#' Generate an alignment file in phylip format for MCMCTree.
#' The option "seqfile" in the control file used by MCMCTree
#' should read the path to the file output by this function.
#'
#' @param proc Array ( p landmarks x k coordinates x n species),
#' matrix (n species x p landmarks), or data.frame (n species x p landmarks)
#' object with ONLY the landmarks of the bones of the species
#' that are included in the alignment. Note that the Procrustes
#' superimposition should have been done before generating
#' this file.
#'
#' @param names (optional) List with the species name included in the
#' morphological alignment. If not provided, the name for each species
#' will be "Species_1", "Species_2", and so on.
#' E.g. sp.list <- list( sp1 = "sp.1", sp2 = "sp.2", ... )
#'
#' @param ages (optional) List of the ages of the species included
#' in the morpholical alignment.
#' E.g. age.list <- list( sp1 = 15, sp2 = 0, ... )
#'
#' @return
#' Alignment file in phylip format for MCMCTree.
#'
#' @author Sandra Alvarez-Carretero
#'
#' @examples
#'
#' # A) Providing only the morphological alignment (proc)
#'
#'      data()
#'      write.phylo( proc = coords.proc )
#'
#' # B) Providing the morphological alignment (proc) and a list
#' #    with the names of the species
#'
#'      data()
#'      names <- list( sp1 = "Sp.011", sp2 = "Sp.23", sp3 = "Sp.1333", sp4 = "Sp.00127",
#'                     sp5 = "Sp.23625", sp6 = "Sp.1304", sp7 = "Sp.543", sp8 = "Sp.02",
#'                     sp9 = "Sp.234" )
#'
#'      write.phylo( proc = coords.proc, names = names )
#'
#' # C) Providing the morphological alignment (proc), a list
#' #    with the names of the species, and a list with the
#' #    ages of the species
#'
#'      data()
#'      names <- list( sp1 = "Sp.011", sp2 = "Sp.23", sp3 = "Sp.1333", sp4 = "Sp.00127",
#'                     sp5 = "Sp.23625", sp6 = "Sp.1304", sp7 = "Sp.543", sp8 = "Sp.02",
#'                     sp9 = "Sp.234" )
#'
#'      ages <- list( sp1 = 15, sp2 = 30, sp3 = 8, sp4 = 0,
#'                    sp5 = 0, sp6 = 0, sp7 = 12, sp8 = 0,
#'                    sp9 = 25.5 )
#'
#'      write.phylo( proc = coords.proc, names = names, ages = ages )
#'
#' # D) Providing an object of class \"array\" after having
#' #    carried out a Procrustes analysis (PA). As an example,
#' #    we use the function geomorph::gpagen, but you can use
#' #    your preferred function meanwhile the format of the
#' #    \"proc\" object for write.phylo is "p landmarks x k coordinates x n species"
#'
#'      data()
#'      df <- coords.raw
#'      mm <- df[,2:( dim( df )[2] )] # 9sp x 144coords (144/3=48 lmk)
#'      rownames( mm ) <- df[,1]
#'
#'      # Create an empty array to store the coordinates in the format
#'      # p x k x n (num.coords x coor3D x ns)
#'
#'      num.coords <- dim(mm)[2]
#'      coords <- 3 # Change according to your data set
#'
#'      ns <- 9
#'      ma <- array( dim = c( num.coords / coords, coords, ns ) ) # 48 lmks, 3D, 9 specimens
#'      dimnames( ma ) <- list( paste( "lmk", seq( 1:( num.coords/coords ) ), sep="" ),
#'                              c("x", "y", "z"),
#'                              df[,1])
#'
#'      # Select x, y, z positions
#'
#'      xi <- seq( from=1, to = num.coords, by = 3)
#'      yi <- seq( from=2, to = num.coords, by = 3)
#'      zi <- seq( from = 3, to = num.coords, by = 3)
#'
#'      # Fill in array p x k x n
#'
#'      for (i in 1:ns) {
#'         ma[,1,i] <- unlist( mm[i,xi] )
#'         ma[,2,i] <- unlist( mm[i,yi] )
#'         ma[,3,i] <- unlist( mm[i,zi] )
#'      }
#'
#'      # Get procrustes analysis done
#'
#'      ma.paln <- geomorph::gpagen(ma)
#'
#'      # Check object is class \"array\"
#'
#'      class( ma.paln$coords )
#'
#'      # Run write.phylo
#'
#'      write.phylo( proc = ma.paln$coords )
#'
#' @export
write.phylo <- function( proc , names=NULL, ages=NULL ) {

  #\\ Check object with lmks is provided

  if ( missing( proc ) ){
         stop( "Please use an object of class \"matrix\", \"data.frame\", or \"array\" with the result of a Procrustes analysis" )
  }

  if ( class( proc ) == "array" | class( proc ) == "matrix" | class( proc ) == "data.frame" ){


    #\\ Check class of object with lmks and
    #\\ get variables

    if( class( proc ) == "array" ){

      chars  <- dim( proc )[1] * dim( proc )[2]
      num.sp <- dim( proc )[3]

    }

    else if( class( proc ) == "data.frame" | class( proc ) == "matrix" ){

      if ( class(proc) == "data.frame" ){
        proc <- as.matrix( proc )
      }

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

        # Put together names, spaces, and transformed ages

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

        # Put together names, spaces, and transformed ages

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

      faln <- matrix( 0, nrow = num.sp, ncol = chars )
      xi <- seq( from = 1, to = chars, by = 3 ); yi <- seq( from = 2, to = chars, by = 3 ); zi <- seq( from = 3, to = chars, by = 3 )

      for ( i in 1:num.sp ) {
        faln[i,xi] <- proc[,1,i]
        faln[i,yi] <- proc[,2,i]
        faln[i,zi] <- proc[,3,i]
      }

      rownames( faln ) <- names
      proc <- faln

    }


    #\\ Generate file

    string <- paste( length( names ), chars, "M", sep = "  " )
    cat( "\n", file = "seqfile.aln" )
    write( paste( "  ", string, sep = "" ), file = "seqfile.aln", append = T )
    #cat(paste("  ", string, sep=""))
    cat( "\n", file = "seqfile.aln", append = T )
    #cat("\n")
    write.table( proc, file = "seqfile.aln", append = T,
                 sep = " ", row.names = T, col.names = F, quote = FALSE )
    cat( "\n", file = "seqfile.aln", append = T )
    #proc

  }

  #\\ If object proc is not a matrix, a data.frame, nor an array...

  else{

    stop( "You need to use an object of class \"matrix\", \"data.frame\", or \"array\" with the result of a Procrustes analysis" )

  }

}

