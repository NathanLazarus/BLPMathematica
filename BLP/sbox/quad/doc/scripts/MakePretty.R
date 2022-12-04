# MakePretty - takes solver output and creates nice LaTeX tables
#
# modification history
# --------------------
# 12jan2012 bss cleaned up captions for final version of thesis.
# 25oct2010 bss added support for Monomial, 5Good, etc. in Short format.
# 25aug2010 bss added support for point estimates & std errs.
# 25aug2010 bss updated for new tree structure.
# 18jun2010 bss added ND case for one start with multiple draws of nodes
# 18may2010 bss updated for CPU Time, general pathnames
# 11may2010 bss Greek column headers, subsets of table
# 10may2010 bss written.
#


#-----------------------------------------------------------
# Setup
library( xtable )

setwd('/Users/bss/sbox/quad/doc/tables')


#-----------------------------------------------------------
# Create nice LaTeX table for Solver Results and f_k for all data sets

MakePretty <- function( szName, szCap ) 
{
  szIn    <- paste( 'SolverResults-', szName, '.txt', sep='' )
  szOut   <- paste( 'SolverResults-', szName, '.inc', sep='' )
  szLabel <- paste( 'TabSolverResults', szName, sep='' )
#  szCap   <- paste( 'MPEC Results:', szName )

  df <- read.table( file = szIn )
  
  # reorder rows
  sol.df <- cbind( df[,1], df[,c(3,4)], df[,2], df[,5] )
  colnames( sol.df ) <- c( 'Dataset', 'EXIT', 'INFORM', 'f_k', 'CPU Time' )
  
  xt <- xtable( sol.df, caption=szCap, label=szLabel, digits=c(0,0,0,0,5,2) )
  
  print( xt, file=szOut, type='latex', math.style.negative=T, 
        include.rownames=F, floating=T, caption.placement='bottom',
        hline.after=c( -1, 0, seq(5,25,5) ) )
  
  return( sol.df )
}

#-----------------------------------------------------------
# Create nice LaTeX table for point estimates for all data sets

MakePointEst <- function( szName, szCap )
{
  szIn    <- paste( 'PointEst-', szName, '.txt', sep='' )
  szOut   <- paste( 'PointEst-', szName, '.inc', sep='' )
  szLabel <- paste( 'TabPointEst', szName, sep='' )
#  szCap   <- paste( 'Point Estimates :', szName )

  info.df <- MakePretty( szName, szCap )

  df <- read.table( file = szIn )
  
  # Get rid of negative sign on scale of random coefficients
  #sol.df <- cbind( df[,1], info.df[,3], df[,2:11] )
  sol.df <- cbind( df[,1], info.df[,3], df[,2:6], abs( df[,7:11] ) )

  colnames( sol.df ) <- c( 'Data', 'INFORM', "$\\theta_{11}$", "$\\theta_{12}$", 
                            "$\\theta_{13}$", "$\\theta_{14}$",
                            "$\\theta_{15}$", "$\\theta_{21}$", 
                            "$\\theta_{22}$", "$\\theta_{23}$",
                            "$\\theta_{24}$", "$\\theta_{25}$" )

  vIxSolDf <- seq( 3, length( colnames( sol.df ) ) )

  xt <- xtable( sol.df, caption=szCap, label=szLabel, digits=4, 
                display=c( 'd', 'd', 'd', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G' ) )
  
  print( xt, file=szOut, type='latex', math.style.negative=T, 
        include.rownames=F, floating=T, caption.placement='bottom',
        floating.environment='sidewaystable', hline.after=c( -1, 0, seq(5,25,5) ), sanitize.text.function=function(x){ x } )
  
#  return( sol.df )     # XXX
}


#-----------------------------------------------------------
# Create nice LaTeX table for Solver Results and f_k for one data set

MakePrettySmall <- function( szName, szCap, nSet ) 
{
  szIn    <- paste( 'SolverResults-', szName, '.txt', sep='' )
  szOut   <- paste( 'SolverResults-', szName, '-Short', nSet, '.inc', sep='' )
  szLabel <- paste( 'TabSolverResults', szName, sep='' )
#  szCap   <- paste( 'MPEC Results:', szName )

  # Create index vector to select correct slice of data
  nStarts <- 5
  vIx     <- seq( nStarts * ( nSet - 1 ) + 1, nStarts * nSet, 1 )

  df <- read.table( file = szIn )
  
  # reorder rows
  sol.df <- cbind( df[vIx,1], df[vIx,c(3,4)], df[vIx,2], df[vIx,5] )
  colnames( sol.df ) <- c( 'Dataset', 'EXIT', 'INFORM', 'f_k', 'CPU Time' )
  
  xt <- xtable( sol.df, caption=szCap, label=szLabel, digits=c(0,0,0,0,5,2) )
  
  print( xt, file=szOut, type='latex', math.style.negative=T, 
        include.rownames=F, floating=T, caption.placement='bottom' )
  
  return( sol.df )
}


#-----------------------------------------------------------
# Create nice LaTeX table for point estimates for all data sets

MakePointEstSmall <- function( szName, szCap, nSet )
{
  szIn    <- paste( 'PointEst-', szName, '.txt', sep='' )
  szOut1  <- paste( 'PointEst-', szName, '-Short', nSet, '.1', '.inc', sep='' )
  szOut2  <- paste( 'PointEst-', szName, '-Short', nSet, '.2', '.inc', sep='' )
  szLabel <- paste( 'TabPointEst', szName, sep='' )
#  szCap   <- paste( 'Point Estimates :', szName )

  # Create index vector to select correct slice of data
  nStarts <- 5
  vIx     <- seq( nStarts * ( nSet - 1 ) + 1, nStarts * nSet, 1 )
  info.df <- MakePretty( szName, szCap )

  df <- read.table( file = szIn )
  
  # Create dataset and get rid of negative signs on variance
  sol.df <- cbind( df[vIx,1], info.df[vIx,3], df[vIx,2:6], abs( df[vIx,7:11] ) )
  colnames( sol.df ) <- c( 'Data', 'INFORM', "$\\theta_{11}$", "$\\theta_{12}$", 
                            "$\\theta_{13}$", "$\\theta_{14}$",
                            "$\\theta_{15}$", "$\\theta_{21}$", 
                            "$\\theta_{22}$", "$\\theta_{23}$",
                            "$\\theta_{24}$", "$\\theta_{25}$" )

#  xt <- xtable( sol.df[,c(1,2,3:7)], caption=szCap, label=szLabel, digits=4, 
#                display=c( 'd', 'd', 'd', 'G', 'G', 'G', 'G', 'G' ) )
  xt <- xtable( sol.df[,c(3:7)], caption=szCap, label=szLabel, digits=2, 
                display=c( 'd', 'G', 'G', 'G', 'G', 'G' ) )
  
  print( xt, file=szOut1, type='latex', math.style.negative=T, 
        include.rownames=F, floating=T, caption.placement='bottom',
        sanitize.text.function=function(x){ x } )
  
#  xt <- xtable( sol.df[,c(1,2,8:12)], caption=szCap, label=szLabel, digits=4, 
#                display=c( 'd', 'd', 'd', 'G', 'G', 'G', 'G', 'G' ) )
  xt <- xtable( sol.df[,c(8:12)], caption=szCap, label=szLabel, digits=2, 
                display=c( 'd', 'G', 'G', 'G', 'G', 'G' ) )
  
  print( xt, file=szOut2, type='latex', math.style.negative=T, 
        include.rownames=F, floating=T, caption.placement='bottom',
        sanitize.text.function=function(x){ x } )
  
}


#===========================================================
# Tables with Point Estimates and Standard Errors
#===========================================================


#-----------------------------------------------------------
# Combine standard errors for multiple seeds
#
# RETURNS (seed x start) x (parameters (10) ) - e.g. five starts
#   stacked on top of each other for each seed.

LoadStdErr <- function( szBaseDir, nQuadType, nNodes )
{
  require( R.matlab )

  # Setup
  szBaseName <- 'Variance-QuadType'
  nSeed      <- 5
  nStart     <- 5
  nParam     <- 10

  vIxCol <- seq( 1, nParam )    # Index for columns of parameters

  # (a seed's start) x (standard errors for the ten parameters)
  mStdErr <- matrix( nrow=nStart*nSeed, ncol=nParam )

  # Read in data for each seed and combine
  for( ixSeed in seq( 1, nSeed ) ) {

    # Build filename of current Variance MATLAB file
    szSeedDir     <- paste( 'Seed', '000', ixSeed, sep='' )
    szVarFilename <- paste( szBaseName, nQuadType, 'N', 
                            formatC( nNodes, format='d', flag='0', width=5 ),
                            '.mat', sep='' )
    szVarPath     <- paste( szBaseDir, szSeedDir, szVarFilename, sep='/' )
  
    mSeedData <- readMat( szVarPath )

    # Store the data
    vIxRow <- seq( (ixSeed-1) * nStart + 1, ixSeed * nStart )

    mStdErr[ vIxRow, vIxCol ] <- t( mSeedData$mStdErr[ vIxCol, 1 : nStart ] )
  }

  return( mStdErr ) 
  
}


#-----------------------------------------------------------
# Create nice LaTeX table for point estimates and standard errors 
# for all data sets

MakePointEstStdErr <- function( szName, szCap, mStdErr )
{

  # Setup
  
  nDigits <- 4    # How many digits
  nWidth  <- 6    # Field width

  szIn    <- paste( 'PointEst-', szName, '.txt', sep='' )
  szOut   <- paste( 'PointEst-SE-', szName, '.inc', sep='' )
  szLabel <- paste( 'Tab.PointEst.SE', szName, sep='' )
#  szCap   <- paste( 'Point Estimates :', szName )


  #------------------
  # Get data
  info.df <- MakePretty( szName, szCap )

  df <- read.table( file = szIn )
  
  # Get rid of negative sign on scale of random coefficients
  #sol.df <- cbind( df[,1], info.df[,3], df[,2:11] )
  sol.df <- cbind( df[,1], info.df[,3], df[,2:6], abs( df[,7:11] ) )


  nSolRow <- nrow( sol.df )
  nSERow  <- nrow( mStdErr )

  nSolCol <- ncol( sol.df )
  nSECol  <- ncol( mStdErr )

  if( ( nSolRow != nSERow ) || ( nSolCol != ( nSECol + 2 ) ) ) {
    stop('Non-conformable point estimates and standard errors')
  }

  #------------------
  # Prepare standard errors

  # Format and add parentheses
  mSEFmt <- matrix (
              paste ( 
                '(', 
                formatC( mStdErr, format='G', digits=nDigits, width=nWidth, flag='#' ), 
#                format( mStdErr, digits=nDigits, width=nDigits, drop0trailing=FALSE ), 
                ')', sep='' 
              ),
              nrow=nSERow
            )

  # Pad standard error matrix with junk for columns 1 & 2
  mPad <- cbind( rep('',nSERow), rep('',nSERow))

  mSEFmt  <- cbind( mPad, mSEFmt )

  # Format Point Estimates
  mPointEstFmt <- cbind (
#                    formatC( sol.df[,1:2], format='d', digits=2, width=4, flag='#' ),
#                    formatC( sol.df[,3:nSolCol], format='G', digits=nDigits, width=nWidth, flag='#' )
                    sol.df[,1:2],
                    format( sol.df[,3:nSolCol], digits=nDigits, width=nWidth ) 
                  )

  # Put Point Estimates and Standard Errors together
  nTotRow <- 2 * nSolRow
  
  mAll <- matrix( nrow=nTotRow, ncol=nSolCol )
  mAll[ seq(1,nTotRow,2), ] <- as.matrix( mPointEstFmt )
  mAll[ seq(2,nTotRow,2), ] <- mSEFmt


  #------------------
  # Create xtable and print LaTeX

  vNames <- vector( length = nTotRow )
  vNames[ seq(1,nTotRow,2) ] <- rownames( sol.df )
  vNames[ seq(2,nTotRow,2) ] <- "s.e."

  colnames( mAll ) <- c( 'Data', 'INFORM', "$\\theta_{11}$", "$\\theta_{12}$", 
                            "$\\theta_{13}$", "$\\theta_{14}$",
                            "$\\theta_{15}$", "$\\theta_{21}$", 
                            "$\\theta_{22}$", "$\\theta_{23}$",
                            "$\\theta_{24}$", "$\\theta_{25}$" )

  xt <- xtable( mAll, caption=szCap, label=szLabel, digits=4 )
#                display=c( 'd', 'd', 'd', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G' ) )
  
  print ( 
    xt,
    file=szOut, type='latex', math.style.negative=T,
#    sanitize.rownames.function=function(x){ return( unlist(vNames)) },
    sanitize.text.function=function(x){x},
    hline.after=c( -1, 0, 2*seq(5,25,5) ), 
    include.rownames=FALSE, floating=T, caption.placement='bottom',
    floating.environment='sidewaystable',
    size='small'
  )

}


#-----------------------------------------------------------
# Create nice LaTeX table for point estimates and standard
# errors for all data sets

MakePointEstStdErrSmall <- function( szName, szCap, nSet, mStdErr )
{

  # Setup
  
  nDigits <- 4    # How many digits
  nWidth  <- 6    # Field width
  nStarts <- 5
  vIx     <- seq( nStarts * ( nSet - 1 ) + 1, nStarts * nSet, 1 )  # Access correct data slice 

  szIn    <- paste( 'PointEst-', szName, '.txt', sep='' )
  szOut1   <- paste( 'PointEst-SE-', szName, '-Short', nSet, '.1', '.inc', sep='' )
  szOut2   <- paste( 'PointEst-SE-', szName, '-Short', nSet, '.2', '.inc', sep='' )
  
  szLabel <- paste( 'Tab.PointEst.SE', szName, sep='' )


  #------------------
  # Get data
  info.df <- MakePretty( szName, szCap )

  df <- read.table( file = szIn )
  
  # Get rid of negative sign on scale of random coefficients
  #sol.df <- cbind( df[,1], info.df[,3], df[,2:11] )
  sol.df <- cbind( df[vIx,1], info.df[vIx,3], df[vIx,2:6], abs( df[vIx,7:11] ) )

  mStdErr <- mStdErr[vIx,]

  nSolRow <- nrow( sol.df )
  nSERow  <- nrow( mStdErr )

  nSolCol <- ncol( sol.df )
  nSECol  <- ncol( mStdErr )

  if( ( nSolRow != nSERow ) || ( nSolCol != ( nSECol + 2 ) ) ) {
    stop('Non-conformable point estimates and standard errors')
  }

  #------------------
  # Prepare standard errors

  # Format and add parentheses
  mSEFmt <- matrix (
              paste ( 
                '(', 
                formatC( mStdErr, format='G', digits=nDigits, width=nWidth, flag='#' ), 
#                format( mStdErr, digits=nDigits, width=nDigits, drop0trailing=FALSE ), 
                ')', sep='' 
              ),
              nrow=nSERow
            )

  # Pad standard error matrix with junk for columns 1 & 2
  mPad <- cbind( rep('',nSERow), rep('',nSERow))

  mSEFmt  <- cbind( mPad, mSEFmt )

  # Format Point Estimates
  mPointEstFmt <- cbind (
#                    formatC( sol.df[,1:2], format='d', digits=2, width=4, flag='#' ),
#                    formatC( sol.df[,3:nSolCol], format='G', digits=nDigits, width=nWidth, flag='#' )
                    sol.df[,1:2],
                    format( sol.df[,3:nSolCol], digits=nDigits, width=nWidth ) 
                  )

  # Put Point Estimates and Standard Errors together
  nTotRow <- 2 * nSolRow
  
  mAll <- matrix( nrow=nTotRow, ncol=nSolCol )
  mAll[ seq(1,nTotRow,2), ] <- as.matrix( mPointEstFmt )
  mAll[ seq(2,nTotRow,2), ] <- mSEFmt


  #------------------
  # Create xtable and print LaTeX

  vNames <- vector( length = nTotRow )
  vNames[ seq(1,nTotRow,2) ] <- rownames( sol.df )
  vNames[ seq(2,nTotRow,2) ] <- "s.e."

  colnames( mAll ) <- c( 'Data', 'INFORM', "$\\theta_{11}$", "$\\theta_{12}$", 
                            "$\\theta_{13}$", "$\\theta_{14}$",
                            "$\\theta_{15}$", "$\\theta_{21}$", 
                            "$\\theta_{22}$", "$\\theta_{23}$",
                            "$\\theta_{24}$", "$\\theta_{25}$" )

#  xt <- xtable( mAll[,c(1,2,3:7)], caption=szCap, label=szLabel, digits=4 )
  xt <- xtable( mAll[,c(3:7)], caption=szCap, label=szLabel, digits=2 )
  
  print( xt, file=szOut1, type='latex', math.style.negative=T, 
        sanitize.text.function=function(x){ x }, 
        include.rownames=F, floating=T, caption.placement='bottom' )


#  xt <- xtable( mAll[,c(1,2,8:12)], caption=szCap, label=szLabel, digits=4, 
#                display=c( 'd', 'd', 'd', 'G', 'G', 'G', 'G', 'G' ) )
  xt <- xtable( mAll[,c(8:12)], caption=szCap, label=szLabel, digits=2, 
                display=c( 'd', 'G', 'G', 'G', 'G', 'G' ) )
  
  print( xt, file=szOut2, type='latex', math.style.negative=T, 
        sanitize.text.function=function(x){ x }, 
        include.rownames=F, floating=T, caption.placement='bottom' )

  # Kludge - make a one line version if processing GH7-5Good
  if( szName == 'GH7-5Good' ){

    szOut3  <- paste( 'PointEst-', szName, '-GH7Cmp.', nSet, '.1', '.inc', sep='' )
    szOut4  <- paste( 'PointEst-', szName, '-GH7Cmp.', nSet, '.2', '.inc', sep='' )

#    xt <- xtable( mAll[c(9,10),c(1,2,3:7)], caption=szCap, digits=4, 
#                  display=c( 'd', 'd', 'd', 'G', 'G', 'G', 'G', 'G' ) )
    xt <- xtable( mAll[c(9,10),c(3:7)], caption=szCap, digits=2, 
                  display=c( 'd', 'G', 'G', 'G', 'G', 'G' ) )
  
    print( xt, file=szOut3, type='latex', math.style.negative=T, 
          include.rownames=F, floating=T, caption.placement='bottom',
          sanitize.text.function=function(x){ x } )
  
#    xt <- xtable( mAll[c(9,10),c(1,2,8:12)], caption=szCap, digits=4, 
#                  display=c( 'd', 'd', 'd', 'G', 'G', 'G', 'G', 'G' ) )
    xt <- xtable( mAll[c(9,10),c(8:12)], caption=szCap, digits=4, 
                  display=c( 'd', 'G', 'G', 'G', 'G', 'G' ) )
    
    print( xt, file=szOut4, type='latex', math.style.negative=T, 
          include.rownames=F, floating=T, caption.placement='bottom',
          sanitize.text.function=function(x){ x } )
  }
}


#-------------------------------------------------------------------------------
# MakeTables - generates all tables

MakeTables <- function()
{

  tmp <- MakePretty( 'pMC', 'MPEC Results: pMC with $R = 1,000$ draws' )
  tmp <- MakePrettySmall( 'pMC', 'MPEC Results: pMC with $R = 1,000$ draws', 1 )
  MakePointEst( 'pMC', 'Point Estimates: pMC with $R = 1,000$ draws' )
  MakePointEstSmall( 'pMC', 'Point Estimates: pMC with $R = 1,000$ draws', 1 )

  tmp <- MakePretty( 'pMC-R10000', 'MPEC Results: pMC with $R = 10,000$ draws' )
  tmp <- MakePrettySmall( 'pMC-R10000', 'MPEC Results: pMC with $R = 10,000$ draws', 1 )
  MakePointEst( 'pMC-R10000', 'Point Estimates: pMC with $R = 10,000$ draws' )
  MakePointEstSmall( 'pMC-R10000', 'Point Estimates: pMC with $R = 10,000$ draws', 1 )


  tmp <- MakePrettySmall( 'pMC-R10000', 'MPEC Results: pMC with $R = 10,000$ draws', 5 )
  MakePointEstSmall( 'pMC-R10000', 'Point Estimates: pMC with $R = 10,000$ draws', 5 )

  tmp <- MakePretty( 'SGI', 'MPEC Results: SGI with $993$ nodes (exact for degree $\\leq 11$).' )
  tmp <- MakePrettySmall( 'SGI', 'MPEC Results: SGI with $993$ nodes (exact for degree $\\leq 11$).', 1 )
  MakePointEst( 'SGI', 'Point Estimates: SGI with $993$ nodes (exact for degree $\\leq 11$).' )
  MakePointEstSmall( 'SGI', 'Point Estimates: SGI with $993$ nodes (exact for degree $\\leq 11$).', 1 )

  tmp <- MakePretty( 'GH5', 'MPEC Results: Gauss-Hermite with $5^5$ nodes.' )
  tmp <- MakePrettySmall( 'GH5', 'MPEC Results: Gauss-Hermite with $5^5$ nodes.', 1 )
  MakePointEst( 'GH5', 'Point Estimates: Gauss-Hermite with $5^5$ nodes.' )
  MakePointEstSmall( 'GH5', 'Point Estimates: Gauss-Hermite with $5^5$ nodes.', 1 )

  #----------------------------------------------------------
  # Handle new runs
  tmp <- MakePretty( 'Mono5Good', 
                     'Point Estimates: Monomial with first $5$ good starts and $983$ nodes (exact for degree $\\leq 11$).' )
  tmp <- MakePrettySmall( 'Mono5Good', 
                     'Point Estimates: Monomial with first $5$ good starts and $983$ nodes (exact for degree $\\leq 11$).', 5 )
  MakePointEst( 'Mono5Good', 
                'Point Estimates: Monomial with first $5$ good starts and $983$ nodes (exact for degree $\\leq 11$).' )
  MakePointEstSmall( 'Mono5Good', 
                'Point Estimates: Monomial with first $5$ good starts and $983$ nodes (exact for degree $\\leq 11$).' , 5 )

  tmp <- MakePretty( 'pMC5Good-R01000', 
        'Point Estimates: pMC with first $5$ good starts and $R = 1,000$ draws.' )
  tmp <- MakePrettySmall( 'pMC5Good-R01000', 
        'Point Estimates: pMC with first $5$ good starts and $R = 1,000$ draws.', 5 )
  MakePointEst( 'pMC5Good-R01000', 
        'Point Estimates: pMC with first $5$ good starts and $R = 1,000$ draws.' )
  MakePointEstSmall( 'pMC5Good-R01000', 
        'Point Estimates: pMC with first $5$ good starts and $R = 1,000$ draws.', 5 )

  tmp <- MakePretty( 'SGI5Good', 
        'Point Estimates: SGI with first $5$ good starts and $993$ nodes (exact for degree $\\leq 11$).' ) 
  tmp <- MakePrettySmall( 'SGI5Good', 
        'Point Estimates: SGI with first $5$ good starts and $993$ nodes (exact for degree $\\leq 11$).', 5 ) 
  MakePointEst( 'SGI5Good', 
        'Point Estimates: SGI with first $5$ good starts and $993$ nodes (exact for degree $\\leq 11$).' ) 
  MakePointEstSmall( 'SGI5Good', 
        'Point Estimates: SGI with first $5$ good starts and $993$ nodes (exact for degree $\\leq 11$).', 5 ) 

  tmp <- MakePretty( 'GH7-5Good', 
        'Point Estimates: Gauss-Hermite with first $5$ good starts and $7^5$ nodes.' ) 
  tmp <- MakePrettySmall( 'GH7-5Good', 
        'Point Estimates: Gauss-Hermite with first $5$ good starts and $7^5$ nodes.', 5 ) 
  MakePointEst( 'GH7-5Good', 
        'Point Estimates: Gauss-Hermite with first $5$ good starts and $7^5$ nodes.' ) 
  MakePointEstSmall( 'GH7-5Good', 
        'Point Estimates: Gauss-Hermite with first $5$ good starts and $7^5$ nodes.', 5 ) 

}

MakeNDTables <- function()
{
  tmp <- MakePretty( 'ND', 
             'MPEC Results: Same starting value, multiple draws of pMC nodes ($R=1000$).' )
  tmp <- MakePrettySmall( 'ND', 
             'MPEC Results: Same starting value, multiple draws of pMC nodes ($R=1000$).', 1 )
  MakePointEst( 'ND', 'Point Estimates: Same starting value, multiple draws of pMC nodes ($R=1000$).' )
  MakePointEstSmall( 'ND', 
      'Point Estimates: Same starting value, multiple draws of nodes ($R=1000$).', 1 )
}

MakeStdErrTables <- function()
{

 # Monomial
 mStdErr <- LoadStdErr('~/sbox/quad/src/DataJ25T50/', 2, 983)
 MakePointEstStdErr( 'Monomial', 'Point Estimates: Stroud Monomial Rule 11-1.', mStdErr )
 MakePointEstStdErrSmall( 'Monomial', 'Point Estimates: Stroud Monomial Rule 11-1.', 5, 
     mStdErr )

#  # SGI
#  mStdErr <- LoadStdErr('~/sbox/quad/src/DataJ25T50/', 4, 993)
#  MakePointEstStdErr( 'SGI', 
#      'Point Estimates: SGI with 993 nodes (exact for degree $\\leq 11$)', mStdErr )
#  MakePointEstStdErrSmall( 'SGI',
#      'Point Estimates: SGI with 993 nodes (exact for degree $\\leq 11$)', 1, mStdErr )
#
#  # pMC with R=1000
#  mStdErr <- LoadStdErr('~/sbox/quad/src/DataJ25T50/', 0, 1000)
#  MakePointEstStdErr( 'pMC-R01000', 'Point Estimates: pMC with $R = 1,000$ draws', mStdErr )
#  MakePointEstStdErrSmall( 'pMC-R01000', 'Point Estimates: pMC with $R = 1,000$ draws', 
#                           1, mStdErr )
#
#  # GH5
#  mStdErr <- LoadStdErr('~/sbox/quad/src/DataJ25T50/', 1, 3125)
#  MakePointEstStdErr( 'GH5', 'Point Estimates: Gauss-Hermite with $5^5$ nodes', mStdErr )
#  MakePointEstStdErrSmall( 'GH5', 'Point Estimates: Gauss-Hermite with $5^5$ nodes', 
#                           1, mStdErr )

  # MonoBox
  mStdErr <- LoadStdErr('~/sbox/quad/src/DataJ25T50/', 2, 983)
  MakePointEstStdErr( 'MonoBox', 'Point Estimates: Monomial with Tight Box and NaN -> Inf', 
                          mStdErr )
  MakePointEstStdErrSmall( 'MonoBox', 'Point Estimates: Monomial with Tight Box and NaN -> Inf', 
                           5, mStdErr )

 # pMC with R=10,000
 mStdErr <- LoadStdErr('~/sbox/quad/src/DataJ25T50/', 0, 10000)
 MakePointEstStdErr( 'pMC-R10000', 
                     'Point Estimates: pMC with $R = 10,000$ draws.', mStdErr )
 MakePointEstStdErrSmall( 'pMC-R10000', 'Point Estimates: pMC with $R = 10,000$ draws.', 
                          5, mStdErr )

 # Monomial - drop starts until 5 converge
 mStdErr <- LoadStdErr('~/sbox/quad/src/DataJ25T50/', 2, 983)
 MakePointEstStdErr( 'Mono5Good', 
                     'Point Estimates: Monomial with first $5$ good starts and $983$ nodes (exact for degree $\\leq 11$).', 
                     mStdErr )
 MakePointEstStdErrSmall( 'Mono5Good', 
                     'Point Estimates: Monomial with first $5$ good starts and $983$ nodes (exact for degree $\\leq 11$).' ,
                     5, mStdErr )


 # pMC-R01000 - drop starts until 5 converge
 mStdErr <- LoadStdErr('~/sbox/quad/src/DataJ25T50/', 0, 1000)
 MakePointEstStdErr( 'pMC5Good-R01000', 
   'Point Estimates: pMC with first $5$ good starts and $R = 1,000$ draws.', mStdErr )
 MakePointEstStdErrSmall( 'pMC5Good-R01000', 
   'Point Estimates: pMC with first $5$ good starts and $R = 1,000$ draws.', 
   5, mStdErr )


 # SGI
 mStdErr <- LoadStdErr('~/sbox/quad/src/DataJ25T50/', 4, 993)
 MakePointEstStdErr( 'SGI5Good', 
     'Point Estimates: SGI with first $5$ good starts and $993$ nodes (exact for degree $\\leq 11$).', 
     mStdErr )
 MakePointEstStdErrSmall( 'SGI5Good',
     'Point Estimates: SGI with first $5$ good starts and $993$ nodes (exact for degree $\\leq 11$).', 
     5, mStdErr )


#  # qMC-R01000 - drop starts until 5 converge
#  mStdErr <- LoadStdErr('~/sbox/quad/src/DataJ25T50/', 0, 1000)
#  MakePointEstStdErr( 'qMC5Good-R01000', 
#    'Point Estimates: qMC with first $5$ good starts and $R = 1,000$ draws', mStdErr )
#  MakePointEstStdErrSmall( 'qMC5Good-R01000', 
#    'Point Estimates: qMC with first $5$ good starts and $R = 1,000$ draws', 
#    1, mStdErr )


  # GH57 - 5Good
  mStdErr <- LoadStdErr('~/sbox/quad/src/DataJ25T50/', 1, 16807)
  MakePointEstStdErr( 'GH7-5Good', 
    'Point Estimates: Gauss-Hermite with first $5$ good starts and $7^5$ nodes.', 
    mStdErr )
  MakePointEstStdErrSmall( 'GH7-5Good', 
    'Point Estimates: Gauss-Hermite with first $5$ good starts and $7^5$ nodes.', 
    1, mStdErr )
  MakePointEstStdErrSmall( 'GH7-5Good', 
    'Point Estimates: Gauss-Hermite with first $5$ good starts and $7^5$ nodes.', 
    5, mStdErr )

}
