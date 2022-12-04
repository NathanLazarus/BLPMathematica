# MakeStdErrTable - Prints a xtable with standard errors and point estimates
#
# modification history
# --------------------
# 24aug2010 bss written.
#

MakeStdErrTable <- function( szOutFile, szCaption, szLabel, 
                             mPointEst, mStdErr, nDigits = 4 )
{

  #--------------------
  # Make sure tables are conformable
  nRow <- nrow( mPointEst )
  nCol <- ncol( mPointEst )

  if( nRow != nrow( mStdErr ) || nCol != ncol( mStdErr ) ) {
    stop( 'mPointEst and mStdErr are non-conformable!' )
  }

  #--------------------
  # Get basics
  vColNames <- colnames( mPointEst )
  vRowNames <- rownames( mPointEst )

  nTotRows <- 2 * nRow


  nTotCol <- nCol
  mAll <- matrix( nrow=nTotRows,ncol=nCol )
  mAll[ seq(1,nTotRows,2), ] <- round( mPointEst, digits=nDigits )
  mAll[ seq(2,nTotRows,2), ] <- matrix(
                  paste( '(', round( mStdErr, digits=nDigits ), ')', sep='' ),
                  nrow=nrow(mStdErr)
                                  )

  #mAll <- data.frame( mAll )
  
  colnames( mAll ) <- vColNames

  vNames <- vector( length = nTotRows )
  vNames[ seq(1,nTotRows,2) ] <- vRowNames
  vNames[ seq(2,nTotRows,2) ] <- "s.e."


  #--------------------
  # Print the xtable

  xt <- xtable( mAll,
                caption=szCaption,
                label=szLabel,
                digits=nDigits )
  print( xt,
         szOutFile,
         type='latex', math.style.negative=T,
         sanitize.rownames.function=function(x){ return( unlist(vNames)) },
         sanitize.text.function=function(x){x},
         hline.after=c( -1, seq( 0, nRow, 2 ) ),
         include.rownames=T, floating=T, caption.placement='bottom',
         floating.environment='sidewaystable' )
 
 
  #--------------------
  # Return data if desired
  rownames( mAll ) <- vNames

  return( mAll )
}
