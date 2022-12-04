# Join lines together for point estimates to avoid editing
#
# modification history
# --------------------
# 18may2010 bss written.
#

#NR % (1 + 3 * 5 ) == 1 { 
#                         nSeed = $1 
#                         nFound++ 
#                       }

BEGIN {
  nStarts        = 5
  nLinesPerStart = 6
  nLinesPerSeed  = 1 + nStarts * nLinesPerStart
}

{
  ixLine   = NR % nLinesPerSeed     # First line with seed of data set
  ixRelPos = ( ixLine - 1 ) % nLinesPerStart 
  
  #print NR, ixLine, ( ixLine - 1 ) % 3 

  if( 1 == ixLine )
  {
    nSeed   = $1
  }
  else
  {

    if( 2 == ixRelPos )
    {
      line1 = $0
    }
    else if( 3 == ixRelPos )
    {
      line2 = $0
    }
    else if( 4 == ixRelPos )
    {
      line3 = $0
    }
    else if( 5 == ixRelPos )
    {
      line4 = $0
    }
    else if( 0 == ixRelPos || 0 == ixLine )
    {
      print nSeed, line1, line2, line3, line4
    }

  }
}
