# Join lines together to avoid editing
#
# modification history
# --------------------
# 17may2010 bss written.
#

#NR % (1 + 3 * 5 ) == 1 { 
#                         nSeed = $1 
#                         nFound++ 
#                       }

{
  ixLine   = NR % (1 + 3 * 5 )     # First line with seed of data set
  ixRelPos = ( ixLine - 1 ) % 3 
  
  #print NR, ixLine, ( ixLine - 1 ) % 3 

  if( 1 == ixLine )
  {
    nSeed   = $1
  }
  else
  {

    if( 1 == ixRelPos )
    {
      f_k     = $1
    }
    else if( 2 == ixRelPos )
    {
      nExit   = $1
      nInform = $2
    }
    else if( 0 == ixRelPos || 0 == ixLine )
    {
      dwTime   = $1
      print nSeed, f_k, nExit, nInform, dwTime
    }

  }
}
