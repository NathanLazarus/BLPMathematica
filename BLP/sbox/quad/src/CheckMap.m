% CheckMap - checks that Jacobian of Berry's map is a contraction
%
% Checks that Berry's mapping satisfies the conditions in BLP (1995) which
% are sufficient to be a contraction mapping.
%
% modification history
% --------------------
% 04jun2010 bss written.
%

function [ ] = CheckMap( x0 )

  global T prods ;
  
  mFullJac = MapJacobian( x0 ) ;

  if any( mFullJac(:) < 0 )
    nBad = sum( mFullJac(:) < 0 ) ;
    fprintf( 2, 'WARNING: mapping jacobian < 0!\n' ) ;
    fprintf( 2, '---> %d of %d violations.\n', nBad, T*prods*prods ) ;
  else
    fprintf( 1, 'mapping jacobian nonnegative.\n' ) ;
  end

  vJacSum = sum( mFullJac, 2 ) ;
  if any( vJacSum >= 1 )
    nBad = sum( vJacSum >= 1 ) ;
    fprintf( 2, 'WARNING: sum of mapping jacobian > 1!\n' ) ;
    fprintf( 2, '---> %d of %d violations.\n', nBad, T*prods ) ;
  else
    fprintf( 2, 'sum of mapping jacobian < 1.\n' ) ;
  end
