% CSDJacobian - Complex step derivative Jacobian of a vector-valued function.
%
% Calculates the Jacobian of a vector-valued function using a complex step derivative.
% Note: functions must not use Matlab's ' operator which performs a
% transpose and complex complement.  Instead, use .' which only performs a
% transpose.
%
% CSD seems to work better with a smaller step size than finite difference.
%
% modification history
% --------------------
% 16jul2009 bss written.
%

function [ mCSDJac ] = CSDJacobian( func, x0, dwStep ) 

  nParams = length( x0 ) ;
  
  tmp   = func( x0 ) ;
  nFunc = length( tmp ) ;
  clear( 'tmp' ) ;
  
  %mCSDJac = spzeros( nFunc, nParams ) ;
  mCSDJac = zeros( nFunc, nParams ) ;
	% Calculate FD derivative by percentage
  %   fdTol = 0.1 ;  % step-size as a percentage of xFree.
  %   dx = fdTol * x0 ;
  
  if nargin < 3
    dx = 1e-5 ;
  else
    dx = dwStep ;
  end

  xPlus = x0 + 1i * dx ;

  for ix = 1 : nParams
    x1 = x0 ;
    x1( ix ) = xPlus( ix ) ;
    [ fval ] = func( x1 ) ;
    mCSDJac( :, ix ) = imag( fval / dx ) ;
  end
    
  