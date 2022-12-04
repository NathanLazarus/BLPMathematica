% TestMonoQuadInit - test quadrature rule with some monomials which should integrate exactly.
%
% modification history
% --------------------
% 15jan2010 bss written.
%

%% Setup

global Q_NODES ;
global Q_WEIGHTS ;
global N_DIM ;

ixRule = 1 ;   % Can be 1 or 2 for the L or R column in Stroud.

MonoQuadInit( ixRule ) ;

% Rescale weights to correct values instead of those for a kernel
% normailzed for a standard normal

Q_WEIGHTS = Q_WEIGHTS * (  pi ^ ( N_DIM / 2 ) ) ;


%% Now let's integrate some monomials
%
% Note:  the integral from -infty to +infty of x^n * exp( -a *x^2 ) is 
% ( n - 1 ) !! / (2 * a ) ^ n/2 * sqrt( pi / a )

% 1. Monomial is 1

fprintf( 1, '---> Diff for monomial( 0, 0, 0, 0, 0 ) is : %10.4e\n', sum( Q_WEIGHTS ) - (  pi ^ ( N_DIM / 2 ) ) ) ;

% 2. Monomial is 2

nDeg = N_DIM ;

% Anonymous function to compute correct answer in one dimension
fCor = @( nDegree )  prod( 1 : 2 : ( nDegree - 1 ) ) / ( 2 ^ ( nDegree / 2 ) ) * sqrt( pi ) ;

dwCor = fCor( 4  ) * fCor( 0 ) ^ 4 ;

fMono  = @( x ) ( x( 1, : ) .^ 4 ) * Q_WEIGHTS ;   % Need to rescale based on number 
                                                                        % of non-degenerate monomial terms
dwMono = fMono( Q_NODES ) ;

fprintf( 1, '---> Diff for monomial( 4, 0, 0, 0, 0 ) is : %10.4e\n', dwCor - dwMono ) ; 


%% 3 Monomial is 2 2 

% % fCor = @( nDegree )  prod( 1 : 2 : ( nDegree - 1 ) ) / ( 2 ^ ( nDegree / 2 ) ) * sqrt( pi ) ;

dwCor = fCor( 2 ) ^ 2 * fCor( 0 ) ^ 3 ;

fMono  = @( x ) prod( x(1:2,:) .^ 2 ) * Q_WEIGHTS ;   % Need to rescale based on number 
                                                                        % of non-degenerate monomial terms
dwMono = fMono( Q_NODES ) ;

fprintf( 1, '---> Diff for monomial( 2, 2, 0, 0, 0 ) is : %10.4e\n', dwCor - dwMono ) ; 


%% 4 Monomial is 6 2 2 

% % fCor = @( nDegree )  prod( 1 : 2 : ( nDegree - 1 ) ) / ( 2 ^ ( nDegree / 2 ) ) * sqrt( pi ) ;

dwCor = fCor( 6 ) * fCor( 2 ) ^2 * fCor( 0 ) ^ 2 ;

fMono  = @( x )  ( ( x(1,:) .^ 6 ) .* ( x(2,:) .^ 2 ) .* ( x(3,:) .^ 2 ) ) * Q_WEIGHTS ;   % Need to rescale based on number 
                                                                        % of non-degenerate monomial terms
dwMono = fMono( Q_NODES ) ;

fprintf( 1, '---> Diff for monomial( 6, 2, 2, 0, 0 ) is : %10.4e\n', dwCor - dwMono ) ; 



%% 5 Monomial is 10 2 2  (This should not be exact) 

% % fCor = @( nDegree )  prod( 1 : 2 : ( nDegree - 1 ) ) / ( 2 ^ ( nDegree / 2 ) ) * sqrt( pi ) ;

dwCor = fCor( 10 ) * fCor( 2 ) ^2 * fCor( 0 ) ^ 2 ;

fMono  = @( x )  ( ( x(1,:) .^ 10 ) .* ( x(2,:) .^ 2 ) .* ( x(3,:) .^ 2 ) ) * Q_WEIGHTS ;   % Need to rescale based on number 
                                                                        % of non-degenerate monomial terms
dwMono = fMono( Q_NODES ) ;

fprintf( 1, '---> Diff for monomial( 10, 2, 2, 0, 0 ) is : %10.4e\n', dwCor - dwMono ) ; 



%% 6 Monomial is 6 2 1  (This should not be exact) 

% % fCor = @( nDegree )  prod( 1 : 2 : ( nDegree - 1 ) ) / ( 2 ^ ( nDegree / 2 ) ) * sqrt( pi ) ;

dwCor = 0 ; 

fMono  = @( x )  ( ( x(1,:) .^ 8 ) .* ( x(2,:) .^ 2 ) .* ( x(3,:) .^ 1 ) ) * Q_WEIGHTS ;   % Need to rescale based on number 
                                                                        % of non-degenerate monomial terms
dwMono = fMono( Q_NODES ) ;

fprintf( 1, '---> Diff for monomial( 8, 2, 1, 0, 0 ) is : %10.4e\n', dwCor - dwMono ) ; 


%% Any monomial vP

vP = [ 2 2 4 2 0 ] ;
if length( vP ) ~= N_DIM
  error( 'Vector of powers must match number of degrees' ) ;
end

dwCor = 1 ;
for ix = 1 : N_DIM
  if 0 == mod( vP( ix ), 2 )
    dwCor = dwCor * fCor( vP( ix ) ) ;
  else
    dwCor = 0 ;
  end
end

fMono  = @( x )  ( ( x(1,:) .^ vP( 1 ) ) .* ( x(2,:) .^ vP( 2 ) ) .* ( x(3,:) .^ vP( 3 ) )        ...
                                         .* ( x(4,:) .^ vP( 4 ) ) .* ( x(5,:) .^ vP( 5 ) ) ) * Q_WEIGHTS ;  


dwMono = fMono( Q_NODES ) ;

fprintf( 1, '---> Diff for monomial( %d %d %d %d %d ) is : %10.4e\n', vP, dwCor - dwMono ) ;                                        