% MapDiag - Diagnostics about Berry's mapping
%
% Compute the Jacobian for Berry's approximate mapping to determine if it
% is a contraction.  
%
% Things to note:
%   - the largest eigenvalue should be less than 1 and determines the speed
%     of convergence.
%   - the smallest eigenvalue(s) determine if the Jacobian is singular in
%     which case there are problems....
%
% To use:
%   1) Run CompareShares.m
%   2) Run this M-file.
%
% modification history
% --------------------
% 30nov2009 bss written.
%

%% Choose Quadrature Method for Jacobian

bSim = 0 ;      % True to use simulation for quadrature
if bSim 
%   v_sim = randn( length( betatrue ), N_DRAWS_MC_INT ) ; 
  v  = v_sim ;
  nn = size( v, 2 ) ;
  Q_WEIGHTS = ones( nn, 1 ) / nn ;
  
  
  oo = ones( 1, nn ) ; 
  szTitle = sprintf( 'Histogram of log10 of Magnitude of Eigenvalues: Simulation N = %d', N_MC_INT ) ;
else
%   GHQuadInit( ) ;
  nStroudRule = 2 ;
  MonoQuadInit( nStroudRule ) ;
  v  = Q_NODES ;
  nn = size( v, 2 ) ;
  
  oo = ones( 1, nn ) ; 
  szTitle = sprintf( 'Histogram of log10 of Magnitude of Eigenvalues: Stroud Monomial Rule %d', nStroudRule ) ;
end
  

%% Compute Jacobian
%
% Compute jacobian at same xi -- that calculated via simulation to compare
% quality of jacobians under different integration rules.

if bSim
  mJ = CSDJacobian( @(x) ComputeMap( x, xOpt ), mResidSim( :, 1 ), 1e-5 ) ;
else
  mJ = CSDJacobian( @(x) ComputeMap( x, xOpt ), mResidSim( :, 1 ), 1e-5 ) ;
  % mJ = CSDJacobian( @(x)ComputeMap( x, xOpt ), mResidStroud(:,1), 1e-5 ) ;
  % mJ = CSDJacobian( @(x)ComputeMap( x, xOpt ), mResidStroud(:,2), 1e-5 )
end


%% Diagnostics

disp( '***Smallest Eigenvalues' ) ;
eigs( mJ, 10, 'sm' )    % 50 smallest eigenvalues

disp( '***Largest Eigenvalues' ) ;
eigs( mJ, 5)            % 5 largest eigenvalues

mJFull = full( mJ ) ;
vEigs = eig( mJFull ) ;

vEdges = -12 : 0.5 : 0 ;
hist( log10( abs( vEigs ) ), vEdges ) ;
title( szTitle ) ;
xlabel( 'log10( abs( \lambda ) )' ) ;
ylabel( 'Frequency' ) ;

disp( '***Condition Number' ) ;
cond( mJFull )

eigMax = eigs( mJ, 1, 'lm' ) ;
eigMin = eigs( mJ, 1, 'sm' ) ;

fprintf( 1, 'Condition Number: %10.5g\n', eigMax / eigMin ) ;

