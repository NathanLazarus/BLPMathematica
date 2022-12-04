function [ delta, mInfo ] = MapMapping(theta2)

% MapMapping
%
% Records diagnositic information from Berry's Mapping to determine rate of
% convergence.
%
% source: Dube, Fox and Su (2009)
% Code Revised: March 2009


%% Setup

global x expmeanval tol_inner share v

nMaxIter    = 2500 ;
ii          = 1;
dwStopCrit  = 1;
expmeanval0 = expmeanval;

expmu = exp(x*diag(theta2)*v);      % exponentiated deviations from mean utilities
                                    % XXX replace v with quadrature nodes
                                    
%% Configure data logging for mapping

dwStep = 25 ;
vStop  = 200 : dwStep : 500 ;        % Which iterations to save
nStop  = length( vStop ) ;
nElem  = length( expmeanval ) ;

mInfo.mStats = zeros( nMaxIter, 10 ) ;   % Stats are ( L-1, L-2, L-Inf norms and 25%, 50%, 75% quantiles of change in shock )
mInfo.mData  = zeros( nElem, nStop ) ;
mInfo.dwStep = dwStep ;

%% Do the mapping until MaxIter or tolerance is satisfied

while dwStopCrit > tol_inner && ii < nMaxIter,
    expmeanval1 = expmeanval0.*share./ind_shnorm(expmeanval0,expmu); 
    t = abs(expmeanval1-expmeanval0);
    dwStopCrit = max(t);
    
    %----------------------------------
    % Compute Stats
    
    vResid = expmeanval1 - expmeanval0 ;
    mInfo.mStats( ii, 1 ) = norm( vResid, 1 ) ;
    mInfo.mStats( ii, 2 ) = norm( vResid, 2 ) ;
    mInfo.mStats( ii, 3 ) = norm( vResid, Inf ) ;
    
    mInfo.mStats( ii, 4:end ) = quantile( abs( vResid ), [ 0.01 ; 0.05 ; 0.10 ; 0.15 ; 0.25 ; 0.5 ; 0.75 ] ) ;
    
    ix = find( vStop == ii ) ;
    if ~isempty( ix )
      mInfo.mData( :, ix ) = expmeanval1 ;
    end
    
    
    %----------------------------------
    
    expmeanval0 = expmeanval1;
    ii = ii + 1;
end;

if ii >= nMaxIter
  fprintf( 2, 'invertshares: mapping failed to converge!\n' ) ;
% else
%   fprintf( 1, '\tinverstshares: converged in %d iterations.\n', ii ) ;
end

mInfo.ii = ii ;

if max(isnan(expmeanval0))<1, expmeanval = expmeanval0; end
delta = log(expmeanval);               % the mean utilities


% Display approximation for beta
for ixCol = 3 : nStop
  beta = ( norm( mInfo.mData( :, ixCol - 2 ) - mInfo.mData( :, ixCol - 1 ), 2 )        ...
        / norm( mInfo.mData( :, ixCol - 1 ) - mInfo.mData( :, ixCol ), 2 ))^( 1 / ( vStop( 2 ) - vStop( 1 ) ) )
end

