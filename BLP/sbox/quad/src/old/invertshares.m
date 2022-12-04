function [ delta, ii ] = invertshares(theta2)

% INVERTSHARES
% GMM objective function for the random coefficients Logit estimated via
% the NFP approach proposed by BLP.
%
% source: Dube, Fox and Su (2009)
% Code Revised: March 2009

global x expmeanval tol_inner share v

ii = 0;
norm = 1;
expmeanval0 = expmeanval;

mMu   = x*diag(theta2)*v ;  % Individual taste-specific utility
expmu = exp( mMu );      % exponentiated deviations from mean utilities
                                    % XXX replace v with quadrature nodes

while norm > tol_inner && ii < 2500,
  vProb = ind_shnorm( expmeanval0, expmu ) ; 
  expmeanval1 = expmeanval0 .* share ./ vProb ;
% %     vDelta      = log( expmeanval0 ) ;    % slower but needed for stability
% %     expmeanval1 = expmeanval0.*share./ind_shnorm_stable( vDelta, mMu ); 
    
  t = abs(expmeanval1-expmeanval0);
  norm = max(t);
  expmeanval0 = expmeanval1;
  ii = ii + 1;
end;

if ii >= 2500
  fprintf( 2, 'invertshares: mapping failed to converge!\n' ) ;
% else
%   fprintf( 1, '\tinverstshares: converged in %d iterations.\n', ii ) ;
end

if max(isnan(expmeanval0))<1, expmeanval = expmeanval0; end
delta = log(expmeanval);               % the mean utilities

Validate( vProb ) ;
Validate( delta ) ;