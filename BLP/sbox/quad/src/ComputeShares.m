% ComputeShares - computes market shares
%
% modification history
% --------------------
% 26oct2009 bss hacked from JP's code.
%

function [ vShares ] = ComputeShares( x0 )

  global x K v numProdsTotal

  % Note:  we don't need theta1 because all of that and more is in delta
  % because of the GMM setup.
  theta2 = x0( K + 2 : 2 * K + 2, 1 ) ;                  % st. deviation of tastes
  delta  = x0( 2 * K + 3 : 2 * K + 2 + numProdsTotal, 1 ) ;         % mean utilities

  expmu = exp(x*diag(theta2)*v);      % exponentiated deviations from mean utilities
  expmeanval = exp(delta);
  vShares    = ind_shnorm( expmeanval, expmu ) ; % constraints on predicted shares
