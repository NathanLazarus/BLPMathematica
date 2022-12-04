% ComputeMap - compute the approximate (contraction) mapping of Berry
% 
% Use parameter estimates and value of xi_jt unobserved product-market
% shock to compute one iteration of Berry's mapping.
%
% To compute the gradient with INTLAB, first compute the values of x0 at
% the optimum and the residual xi_jt.  Then initialize an INTLAB AD object:
%
%   g_vShock = gradientinit( vShock ) ;
%
% and compute the gradient as
%
%   ff = ComputeMap( g_vShock, x0 ) ;
%
% ff.dx is the Jacobian.
%
%
% modification history
% --------------------
% 23nov2009 bss written.
%

function [ vShock1 ] = ComputeMap( vShock, x0 )

  global share ;
  global x K v ;


  %% Marshall arguments
  theta1 = x0( 1 : K + 1, 1 ) ;
  theta2 = x0( K + 2 : 2 * K + 2, 1 ) ;                  % st. deviation of tastes  
  
  delta       = x * theta1 + vShock ;
  expmeanval0 = exp( delta ) ;
  expmu       = exp( x * diag( theta2 ) * v ) ;  
  
  %% Perform the mapping -- ideally this is the contraction mapping
  expmeanval1 = expmeanval0 .* share ./ ind_shnorm( expmeanval0, expmu ) ;
  
  %% Package the result
  delta1      = log( expmeanval1 ) ;
  vShock1     = delta1 - x * theta1 ;
  
  