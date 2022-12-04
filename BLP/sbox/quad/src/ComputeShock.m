% ComputeShock - computes product-market shock xi_jt
%
% modification history
% --------------------
% 14jan2010 bss fixed missing initialization of expmeanval
% 30nov2009 bss written.
%

function [ vShock, nIter ] = ComputeShock( x0 )

  global x K prods T expmeanval ;


  theta1 = x0( 1 : K + 1, 1 ) ;
  theta2 = x0( K + 2 : 2 * K + 2, 1 ) ;                  % st. deviation of tastes
  deltaSolver  = x0( 2 * K + 3 : 2 * K + 2 + T * prods, 1 ) ;         % mean utilities
  
  expmeanval = exp( deltaSolver ) ;
  
  [ delta, nIter ]  = invertshares( theta2 );    % mean utilities
  
  vShock = delta - x*theta1;  % xi
  
