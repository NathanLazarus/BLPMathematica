% ComputeMonomial - computes a monomial using an arbitrary quadrature rule
%
% modification history
% --------------------
% 27jan2011 bss written.
%

function [ vResult ] = ComputeMonomial( mCases, Rule )

  [ nCases, nRandCoef ] = size( mCases ) ;
  
  vResult = zeros( nCases, 1 ) ;
  
  nNodes  = length( Rule.vWeights ) ;
  vOnes   = ones( nNodes, 1 ) ;
  
  for ix = 1 : nCases
    mExponent     = vOnes * mCases( ix, : ) ;
    vResult( ix ) = Rule.vWeights' * prod( Rule.mNodes .^ mExponent, 2 ) ;
  end
  