% NiedSeq - Calculate the Niederreiter sequence of quasi-Monte Carlo points
%
% Args:
%   nPoints - how many points to generate
%   nDim    - how many dimensions
% mNS: an ( nPoints x nDim ) matrix of points
%
% modification history
% --------------------
% 30jun2009 bss written.
%

function [ mNS ] = NiedSeq( nPoints, nDim )

  % Generate numbers
  vNumList = ( 1 : nPoints )' ;
  
  % Generate powers of 2
  vTwos = 2 .^ ( ( 1 : nDim ) ./ ( nDim + 1 ) ) ;
  
  % Generate Niederreiter sequence.
  mNS = rem( vNumList * vTwos, 1 ) ;