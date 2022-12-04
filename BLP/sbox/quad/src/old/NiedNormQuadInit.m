% NiedNormQuadInit - create a set of Niedereiter quasi-Monte Carlo points
%
% NOTE: these points are normally distributed ~ N( 0, eye( nDim ) )
%
% PARAMS
%
%   nPoints     - how many points
%   nDim        - number of dimensions
%   nBurnIn     - how many points to burn in
%
% modification history
% --------------------
% 19jan2010 bss written.
%


function [ mNied ] = NiedNormQuadInit( nPoints, nDim, nBurnIn )

% addpath( '../SimpleModel/Tools' ) ;   % Path to Niedereiter code

% Get raw points in [ -1, 1 ]
mRaw = NiedSeq( nPoints + nBurnIn, nDim ) ;

% Discard burnin points and convert to normal distribution
mNied = norminv( mRaw( nBurnIn + 1 : end, : )', 0, 1 ) ;

