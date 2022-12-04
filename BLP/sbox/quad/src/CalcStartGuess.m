% CalcStartGuess - computes starting guess for optimization
%
% Loads point estimates for all GH product rules and then takes the average
% of those which converged.  ComputeGMMEstimates.m uses this as the mean
% for computing the start values.
%
% modification history
% --------------------
% 29aug2010 bss written.
%

function [ vStart ] = CalcStartGuess( szRootDataDir )


%% Setup

szPointEstFile = 'PointEst-QuadType1N03125.mat' ;

% Valid optima for each dataset
mValidOpt{1} = [ 1,    3, 4, 5 ] ;
mValidOpt{2} = 1 : 5 ;
mValidOpt{3} = 1 : 5 ;
mValidOpt{4} = 1 : 5 ;
mValidOpt{5} = [ 1, 2,    4, 5 ] ;

%% Load each optimum by dataset and start value

szResultFile = fullfile( szRootDataDir, 'Seed0001', szPointEstFile ) ;
mSeed1       = load( szResultFile, 'mXOpt' ) ;

szResultFile = fullfile( szRootDataDir, 'Seed0002', szPointEstFile ) ;
mSeed2       = load( szResultFile, 'mXOpt' ) ;

szResultFile = fullfile( szRootDataDir, 'Seed0003', szPointEstFile ) ;
mSeed3       = load( szResultFile, 'mXOpt' ) ;

szResultFile = fullfile( szRootDataDir, 'Seed0004', szPointEstFile ) ;
mSeed4       = load( szResultFile, 'mXOpt' ) ;

szResultFile = fullfile( szRootDataDir, 'Seed0005', szPointEstFile ) ;
mSeed5       = load( szResultFile, 'mXOpt' ) ;


%% compute mean start value

nParam = 5 ;
vIx    = 1 : (2 * nParam ) ;
vIxVar = (nParam + 1 ) : (2 * nParam ) ;

mOpt1 = mSeed1.mXOpt( vIx, mValidOpt{ 1 } ) ;
mOpt2 = mSeed2.mXOpt( vIx, mValidOpt{ 2 } ) ;
mOpt3 = mSeed3.mXOpt( vIx, mValidOpt{ 3 } ) ;
mOpt4 = mSeed4.mXOpt( vIx, mValidOpt{ 4 } ) ;
mOpt5 = mSeed5.mXOpt( vIx, mValidOpt{ 5 } ) ;

% Fix sign on variance terms
mOpt1( vIxVar, : ) = abs( mOpt1( vIxVar, : ) ) ;
mOpt2( vIxVar, : ) = abs( mOpt2( vIxVar, : ) ) ;
mOpt3( vIxVar, : ) = abs( mOpt3( vIxVar, : ) ) ;
mOpt4( vIxVar, : ) = abs( mOpt4( vIxVar, : ) ) ;
mOpt5( vIxVar, : ) = abs( mOpt5( vIxVar, : ) ) ;

% Compute averate
vOpt1 = mean( mOpt1, 2 ) ;
vOpt2 = mean( mOpt2, 2 ) ;
vOpt3 = mean( mOpt3, 2 ) ;
vOpt4 = mean( mOpt4, 2 ) ;
vOpt5 = mean( mOpt5, 2 ) ;

vStart = mean( [ vOpt1, vOpt2, vOpt3, vOpt4, vOpt5 ], 2 ) ;

