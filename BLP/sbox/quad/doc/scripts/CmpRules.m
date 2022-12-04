% CmpRules - compares integration rules performance for monmials.
%
% Computes various monomial moments using different rules to
% compare their performance.
%
% modification history
% --------------------
% 27jan2011 bss written.
%

%% Setup

% % szDataBase = '/Users/bss/sbox/blpcpp/data/quad' ;
% % szSGI   = fullfile( szDataBase, 'SGIRule.5-6.csv' ) ;
% % szMonoL = fullfile( szDataBase, 'Mono.11-1.L.csv' ) ;
% % szMonoR = fullfile( szDataBase, 'Mono.11-1.R.csv' ) ;
% % szGH3   = fullfile( szDataBase, 'GHProdRule.3-5.csv' ) ;
% % szGH5   = fullfile( szDataBase, 'GHProdRule.5-5.csv' ) ;
% % szGH7   = fullfile( szDataBase, 'GHProdRule.7-5.csv' ) ;
% % szPMC   = fullfile( szDataBase, 'pMCRule.10000.csv' ) ;

Config.szOutDir  = '/Users/bss/sbox/quad/doc/tables' ;
Config.szOutFile = fullfile( Config.szOutDir, 'CmpRules.mat' ) ;

Config.vGHNodes  = [ 3, 5, 7 ] ; % Number of nodes in each dimension for GH product rule
Config.nRandCoef = 5 ;      % number of random coefficients
Config.nSGIAccuracy = 6 ;   % K parameter for setting sparse grid accuracy 2K-1
Config.nSeed     = 9310 ;   % Seed for randn
Config.nMCDraws  = 10000 ;  % Number of pMC draws
Config.nMCRep    = 100 ;    % Number of pMC replications


% Set seed for pMC draws
randn( 'seed', Config.nSeed ) ;


% specify path to relevant files
addpath( '/Users/bss/sbox/quad/src' ) ;
addpath( '/Users/bss/Tools/SparseGrids' ) ;


%% Generate nodes and weights for each rule

Rules = [ ] ;     % Structure to manage all rules

global Q_WEIGHTS
global Q_NODES

%---------------------
% Gauss-Hermite

GHQuadInit( ones( Config.nRandCoef, 1 ), 0, 3 ) ;
Rules.GH3.vWeights = Q_WEIGHTS ;
Rules.GH3.mNodes   = Q_NODES' ;

GHQuadInit( ones( Config.nRandCoef, 1 ), 0, 5 ) ;
Rules.GH5.vWeights = Q_WEIGHTS ;
Rules.GH5.mNodes   = Q_NODES' ;

GHQuadInit( ones( Config.nRandCoef, 1 ), 0, 7 ) ;
Rules.GH7.vWeights = Q_WEIGHTS ;
Rules.GH7.mNodes   = Q_NODES' ;


%---------------------
% Sparse Grids

[ Rules.SGI.mNodes, Rules.SGI.vWeights ] = nwspgr( 'KPN', Config.nRandCoef, Config.nSGIAccuracy ) ;      

%---------------------
% Monomial Rule 11-1 Left

MonoQuadInit( 1 ) ;   % Choose left (1) or right (2) column
Rules.MonoL.vWeights = Q_WEIGHTS ;
Rules.MonoL.mNodes   = Q_NODES' ;

%---------------------
% Monomial Rule 11-1 Right

MonoQuadInit( 2 ) ;   % Choose left (1) or right (2) column
Rules.MonoR.vWeights = Q_WEIGHTS ;
Rules.MonoR.mNodes   = Q_NODES' ;

%---------------------
% pMC Rule

nn          = Config.nMCDraws ;  
vPMCWeights = ones( nn, 1 ) / nn ;

Rules.pMC = cell( Config.nMCRep, 1 ) ;
for ix = 1 : Config.nMCRep 
  Rules.pMC{ ix }.vWeights = vPMCWeights ;
  Rules.pMC{ ix }.mNodes   = randn( nn, Config.nRandCoef ) ; 
end


%---------------------
% Halton Rule

nn          = Config.nMCDraws ;  
vPMCWeights = ones( nn, 1 ) / nn ;

qrsHalton = qrandstream( 'halton', Config.nRandCoef, 'Skip', 1e3, 'Leap', 1e2 ) ;

Rules.qMC = cell( Config.nMCRep, 1 ) ;
for ix = 1 : Config.nMCRep 
  Rules.qMC{ ix }.vWeights = vPMCWeights ;
  Rules.qMC{ ix }.mNodes   = norminv( qrand( qrsHalton, nn ), 0, 1 ); 
end



%% Cases to compute

Config.nCases = 20 ;

Cases  = zeros( Config.nCases, Config.nRandCoef ) ;     % exponents of terms in monomial
vTruth = zeros( Config.nCases, 1 ) ;                    % Correct answer of integral of monomial

% Correct values of 1-d moments for n = 0 : 16
vMoments = [ 1, 0, 1, 0, 3, 0, 15, 0, 105, 0, 945, 0, 10395, 0, 135135, 0, 2027025, 0 ] ;

ixCase = 1 ;
Cases( ixCase, : )  = [ 0, 0 , 0, 0, 0 ] ;
vTruth( ixCase )    = 1 ;

ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 1, 0, 0, 0, 0 ] ;
vTruth( ixCase )    = 0 ;

ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 1, 1, 0, 0, 0 ] ;
vTruth( ixCase )    = 0 ;

ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 1, 1, 1, 1, 1 ] ;
vTruth( ixCase )    = 0 ;

ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 2, 0, 0, 0, 0 ] ;
vTruth( ixCase )    = 1 ;

ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 4, 0, 0, 0, 0 ] ;
vTruth( ixCase )    = 3 ;

ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 6, 0, 0, 0, 0 ] ;
vTruth( ixCase )    = 15 ;

ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 0, 6, 0, 4, 0 ] ;
vTruth( ixCase )    = 45 ;

ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 10, 0, 0, 0, 0 ] ;
vTruth( ixCase )    = prod( vMoments( Cases( ixCase, : ) + 1 ) ) ;

% degree p = 11 Monomial
ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 5, 4, 2, 0, 0 ] ;
vTruth( ixCase )    = prod( vMoments( Cases( ixCase, : ) + 1 ) ) ;

% degree p = 12 Monomial
ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 12, 0, 0, 0, 0 ] ;
vTruth( ixCase )    = prod( vMoments( Cases( ixCase, : ) + 1 ) ) ;

% p = 13
ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 13, 0, 0, 0, 0 ] ;
vTruth( ixCase )    = prod( vMoments( Cases( ixCase, : ) + 1 ) ) ;

% p = 14
ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 14, 0, 0, 0, 0 ] ;
vTruth( ixCase )    = prod( vMoments( Cases( ixCase, : ) + 1 ) ) ;

% p = 15
ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 15, 0, 0, 0, 0 ] ;
vTruth( ixCase )    = prod( vMoments( Cases( ixCase, : ) + 1 ) ) ;

% p = 16
ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 16, 0, 0, 0, 0 ] ;
vTruth( ixCase )    = prod( vMoments( Cases( ixCase, : ) + 1 ) ) ;

% p = 20
ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 6, 6, 4, 2, 2 ] ;
vTruth( ixCase )    = prod( vMoments( Cases( ixCase, : ) + 1 ) ) ;

% p = 22
ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 8, 6, 4, 2, 2 ] ;
vTruth( ixCase )    = prod( vMoments( Cases( ixCase, : ) + 1 ) ) ;

% p = 23
ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 10, 5, 4, 2, 2 ] ;
vTruth( ixCase )    = prod( vMoments( Cases( ixCase, : ) + 1 ) ) ;

% p = 32
ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 10, 10, 6, 4, 2 ] ;
vTruth( ixCase )    = prod( vMoments( Cases( ixCase, : ) + 1 ) ) ;

% p = 38
ixCase = ixCase + 1 ;
Cases( ixCase, : )  = [ 16, 12, 4, 4, 2 ] ;
vTruth( ixCase )    = prod( vMoments( Cases( ixCase, : ) + 1 ) ) ;


%% Compute for each rule

nPolyRules = 1 + 2 + length( Config.vGHNodes ) ;
nTotRules  = Config.nMCRep + nPolyRules ;
global mErr
mErr       = zeros( Config.nCases, nTotRules ) ;

mErr( :, 1 ) = ComputeMonomial( Cases, Rules.GH3 ) - vTruth ; 
mErr( :, 2 ) = ComputeMonomial( Cases, Rules.GH5 ) - vTruth ;
mErr( :, 3 ) = ComputeMonomial( Cases, Rules.GH7 ) - vTruth ;

mErr( :, 4 ) = ComputeMonomial( Cases, Rules.SGI ) - vTruth ;

mErr( :, 5 ) = ComputeMonomial( Cases, Rules.MonoL ) - vTruth ;
mErr( :, 6 ) = ComputeMonomial( Cases, Rules.MonoR ) - vTruth ;

for ix = 1 : Config.nMCRep
  mErr( :, ix + nPolyRules ) = ComputeMonomial( Cases, Rules.pMC{ix} ) - vTruth ;
end


%% Compute Statistics on pMC Errors

mErrPMC = mErr( :, 1+nPolyRules : end ) ;

vMean   = mean( abs( mErrPMC ), 2 ) ;
% % vMedian = median( abs( mErrPMC ), 2 ) ;
% % vQ25    = quantile( abs( mErrPMC ), .25, 2 ) ;
% % vQ75    = quantile( abs( mErrPMC ), .75, 2 ) ;

vStd    = std( mErrPMC, 0, 2 ) ;


%% qMC

mQErr = zeros( Config.nCases, Config.nMCRep ) ;
for ix = 1 : Config.nMCRep
  mQErr( :, ix ) = ComputeMonomial( Cases, Rules.qMC{ix} ) - vTruth ;
end

vQMean = mean( abs( mQErr ), 2 ) ;
vQStd  = std( mQErr, 0, 2 ) ;

mQCmpErr = [ mErr(:,1:nPolyRules) vQMean vQStd ] ;


%% Save

% mCmpErr = [ mErr(:,1:nPolyRules) vMean vQ25 vMedian vQ75 ] ;

mCmpErr = [ mErr(:,1:nPolyRules) vMean vStd ] ;

save( Config.szOutFile, 'mCmpErr', 'mQCmpErr', 'Cases' ) ;