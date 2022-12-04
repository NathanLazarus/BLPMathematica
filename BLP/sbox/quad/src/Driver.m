% Driver.m - Run several cases of model
%
% Runs BLP quadrature calculations.
%
% modification history
% --------------------
% 01feb2012 bss relocated data to data dir.
% 05may2010 bss improved configurability.
% 16jan2010 bss written.
%

%% Access latest KNITRO
addpath('/opt/knitro/12.3.0/knitromatlab')  % BSS XXX

%% Initial Configuration

nSeed = 5329 ;          % Used to generate all other seeds and datasets
if false
    % Debug
    nSim = 2;
    nProducts = 25;
    nMarkets = 10;
else
    nSim  = 5 ;             % Number of simulated datasets

    nProducts = 25 ;                % Number of products J
    nMarkets = 50 ;                 % Number of markets T
end

szRootDataDir = sprintf( '../data/DataJ%dT%d', nProducts, nMarkets ) ;


%% Create Data Sets

BLPSetup( szRootDataDir, nProducts, nMarkets, nSeed, nSim ) ;   % Create datasets and ./Data subdirectory of output


%% Get Point Estimates
%
% Solves each model via GMM + MPEC to get initial points
%
% NOTE: You must generate point estimates for each Seed0000 directory
% before computing share integrals.
%

nQuadType = 1 ;       % 0 -> pseudo-MC, 1 -> Gauss-Hermite, 2 -> Monomial, 3 -> quasi-MC, 4 -> Sparse Grids
nMCDraws  = 1000 ;    % number of draws for pMC or qMC

for ixSim = 1 : nSim
  
%   szDataDir = sprintf( './Data/Seed%04d', ixSim ) ;
  szDataDir = sprintf( '%s/Seed%04d', szRootDataDir, ixSim ) ;
  
  ComputeGMMEstimates( szRootDataDir, szDataDir, nQuadType, nMCDraws ) ;       
% %   ComputeGMM1StartNDraws( szRootDataDir, szDataDir, nQuadType, nMCDraws ) ;   
end


%% Compute Standard Errors

for ixSim = 1 : nSim
  
%   szDataDir = sprintf( './Data/Seed%04d', ixSim ) ;
  szDataDir = sprintf( '%s/Seed%04d', szRootDataDir, ixSim ) ;
  
  ComputeStdErr( szRootDataDir, szDataDir, nQuadType, nMCDraws ) ;       
end


%% Compute Share Integrals
%
% Compute share integrals for a variety of scenarios
%

for ixSim = 1 : nSim
% % for ixSim = 1 : 5
  
%   szDataDir = sprintf( './Data/Seed%04d', ixSim ) ;
  szDataDir = sprintf( '%s/Seed%04d', szRootDataDir, ixSim ) ;
  
  nQuadType   = 1 ;       % Must agree with parameters to ComputeGMMEstimates
% %   nSolMCDraws = 1000 ;    % Must agree with parameters to ComputeGMMEstimates
  nSolMCDraws = 7^5 ;
  % nQMCBurnIn  = 10000 ;   % Number of points to discard when generating quasi-MC points
  nQMCBurnIn  = 0 ;   % Number of points to discard when generating quasi-MC points
  
  % Compute share integrals and compare them.  
  % Warning: this will produce graphs which are stored in a PostScript file
    
  nMCDraws    = 100 ;
  nQMCDraws   = 100 ;
% %   CompareShares( szRootDataDir, szDataDir, nQuadType, nSolMCDraws, nMCDraws, nQMCDraws, nQMCBurnIn ) ;
  
  nMCDraws    = 1000 ;
  nQMCDraws   = 1000 ;
% %   CompareShares( szRootDataDir, szDataDir, nQuadType, nSolMCDraws, nMCDraws, nQMCDraws, nQMCBurnIn ) ;
  
  nMCDraws    = 10000  ;   % Number of pseudo-MC points to use in integration
  nQMCDraws   = 10000  ;   % Number of quasi-MC points to use in integration
  CompareShares( szRootDataDir, szDataDir, nQuadType, nSolMCDraws, nMCDraws, nQMCDraws, nQMCBurnIn ) ;
  
end

