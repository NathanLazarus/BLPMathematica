% BLPSetup - Configures Monte Carlo Data Sets for Dube-Fox-Su BLP
%
% This function generates Monte Carlo datasets using the BLP model.  To facilite running
% multiple simulations, data is stored for different seeds under the `Data' subdirectory.
% Data for each seed is stored in ./Data/Seed<XX>.
%
% Configuration information is stored in a BLSConfig structure in ./Data/BLPConfig.mat.
%
% Shamelessly based on Dube-Fox-Su example code.
%
% modification history
% --------------------
% 05may2010 bss added parameters for more configurability.
% 25jan2010 bss written.
%


function [] = BLPSetup( szRootDataDir, nProducts, nMarkets, nOrigSeed, nSim )


%% SETTINGS


%------------------
% Configuration for Simulation

BLPConfig.nOrigSeed = nOrigSeed ;    % Original seed to generate all others  (5329)
BLPConfig.nSim      = nSim ;         % Number of simulated datasets
BLPConfig.nSeeds    = 3 ;            % Number of seeds needed per dataset (2 for data generation, 1 for pMC )

randn( 'seed', BLPConfig.nSim ) ;

BLPConfig.vNormSeeds = floor( rand( BLPConfig.nSeeds, BLPConfig.nSim ) * 1e4 );
BLPConfig.vUnifSeeds = floor( rand( BLPConfig.nSeeds, BLPConfig.nSim ) * 1e4 );

BLPConfig.tol_inner = 1.e-14 ;
BLPConfig.tol_outer = 1.e-6 ;
  
BLPConfig.N_DRAWS_GEN_DATA = 100 ;        % How many draws to generate the data
BLPConfig.N_DRAWS_MC_INT   = 5000 ;      % How many draws to simulate the integral

% BLPConfig.prods = 15 ;    % Number of producs (25)
% BLPConfig.T     = 20 ;    % Number of markets (50)
BLPConfig.prods = nProducts ;    % Number of producs (25)
BLPConfig.T     = nMarkets ;    % Number of markets (50)

BLPConfig.nStarts = 5; % XXX 1 ;   % Number of random starts for solver

% BLPConfig.szDataDir = './Data' ;
BLPConfig.szDataDir = szRootDataDir ;

% Prepare subdirectory
if 7 ~= exist( BLPConfig.szDataDir, 'dir' ) 
  disp( 'Creating BLP Data Dir' ) ;
  mkdir( BLPConfig.szDataDir ) ;
end


% Save basic configuration for all simulations
disp( 'Saving BLP Config Data' ) ;
save( sprintf( '%s/BLPConfig.mat', BLPConfig.szDataDir ), 'BLPConfig' ) ;


%------------------
% Setup variables

global share W PX x IV rc expmeanval
global datasort marketForProducts sharesum oo
global T prods nn K v rc rctrue oo tol_inner tol_outer
global prodsMarket marketStarts marketEnds numProdsTotal

szDate      = date ;
szDiaryName = [ 'Log' szDate(8:end) szDate(4:6) szDate(1:2) '.out' ] ;


%% Create Simulated Datasets

for ixSim = 1 : BLPConfig.nSim

  % randn('seed',155)                                   % reset normal random number generator
  % rand('seed',155)                                    % reset uniform random number generator
  % randn('seed',231)                                   % reset normal random number generator
  % rand('seed',159)                                    % reset uniform random number generator
  
  randn( 'seed', BLPConfig.vNormSeeds( 1, ixSim ) ) ;
  rand( 'seed', BLPConfig.vUnifSeeds( 1, ixSim ) ) ;

  nn = BLPConfig.N_DRAWS_GEN_DATA ;                                           % # draws to simulate shares

  tol_inner = BLPConfig.tol_inner ;
  tol_outer = BLPConfig.tol_outer ;

  prods = BLPConfig.prods ;                                         % # products
  T     = BLPConfig.T ;                                         % # Markets
  numProdsTotal = prods * T

  starts = BLPConfig.nStarts ;                                         % # random start values to check during estimation


  %% Number of products per market
  % Discrete uniform between 1 and prods

  % Set to true to use a variable number of products per market
  if false
      prodsMarket = randi([1 prods],T,1);  % Number of products per market,
  else
      prodsMarket = repelem(prods, T);
  end
  
  
  % Use two dimensional matrices for speed reasons (vectorization)
  % Market indices tell us the first and last spots of records for each
  % market in list of all products
  marketStarts = zeros(T,1);
  marketEnds = zeros(T,1);
  marketStarts(1) = 1;
  marketEnds(1) = prodsMarket(1);
  for t=2:T
      marketStarts(t) = marketEnds(t-1) + 1;
      marketEnds(t) = marketStarts(t) + prodsMarket(t) - 1;
  end
  numProdsTotal = marketEnds(T);  % Total number of products in all markets

  %% TRUE PARAMETER VALUES

  costparams = ones(6,1)/2;
  betatrue = [2 1.5 1.5 .5 -3]';                     % true mean tastes
  betatrue0 = betatrue;
  K = size(betatrue,1)-1;                             % # attributes = length( beta ) - 1 (for constant)
  covrc = diag( [ .5 .5 .5 .5 .2] );                  % true co-variances in tastes
  rctrue = covrc(find(covrc));
  thetatrue = [betatrue;sqrt(rctrue)];

  BLPConfig.vThetaTrue = thetatrue ;

  v = randn(length(betatrue),nn);                     % draws for share integrals during estimation
  rc = chol(covrc)'*v;                                % draws for share integrals for data creation
  sigmaxi=1;
  covX = -.8*ones(K-1)+1.8*eye(K-1);                  % allow for correlation in X's
  covX(1,3) = .3;
  covX(2,3) = .3;
  covX(3,1) = .3;
  covX(3,2) = .3;

  global Q_WEIGHTS ;
  Q_WEIGHTS = ones( nn, 1 ) / nn ;


  %% Create indices to speed share calculations
  oo = ones(1, nn);                                    % index for use in simulation of shares
  
  sharesum = sparse(zeros(T,numProdsTotal));  % used to create denominators in logit predicted shares (i.e. sums numerators)
  for t=1:T
      sharesum(t,marketStarts(t):marketEnds(t)) = 1;
  end

  marketForProducts = zeros(numProdsTotal,1);  % market indices for each product
  % used to expand one logit denominator per consumer to one for each product
  for t=1:T
      marketForProducts(marketStarts(t):marketEnds(t)) = t;
  end


  %% SIMULATE DATA

% %   randn('seed',3837)                                   % reset normal random number generator
% %   rand('seed',2234)                                    % reset uniform random number generator

  randn( 'seed', BLPConfig.vNormSeeds( 2, ixSim ) ) ;
  rand( 'seed', BLPConfig.vUnifSeeds( 2, ixSim ) ) ;


  
  xi = randn( numProdsTotal,1)*sigmaxi;                % draw demand shocks
  A = randn(numProdsTotal,K-1)*chol(covX);             % product attributes
  prand = rand(numProdsTotal,1)*5;
  price = 3 +   xi*1.5 +  prand + sum(A,2);
  z = rand(numProdsTotal,length(costparams)) + 1/4*repmat( abs( prand +  sum(A,2)*1.1 ) ,1,length(costparams));
  x = [ones(numProdsTotal,1) A price];
  [share,nopurch] = mksharesim(betatrue,x,xi,rc);
  y = log(share) - log(nopurch(marketForProducts,:) );     % log-odds ratio for no RC logit estimation
  % matrix of instruments
  % XXX BSS Orig XXX %iv = [ones(prods*T,1) A z];                         % matrix of instruments
  iv = [ones(numProdsTotal,1) A z A.^2 A.^3 z.^2 z.^3 prod(A,2) prod(z,2) kron(A(:,1),ones(1,size(z,2))).*z  kron(A(:,2),ones(1,size(z,2))).*z];
  IV = iv;


  %% Store data
  
  szMCDataDir  = sprintf( '%s/Seed%04d', BLPConfig.szDataDir, ixSim ) ;
  
  if 7 ~= exist( szMCDataDir, 'dir' ) 
    disp( strcat( 'Creating ', szMCDataDir ) ) ;
    mkdir( szMCDataDir ) ;
  end
  
  szMCDataFile = strcat( szMCDataDir, '/MonteCarloBLPData.mat' ) ;
  disp( strcat( 'Saving ', szMCDataFile ) ) ;
  
  MySeedId = ixSim ;      % Which seeds were used to create the data
  save( szMCDataFile, 'share', 'W', 'PX', 'x', 'IV', 'rc', 'expmeanval',                ...
        'datasort', 'marketForProducts', 'sharesum', 'oo', 'T', 'prods', 'nn', 'K',     ...
        'v', 'rctrue', 'tol_inner', 'tol_outer','betatrue', 'thetatrue',                ...
        'covrc', 'iv', 'y', 'A', 'z', 'MySeedId',                                       ...
        'prodsMarket', 'marketStarts', 'marketEnds', 'numProdsTotal');

end


disp( 'Finished creating datasets!' ) ;


