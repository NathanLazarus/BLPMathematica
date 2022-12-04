% ComputeStdErr - Computes standard errors for point estimates
%
% Uses GMM formula to compute standard errors for point estimates
%
% modification history
% --------------------
% 24aug2010 bss changed so SGI has same degree of exactness as 
%               monomial rule.
% 23aug2010 bss written.
%

function [ ] = ComputeStdErr( szRootDataDir, szDataDir, nQuadType, nMCDraws )

%% Startup

  iv = [] ;     % kludge to workaround file iv.m which clashes with variable iv in new version of Matlab

  szDate      = date ;
  szCurDate   = [ szDate(8:end) szDate(4:6) szDate(1:2) ] ;
  szDiaryName = sprintf( '%s/StdErr-%s-QuadType%1dN%05d.log', szDataDir, szCurDate, nQuadType, nMCDraws ) ;

  diary( szDiaryName ) ;
  
  disp( 'Computing GMM Standard Errors for BLP with Monte Carlo Data' ) ;
  
  % 0 -> pseudo-MC, 1 -> Gauss-Hermite, 2 -> Monomial, 3 -> quasi-MC, 4 -> Sparse Grids
  vQuadNames = [ ] ;
  vQuadNames{1} = 'pMC' ;
  vQuadNames{2} = 'Gauss-Hermite' ;
  vQuadNames{3} = 'Monomial' ;
  vQuadNames{4} = 'qMC' ;
  vQuadNames{5} = 'Sparse Grids' ;
  
  fprintf( 1, '\n============================================================\n\n' ) ;
  fprintf( 1, '\tData Directory : %s\n', szDataDir ) ;
  fprintf( 1, '\tQuadrature Type: %s (%d)\n', vQuadNames{nQuadType+1}, nQuadType ) ;
  fprintf( 1, '\tNumber MC Draws: %d\n', nMCDraws ) ;
  fprintf( 1, '\n============================================================\n\n' ) ;


%% Globals Variables -- yuck :-(

  global share W PX x IV rc expmeanval delta
  global datasort denomexpand sharesum oo
  global T prods nn K v rc rctrue oo oo1 tol_inner tol_outer

  global Q_WEIGHTS ;
  
  % Determine which quadrature type to use
  bUseGHQuad      = 0 ;
  bUseQuasiMC     = 0 ;
  bUseMonomial    = 0 ;
  bUseSparseGrids = 0 ;
  
  switch( nQuadType )
    case 1 
      bUseGHQuad = 1 ;
    case 2
      % Use monomial rule
      bUseMonomial = 1 ;
    case 3
      % Use qMC
      bUseQuasiMC = 1 ;
    case 4
      % Use Sparse Grids
      bUseSparseGrids = 1 ;
    otherwise
      % Use pseudo-Monte Carlo
  end
  
  
%% Load Data and setup parameters

  szBLPConfigData = fullfile( szRootDataDir, 'BLPConfig.mat' ) ;
  load( szBLPConfigData ) ;   % Load general configuration information
  load( fullfile( szDataDir, 'MonteCarloBLPData.mat' ) ) ;         % Load data for a specific simulated dataset
  
  
  % Set seeds
  randn( 'seed', BLPConfig.vNormSeeds( 3, MySeedId ) ) ;
  rand( 'seed', BLPConfig.vUnifSeeds( 3, MySeedId ) ) ;
  
  
  nn = nMCDraws ;                                           % # draws to simulate shares

  tol_inner = BLPConfig.tol_inner ;
  tol_outer = BLPConfig.tol_outer ;
% %   tol_outer = 1e-4 ;    % XXX See if this improves convergence
  
  prods = BLPConfig.prods ;                                         % # products
  T     = BLPConfig.T ;                                         % # Markets

  nStarts = BLPConfig.nStarts ;                                         % # random start values to check during estimation

  Q_WEIGHTS = ones( nn, 1 ) / nn ;


  %% Create Quadrature Nodes and Weights

  % XXX Rework this as a switch statement sometime
  if bUseGHQuad
    % Setup Gaussian Quadrature nodes and weights
    disp( '***Using GH Quadrature' ) ;

    GHQuadInit( betatrue, covrc, 7 ) ;  % Change the last argument to change the number of nodes

    global Q_NODES ;
    v  = Q_NODES ;
  %   rc = chol(covrc)' * v ;
    nn = size( v, 2 ) ;
    oo = ones(1,nn);                    % Redefine index for use in simulation of shares

  elseif bUseMonomial
    disp( '***Using Monomial Rule 11-1' ) ;

  % %   MonoQuadInit( 1 ) ;   % Choose left (1) or right (2) column

    % The right rule is more reliable than the left
    MonoQuadInit( 2 ) ;   % Choose left (1) or right (2) column

    global Q_NODES ;
    v  = Q_NODES ;    % Use Gaussian-Hermite nodes

    nn = size( v, 2 ) ;
    oo = ones( 1, nn ) ; 

  elseif bUseQuasiMC
    disp( '***Using qMC' ) ;

    nPoints = nMCDraws ;
    % nBurnIn = 10000 ;
    nBurnIn = 0 ;
    nDim    = size( v, 1 ) ;

    mNied = NiedNormQuadInit( nPoints, nDim, nBurnIn ) ;

    v  = mNied ;
    nn = size( v, 2 ) ;
    oo = ones( 1, nn ) ;
    Q_WEIGHTS = ones( nn, 1 ) / nn ;

  elseif bUseSparseGrids

    disp( '***Using Sparse Grids' ) ;

    addpath( '~/Tools/SparseGrids' ) ;

    % Konrad-Patterson Sparse Grid for Normal Kernel in 5 dim with K = 6
    % K = 6 => accuracy for total order 2*K-1
    % So accuracy is same as monomial rule
    [ Q_NODES, Q_WEIGHTS ] = nwspgr( 'KPN', 5, 6 ) ;      

    v  = Q_NODES' ;
    nn = size( v, 2 ) ;
    oo = ones( 1, nn ) ;

  else
    % configure number of simulation draws for solution of MPEC
    disp( '***Using MC' ) ;

  % %   global Q_WEIGHTS ;
    nn = nMCDraws ;       % BLPConfig.N_DRAWS_MC_INT ;
    oo = ones(1,nn);                    % Redefine index for use in simulation of shares
    Q_WEIGHTS = ones( nn, 1 ) / nn ;

    % resimulate the shock, because the econometrician does not observe it.
    v = randn(length(betatrue),nn);                     % draws for share integrals during estimation
  end


  %% Load Solver Results in order to compute standard errors

  szSolverFile = sprintf( '%s/PointEst-QuadType%1dN%05d.mat', szDataDir, nQuadType, nn ) ;
  load( szSolverFile ) ;

  %% Compute Standard Errors

  if nStarts ~= size( mXOpt, 2 )
    error( 'nStarts nonconformable' ) ;
  end

  nParams    = K + 1 ; % number of dimensions.  Total number of parameters == 2 * nParams
  nTotParams = 2 * nParams ;
  mVar    = zeros( nTotParams, nTotParams, nStarts ) ;
  mBLPVar = zeros( nTotParams, nTotParams, nStarts ) ;

  
  mStdErr    = zeros( nTotParams, nStarts ) ;
  mBLPStdErr = zeros( nTotParams, nStarts ) ;
  
  for ixStart = 1 : nStarts
    
    % Compute standard errors two ways:  using mean utilities from optimal
    % solution and from Berry's mapping
    
    fprintf( 1, '\n==============================\n' ) ;
    fprintf( 1, '\tComputing Standard errors for Start = %d\n', ixStart ) ;
    
    vTheta1 = mXOpt( 1 : nParams, ixStart ) ;
    vTheta2 = mXOpt( ( nParams + 1 ) : nTotParams, ixStart ) ;
    vDelta  = mXOpt( (nTotParams + 1) : (nTotParams + T * prods ), ixStart ) ;
    
    
    %------------------------
    % 1: Use solver mean utilities    
    
    delta = vDelta ;   % Use mean utilities
    resid = delta - x * vTheta1 ;  % xi
    Ddelta = jacobian( vTheta2 ) ;       % jacobian matrix of mean utilities
    covg = zeros( size( IV, 2) ) ;
    for ii = 1 : length( IV )
        covg = covg + IV(ii,:)' * IV(ii,:) * ( resid(ii)^2 ) ;
    end
    Dg = [x Ddelta]' * IV ;            % gradients of moment conditions
    covMPEC = inv( Dg*W*Dg' ) * Dg * W * covg * W * Dg' * inv( Dg*W*Dg' ) ;

    mVar( :, :, ixStart ) = covMPEC ;
    mStdErr( :, ixStart ) = sqrt( diag( covMPEC ) ) ;
    
    fprintf( 1, '\tStandard Errors :\n' ) ;
    mStdErr( :, ixStart )'
    if any( isnan( covMPEC(:) ) )
      fprintf( 2, '\tVariance matrix is NaN\n' ) ;
    else
      fprintf( 1, '\tCondition Number :\t\t%g\n', cond( covMPEC ) ) ;
      fprintf( 1, '\tEigen Values :\n' ) ;
      eig( covMPEC )' 
    end
    
    
    fprintf( 1, '\n------------------------------\n\n' ) ;
    
    %------------------------
    % 2: Use Berry map to compute mean utilities
 bComputeBLPVar = 0 ;    
 if bComputeBLPVar
    delta = invertshares(thetaMPEC2);    % mean utilities
    
    resid = delta - x * vTheta1 ;  % xi
    Ddelta = jacobian( vTheta2 ) ;       % jacobian matrix of mean utilities
    covg = zeros( size( IV, 2) ) ;
    for ii = 1 : length( IV )
        covg = covg + IV(ii,:)' * IV(ii,:) * ( resid(ii)^2 ) ;
    end
    Dg = [x Ddelta]' * IV ;            % gradients of moment conditions
    covMPEC = inv( Dg*W*Dg' ) * Dg * W * covg * W * Dg' * inv( Dg*W*Dg' ) ;

    mBLPVar( :, :, ixStart ) = covMPEC ;
    mBLPStdErr( :, ixStart ) = sqrt( diag( covMPEC ) ) ;
    
    fprintf( 1, '\tBLP Standard Errors :\n' ) ;
    mBLPStdErr( :, ixStart )'
    fprintf( 1, '\tBLP Condition Number :\t\t%g\n', cond( covMPEC ) ) ;
    fprintf( 1, '\tBLP Eigen Values :\n' ) ;
    eig( covMPEC )' 
    
    fprintf( 1, '\n==============================\n\n\n' ) ;
 else
   disp( 'Skipping BLP Variance Computation' ) ;
   mBLPVar    = [ ] ;   
   mBLPStdErr = [ ] ;
 end
 
  end

  %% Save Results

  szStdErrFile = sprintf( '%s/Variance-QuadType%1dN%05d.mat', szDataDir, nQuadType, nn ) ;
  save( szStdErrFile, 'mVar', 'mBLPVar', 'mStdErr', 'mBLPStdErr' ) ;
  
  fprintf( 1, 'MANEUVER COMPLETE\n' ) ;
  
  
