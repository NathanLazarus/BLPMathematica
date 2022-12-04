% ComputeGMMEstimates - Computes GMM estimates for a BLP MC dataset
%
% This function computes the GMM estimates for a Monte Carlo BLP dataset
% from Dube, Fox, and Su using MPEC.  First generate the data with
% BLPSetup.m.
%
% modification history
% --------------------
% 29aug2010 bss starting guess based on best point estimates.
% 24aug2010 bss changed so SGI has same degree of exactness as 
%               monomial rule.
% 19aug2010 bss set qMC burn-in to 0.
% 05jul2010 bss kludge to fix namespace clash over iv vs. iv.m in new
%               version of Matlab.
% 25jan2010 bss written.
%

function [ ] = ComputeGMMEstimates( szRootDataDir, szDataDir, nQuadType, nMCDraws ) 

  iv = [] ;     % kludge to workaround file iv.m which clashes with variable iv in new version of Matlab

  szDate      = date ;
  szCurDate   = [ szDate(8:end) szDate(4:6) szDate(1:2) ] ;
  szDiaryName = sprintf( '%s/Solver-%s-QuadType%1dN%05d.log', szDataDir, szCurDate, nQuadType, nMCDraws ) ;

  % 0 -> pseudo-MC, 1 -> Gauss-Hermite, 2 -> Monomial, 3 -> quasi-MC, 4 -> Sparse Grids
  vQuadNames = [ ] ;
  vQuadNames{1} = 'pMC' ;
  vQuadNames{2} = 'Gauss-Hermite' ;
  vQuadNames{3} = 'Monomial' ;
  vQuadNames{4} = 'qMC' ;
  vQuadNames{5} = 'Sparse Grids' ;
  
  diary( szDiaryName ) ;
  
  disp( 'Computing GMM Estimates for BLP with Monte Carlo Data' ) ;
  
  fprintf( 1, '\n============================================================\n\n' ) ;
  fprintf( 1, '\tData Directory : %s\n', szDataDir ) ;
  fprintf( 1, '\tQuadrature Type: %s (%d)\n', vQuadNames{nQuadType+1}, nQuadType ) ;
  fprintf( 1, '\tNumber MC Draws: %d\n', nMCDraws ) ;
  fprintf( 1, '\n============================================================\n\n' ) ;


%% Globals Variables -- yuck :-(

  global share W PX x IV rc expmeanval delta
  global datasort marketForProducts sharesum oo
  global T prods nn K v rc rctrue oo oo1 tol_inner tol_outer
  global prodsMarket marketStarts marketEnds numProdsTotal

  global Q_NODES Q_WEIGHTS
  
  % Determine which quadrature type to use
  bUseGHQuad      = 0 ;
  bUseQuasiMC     = 0 ;
  bUseMonomial    = 0 ;
  bUseSparseGrids = 0 ;
  
  switch( nQuadType )
    case 0
      % Use pseudo-Monte Carlo
      true ;
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
      error('Unsupported quadrature type');
  end
  
  
%% Load Data and setup parameters

  szBLPConfigData = sprintf( '%s/BLPConfig.mat', szRootDataDir ) ;
  load( szBLPConfigData ) ;   % Load general configuration information
  load( strcat( szDataDir, '/MonteCarloBLPData.mat' ) ) ;         % Load data for a specific simulated dataset
  
  
  % Set seeds
  randn( 'seed', BLPConfig.vNormSeeds( 3, MySeedId ) ) ;
  rand( 'seed', BLPConfig.vUnifSeeds( 3, MySeedId ) ) ;
  
  
  nn = nMCDraws ;                                           % # draws to simulate shares

  tol_inner = BLPConfig.tol_inner ;
  tol_outer = BLPConfig.tol_outer ;
% %   tol_outer = 1e-4 ;    % XXX See if this improves convergence
  
  prods = BLPConfig.prods ;                                         % # products
  T     = BLPConfig.T ;                                         % # Markets

  starts = BLPConfig.nStarts ;                                         % # random start values to check during estimation

  Q_WEIGHTS = ones( nn, 1 ) / nn ;


%% Create Quadrature Nodes and Weights

% XXX Rework this as a switch statement sometime
if bUseGHQuad
  % Setup Gaussian Quadrature nodes and weights
  disp( '***Using GH Quadrature' ) ;
  
  GHQuadInit( length(betatrue), 7 ) ;  % Change the last argument to change the number of nodes

  v  = Q_NODES ;
  nn = size( v, 2 ) ;
  oo = ones(1,nn);                    % Redefine index for use in simulation of shares
  
elseif bUseMonomial
  disp( '***Using Monomial Rule 11-1' ) ;
  
% %   MonoQuadInit( 1 ) ;   % Choose left (1) or right (2) column

  % The right rule is more reliable than the left
  MonoQuadInit( 2 ) ;   % Choose left (1) or right (2) column

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


%% Back to Dube's code ...


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2SLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xhat = iv*inv(iv'*iv)*iv'*x;
PX2sls = inv(xhat'*xhat)*xhat';                     % project. matrix on weighted x-space for 2SLS
beta2sls = PX2sls*y;                                % y is y = log(share) - log(repmat(nopurch,prods,1));      % log-odds ratio
se2sls = sqrt(diag( mean((y-x*beta2sls).^2)*inv(xhat'*xhat) ));

expmeanval0 = exp(y);
% IV = [ones(T*prods,1) A z A.^2 A.^3 z.^2 z.^3 prod(A,2) prod(z,2) kron(A(:,1),ones(1,size(z,2))).*z  kron(A(:,2),ones(1,size(z,2))).*z];
IV = [ones(numProdsTotal,1) A z A.^2 A.^3 z.^2 z.^3 prod(A,2) prod(z,2) kron(A(:,1),ones(1,size(z,2))).*z  kron(A(:,2),ones(1,size(z,2))).*z];

W = inv(IV'*IV);
PX = inv(x'*IV*W*IV'*x)*x'*IV*W*IV';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPEC ESTIMATION OF RANDOM COEFFICIENTS LOGIT MODEL %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


bStartFromPointEst = 0 ;  % True to use starting values based on mean point estimate

if bStartFromPointEst
  nStarts = 10 ;
  vStart = CalcStartGuess( szRootDataDir ) ;
  
  theta20     = vStart( ( K + 2 ) : end ) ;
  startvalues = repmat( vStart( (K + 2) : (2 * K + 2) )', nStarts, 1 ) .* cat( 1, ones( 1, K + 1 ),      ...
                abs( ones( nStarts - 1, K + 1 ) + 0.5 * randn( nStarts - 1, K + 1 ) ) ) ;
else
  nStarts     = starts ;
  theta20     = 0.5*abs(beta2sls);
  startvalues = repmat(theta20',starts,1) .* cat( 2, ones( size(theta20) ), rand([ starts-1, size(theta20',2)] )' *1 )';
end


nIV = size(IV,2);       % # instrumental variables

GMPEC = 1.0e20;
CPUtMPEC = 0;
FuncEvalMPEC = 0;
GradEvalMPEC = 0;
expmeanval = expmeanval0;

mTheta1 = zeros( size( PX, 1 ), starts ) ;
mTheta2 = zeros( ( K + 1 ), starts ) ;
nFree   = sum( size( IV ) ) + size( PX, 1 ) + ( K + 1 ) ;
mXOpt   = zeros( nFree, starts ) ;    % Optimum x found by solver for each start
mFOpt = zeros( starts, 1 ) ;          % Value of GMM objective function at optimum for each start
mExitCode = zeros(starts, 1) ;        % Solver EXIT/INFO code

%% Compute point estimates for several starting values

nSaved = 0 ;  % number of optima which have been saved
for reps = 1 : nStarts
% % for reps=1:starts
  
    %%% NOTE: START VALUES SET TO GIVE SAME INITIAL GMM OBJECTIVE FUNCTIONAL
    %%% VALUE AS WITH NFP
    

    theta20 = startvalues( reps, : )' ;
    delta0  = invertshares( theta20 ) ;
    theta10 = PX*delta0 ;
    
    resid0 = delta0 - x * theta10 ;
    g0 = IV' * resid0 ;
    x0 = [theta10; theta20; delta0; g0];
    [ theta1MPEC_rep, theta2MPEC_rep, GMPEC_rep, INFOMPEC_rep, CPUtMPEC_rep, FuncEvalMPEC_rep, GradEvalMPEC_rep, xOpt ] =   ...
          runGMMMPECKnitro( x0, 1 ) ;
    
    CPUtMPEC = CPUtMPEC + CPUtMPEC_rep;
    FuncEvalMPEC = FuncEvalMPEC + FuncEvalMPEC_rep;
    GradEvalMPEC = GradEvalMPEC + GradEvalMPEC_rep;
    
    
    if ~bStartFromPointEst
      nSaved               = nSaved + 1 ;     % Number of optima which converged and have been saved
      mTheta1( :, nSaved ) = theta1MPEC_rep ;
      mTheta2( :, nSaved ) = theta2MPEC_rep ;
      mXOpt( :, nSaved )   = xOpt ;
      mFOpt( nSaved ) = GMPEC_rep ;
      mExitCode( nSaved) = INFOMPEC_rep ;
    end

    
    %% Store best solution found via MPEC
    if INFOMPEC_rep == 0 || INFOMPEC_rep == 1    % Valid solution

      if GMPEC_rep < GMPEC        % Best solution found so far (this compares f_k)        
        thetaMPEC1 = theta1MPEC_rep;
        thetaMPEC2 = theta2MPEC_rep;
        GMPEC      = GMPEC_rep;
        INFOMPEC   = INFOMPEC_rep;                
      end
      
      if bStartFromPointEst
        fprintf( 1, '===>>> Keeping iteration %d\n', reps ) ;
        nSaved               = nSaved + 1 ;     % Number of optima which converged and have been saved
        mTheta1( :, nSaved ) = theta1MPEC_rep ;
        mTheta2( :, nSaved ) = theta2MPEC_rep ;
        mXOpt( :, nSaved )   = xOpt ;
        mFOpt( nSaved ) = GMPEC_rep ;
        mExitCode( nSaved) = INFOMPEC_rep ;
      end
      
    end
    
    if nSaved == starts
      break ;
    end
    
end

if nSaved ~= starts
  fprintf( 2, '*** nSaved != starts: %d %d\n', nSaved, starts ) ;
end

% COVARIANCE MATRIX FOR MPEC STRUCTURAL PARAMETERS
if exist( 'thetaMPEC1', 'var' )
  delta = invertshares(thetaMPEC2);    % mean utilities
  resid = delta - x*thetaMPEC1;  % xi
  Ddelta = jacobian(thetaMPEC2);       % jacobian matrix of mean utilities : D delta_jt / D theta_2
  covg = zeros(size(IV,2));
  for ii =1:length(IV),
      covg = covg + IV(ii,:)'*IV(ii,:)*(resid(ii)^2);
  end
  Dg = [x Ddelta]'*IV;            % gradients of moment conditions
  covMPEC = inv( Dg*W*Dg')*Dg*W*covg*W*Dg'*inv( Dg*W*Dg');

  if true
    disp("\tCondition number and eigenvalues for covMPEC")
    cond( covMPEC )
    eig( covMPEC )
  end
  

  % Save results
  szOutFile = sprintf( '%s/PointEst-QuadType%1dN%05d.mat', szDataDir, nQuadType, nn ) ;
  save( szOutFile, 'share', 'W', 'PX', 'x', 'IV', 'rc', 'expmeanval',                         ...
        'datasort', 'marketForProducts', 'sharesum', 'oo', 'T', 'prods', 'nn', 'K',           ...
        'prodsMarket', 'marketStarts', 'marketEnds', 'numProdsTotal',                         ...
        'v', 'rctrue', 'oo1', 'tol_inner', 'tol_outer', 'mXOpt', 'mFOpt', 'mExitCode',        ...
        'covMPEC', 'betatrue', 'thetatrue', 'covrc', 'thetaMPEC1', 'thetaMPEC2' ) ;
    
else
  
  fprintf( 2, '\n\n********************\nDANGER: None of random starts solved!\n\n********************\n' ) ;
  
  szOutFile = sprintf( '%s/PointEst-QuadType%1dN%05d.mat', szDataDir, nQuadType, nn ) ;
  save( szOutFile, 'share', 'W', 'PX', 'x', 'IV', 'rc', 'expmeanval',                         ...
        'datasort', 'marketForProducts', 'sharesum', 'oo', 'T', 'prods', 'nn', 'K',           ...
        'prodsMarket', 'marketStarts', 'marketEnds', 'numProdsTotal',                         ...
        'v', 'rctrue', 'oo1', 'tol_inner', 'tol_outer', 'mXOpt', 'mFOpt', 'mExitCode',        ...
        'betatrue', 'thetatrue', 'covrc' ) ;
end

disp( 'Finished computing GMM Estimates for BLP with Monte Carlo Data' ) ;

diary off;


