% ComputeGMM1StartNDraws - Computes GMM estimates for different draws
%
% This function computes the GMM estimates for a Monte Carlo BLP dataset
% from Dube, Fox, and Su using MPEC.  This code investigates how, for a 
% given starting value, different pMC draws affects convergence of the solver
% First generate the data with BLPSetup.m.
%
% modification history
% --------------------
% 05jul2010 bss kludge to fix namespace clash over iv vs. iv.m in new
%               version of Matlab.
% 18jun2010 bss hacked from ComputeGMMEstimates.m.
%

function [ ] = ComputeGMM1StartNDraws( szRootDataDir, szDataDir, nQuadType, nMCDraws ) 
 
  iv = [] ;     % kludge to workaround file iv.m which clashes with variable iv in new version of Matlab
  
  szDate      = date ;
  szCurDate   = [ szDate(8:end) szDate(4:6) szDate(1:2) ] ;
  szDiaryName = sprintf( '%s/NDraws-%s-QuadType%1dN%05d.log', szDataDir, szCurDate, nQuadType, nMCDraws ) ;

  diary( szDiaryName ) ;
  
  disp( 'Computing GMM Estimates for BLP with Monte Carlo Data' ) ;
  
  fprintf( 1, '\n============================================================\n\n' ) ;
  fprintf( 1, '\tData Directory : %s\n', szDataDir ) ;
  fprintf( 1, '\tQuadrature Type: %d\n', nQuadType ) ;
  fprintf( 1, '\tNumber MC Draws: %d\n', nMCDraws ) ;
  fprintf( 1, '\n============================================================\n\n' ) ;


%% Globals Variables -- yuck :-(

  global share W PX x IV rc expmeanval
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
      error( 'GH Quad unsupported' ) ;
    case 2
      % Use monomial rule
      bUseMonomial = 1 ;
      error( 'Monomial Rule unsupported' ) ;
    case 3
      % Use qMC
      bUseQuasiMC = 1 ;
    case 4
      % Use Sparse Grids
      bUseSparseGrids = 1 ;
      error( 'Sparse Grids unsupported' ) ;
    otherwise
      % Use pseudo-Monte Carlo
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
  
  prods = BLPConfig.prods ;                                         % # products
  T     = BLPConfig.T ;                                         % # Markets

  starts = BLPConfig.nStarts ;                                         % # random start values to check during estimation

  
  % Compute nodes for inverting shares for starting value
  nn = nMCDraws ;       % BLPConfig.N_DRAWS_MC_INT ;
  oo = ones(1,nn);                    % Redefine index for use in simulation of shares
  Q_WEIGHTS = ones( nn, 1 ) / nn ;

  % resimulate the shock, because the econometrician does not observe it.
  v = randn( length( betatrue ), nn ) ;                     % draws for share integrals during estimation

%% Back to Dube's code ...


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2SLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xhat = iv * inv( iv' * iv ) * iv' * x ;
PX2sls = inv(xhat'*xhat)*xhat';                     % project. matrix on weighted x-space for 2SLS
beta2sls = PX2sls*y;
se2sls = sqrt(diag( mean((y-x*beta2sls).^2)*inv(xhat'*xhat) ));

expmeanval0 = exp(y);
IV = [ones(T*prods,1) A z A.^2 A.^3 z.^2 z.^3 prod(A,2) prod(z,2) kron(A(:,1),ones(1,size(z,2))).*z  kron(A(:,2),ones(1,size(z,2))).*z];
W = inv(IV'*IV);
PX = inv(x'*IV*W*IV'*x)*x'*IV*W*IV';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPEC ESTIMATION OF RANDOM COEFFICIENTS LOGIT MODEL %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta20 = 0.5*abs(beta2sls);
startvalues = repmat(theta20',starts,1).* cat(2,ones(size(theta20)),rand([starts-1,size(theta20',2)] )' *1 )';

nIV = size(IV,2);       % # instrumental variables

GMPEC = 1.0e20;
CPUtMPEC = 0;
FuncEvalMPEC = 0;
GradEvalMPEC = 0;
expmeanval = expmeanval0;

mTheta1 = zeros( size( PX, 1 ), starts ) ;
mTheta2 = zeros( length( startvalues( 1, : ) ), starts ) ;
nFree   = sum( size( IV ) ) + size( PX, 1 ) + length( startvalues( 1, : ) ) ;
mXOpt   = zeros( nFree, starts ) ;

%%% NOTE: START VALUES SET TO GIVE SAME INITIAL GMM OBJECTIVE FUNCTIONAL
%%% VALUE AS WITH NFP
theta20 = startvalues( 1, : )';
delta0 = invertshares( theta20 );
theta10 = PX*delta0;
resid0 = delta0 - x*theta10;
g0 = IV'*resid0;
x0 = [theta10; theta20; delta0; g0];
    
    
for reps=1:starts
  
    %% Create Quadrature Nodes and Weights

    % XXX Rework this as a switch statement sometime
    if bUseQuasiMC
      disp( '***Using qMC' ) ;

      nPoints = nMCDraws ;
      nBurnIn = 10000 ;
      nDim    = size( v, 1 ) ;

      mNied = NiedNormQuadInit( nPoints, nDim, nBurnIn ) ;

      v  = mNied ;
      nn = size( v, 2 ) ;
      oo = ones( 1, nn ) ;
      Q_WEIGHTS = ones( nn, 1 ) / nn ;

    else
      % configure number of simulation draws for solution of MPEC
      disp( '***Using MC' ) ;

    % %   global Q_WEIGHTS ;
      nn = nMCDraws ;       % BLPConfig.N_DRAWS_MC_INT ;
      oo = ones(1,nn);                    % Redefine index for use in simulation of shares
      Q_WEIGHTS = ones( nn, 1 ) / nn ;

      % resimulate the shock, because the econometrician does not observe it.
      v = randn( length( betatrue ), nn ) ;                     % draws for share integrals during estimation
    end


    %% Compute GMM Estimate
    
    [theta1MPEC_rep,theta2MPEC_rep,GMPEC_rep,INFOMPEC_rep,CPUtMPEC_rep,FuncEvalMPEC_rep,GradEvalMPEC_rep,xOpt] = runGMMMPECKnitro(x0,1) ;
    
    CPUtMPEC = CPUtMPEC + CPUtMPEC_rep;
    FuncEvalMPEC = FuncEvalMPEC + FuncEvalMPEC_rep;
    GradEvalMPEC = GradEvalMPEC + GradEvalMPEC_rep;
    
    mTheta1( :, reps ) = theta1MPEC_rep ;
    mTheta2( :, reps ) = theta2MPEC_rep ;
    mXOpt( :, reps )   = xOpt ;
    
    % Store best solution found via MPEC
    if (GMPEC_rep < GMPEC && ( INFOMPEC_rep == 0 || INFOMPEC_rep == 1 ))    % May need to add INFOMPEC_rep == 1 for SNOPT...
        thetaMPEC1 = theta1MPEC_rep;
        thetaMPEC2 = theta2MPEC_rep;
        GMPEC = GMPEC_rep;
        INFOMPEC = INFOMPEC_rep;
    end
end


% COVARIANCE MATRIX FOR MPEC STRUCTURAL PARAMETERS
if exist( 'thetaMPEC1', 'var' )
  delta = invertshares(thetaMPEC2);    % mean utilities
  resid = delta - x*thetaMPEC1;  % xi
  Ddelta = jacobian(thetaMPEC2);       % jacobian matrix of mean utilities
  covg = zeros(size(IV,2));
  for ii =1:length(IV),
      covg = covg + IV(ii,:)'*IV(ii,:)*(resid(ii)^2);
  end
  Dg = [x Ddelta]'*IV;            % gradients of moment conditions
  covMPEC = inv( Dg*W*Dg')*Dg*W*covg*W*Dg'*inv( Dg*W*Dg');

  cond( covMPEC )

  eig( covMPEC ) 
  

  % Save results
% %   szOutFile = sprintf( '%s/PointEst-QuadType%1dN%05d.mat', szDataDir, nQuadType, nn ) ;
% %   save( szOutFile, 'share', 'W', 'PX', 'x', 'IV', 'rc', 'expmeanval',                         ...
% %         'datasort', 'denomexpand', 'sharesum', 'oo', 'T', 'prods', 'nn', 'K',                 ...
% %         'v', 'rctrue', 'oo1', 'tol_inner', 'tol_outer', 'mXOpt', 'covMPEC',                   ...
% %         'betatrue', 'thetatrue', 'covrc', 'thetaMPEC1', 'thetaMPEC2' ) ;
    
else
  
  fprintf( 2, '\n\n********************\nDANGER: None of random starts solved!\n\n********************\n' ) ;
  
% %   szOutFile = sprintf( '%s/PointEst-QuadType%1dN%05d.mat', szDataDir, nQuadType, nn ) ;
% %   save( szOutFile, 'share', 'W', 'PX', 'x', 'IV', 'rc', 'expmeanval',                         ...
% %         'datasort', 'denomexpand', 'sharesum', 'oo', 'T', 'prods', 'nn', 'K',                 ...
% %         'v', 'rctrue', 'oo1', 'tol_inner', 'tol_outer', 'mXOpt',                              ...
% %         'betatrue', 'thetatrue', 'covrc' ) ;
end

disp( 'Finished computing GMM Estimates for BLP with Monte Carlo Data' ) ;

diary off;


