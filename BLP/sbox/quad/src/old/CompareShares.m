% CompareShares - compares different methods for computing market share integral
%
% CompareShares computes market shares for a variety of starting values and
% integration methods (simulation, Gauss-Hermite quadrature) to compare
% the results.  nParamValues is used for sensitivity testing and specifies
% how many other paramter values to generate and then use to compute share
% integrals.
%
% This code is shamelessly plundered from JP Dube.
%
% modification history
% --------------------
% 06sep2010 bss resolved plotting of two figures on one figure.
% 04sep2010 bss removed qMC.
% 19aug2010 bss added computation of max ave abs residual.
% 05may2010 bss improved configurability.
% 15apr2010 bss replaced variance with std dev to use comparable dimensions.
%               Moved `hold on' until after first plot because otherwise
%               plots were rubbish in 64-bit. Converted log to log10.
% 28mar2010 bss computation of own-price elasticity.
% 29jan2010 bss added support for Sparse Grids
% 25jan2010 bss modified to support multiple MC datasets
% 20jan2010 bss added Niederreiter QMC integration.
% 01dec2009 bss modified to take all settings from variables loaded from main.m.
% 27oct2009 bss written.
%

%
% szDataDir     - directory containing Monte Carlo dataset
% nQuadType     - quadrature type used by solver (0 == pMC)
% nSolMCDraws   - number of draws used for pMC integration with solver
% nMCDraws      - how many draws to use per simulation for computing shares via MPEC/GMM
% nQMCDraws     - number of qMC draws for integration
% nQMCBurnIn    - number of points to use for burn in when generating qMC points
% nParamValues  - number of parameter values at which to calculate shares (optional)
%

function [ ] = CompareShares( szRootDataDir, szDataDir, nQuadType, nSolMCDraws, nMCDraws, nQMCDraws, nQMCBurnIn, nParamValues ) 

  % Handle optional parameter to specify number of parameters at which to
  % evaluate share integrals.
  if nargin ~= 8
    nParamValues = 1 ;
  end
  
  szDate      = date ;
  szCurDate   = [ szDate(8:end) szDate(4:6) szDate(1:2) ] ;
  szDiaryName = sprintf( '%s/ShareIntegrals-%s-QuadType%1dN%05d.log', szDataDir, szCurDate, nQuadType, nQMCDraws ) ;

  diary( szDiaryName ) ;
  
  
  disp( 'Computing Share Integrals for BLP with Monte Carlo Data' ) ;
  
  fprintf( 1, '\n============================================================\n\n' ) ;
  fprintf( 1, '\tData Directory    : %s\n', szDataDir ) ;
  fprintf( 1, '\tQuadrature Type   : %d\n', nQuadType ) ;
  fprintf( 1, '\tNumber MC Draws   : %d\n', nMCDraws ) ;
  fprintf( 1, '\tNumber qMC Draws  : %d\n', nQMCDraws ) ;
  fprintf( 1, '\tNumber qMC Burn In: %d\n', nQMCBurnIn ) ;  
  fprintf( 1, '\tNumber of Params  : %d\n', nParamValues ) ;  
  fprintf( 1, '\n============================================================\n\n' ) ;
  

%% SETTINGS

  global share W PX x IV rc expmeanval
  global datasort denomexpand sharesum oo
  global T prods nn K v rctrue oo1 tol_inner tol_outer
  

  %% Load Data
  
%   szBLPConfigData = './Data/BLPConfig.mat' ;
%   load( szBLPConfigData ) ;   % Load general configuration information
%   load( strcat( szDataDir, '/MonteCarloBLPData.mat' ) ) ;         % Load data for a specific simulated dataset

  szBLPConfigData = sprintf( '%s/BLPConfig.mat', szRootDataDir ) ;
  load( szBLPConfigData ) ;   % Load general configuration information
  load( strcat( szDataDir, '/MonteCarloBLPData.mat' ) ) ;         % Load data for a specific simulated dataset
  
%   szPointEstFile = sprintf( '%s/PointEst_Quad%02_dN%05d.mat', szDataDir,
% %   nQuadType, nSolMCDraws ) ;
  szPointEstFile = sprintf( '%s/PointEst-QuadType%1dN%05d.mat', szDataDir, nQuadType, nSolMCDraws ) ;
  load( szPointEstFile ) ;
  
  % Files to store figures and market shares
  szFigFile   = sprintf( '%s/Figures-Quad%dN%05d.ps', szDataDir, nQuadType, nMCDraws ) ;
  szShareFile = sprintf( '%s/ShareIntegrals-Quad%dN%05d.mat', szDataDir, nQuadType, nSolMCDraws ) ;
  
  % Set seeds
  randn( 'seed', BLPConfig.vNormSeeds( 3, MySeedId ) ) ;
  rand( 'seed', BLPConfig.vUnifSeeds( 3, MySeedId ) ) ;
 
  
  nn = nMCDraws ;                                           % # draws to simulate shares

  tol_inner = BLPConfig.tol_inner ;
  tol_outer = BLPConfig.tol_outer ;
  
  prods = BLPConfig.prods ;                                         % # products
  T     = BLPConfig.T ;                                         % # Markets

  
  N_DRAWS_MC_INT   = nMCDraws ;      % How many draws to simulate the integral
  N_MC_INT         = 100 ;        % How many different MC draws to try for integral

  % % szPointEstimates = 'xOptN1000.mat' ;    % Where point estimates xOptN1000 are stored
  % % szPointEstimates = './Data/SeedB/mpecVarsN5000.mat' ;


  % aa=maxNumCompThreads(1);                            % disable multi-thread processing b/c Knitro cannot exploit this so it slows the code down


  randn('seed',155)                                   % reset normal random number generator
  rand('seed',155)                                    % reset uniform random number generator


  %% Compute Standard Errors, variance-covariance matrix 

  % Use best point estimate 
  % % xOpt = xOptN1000( :, 1 ) ;
  xOpt = mXOpt( :, 3 ) ;

  % Used only to invert shares so that all integrals are calculated with
  % the same product-market shock xi.
  global Q_WEIGHTS ;
  Q_WEIGHTS = ones( nn, 1 ) / nn ;


  %  Uses the shocks which generated the data

  nTheta     = length( betatrue ) ;
  thetaMPEC1 = xOpt( 1 : nTheta ) ;
  thetaMPEC2 = xOpt( nTheta + 1 : 2 * nTheta ) ;

  bComputeStdErrors = 0 ;
  if bComputeStdErrors
    % COVARIANCE MATRIX FOR MPEC STRUCTURAL PARAMETERS
    delta = invertshares(thetaMPEC2);    % mean utilities
  %   delta  = xOptN100( 2 * K + 3 : 2 * K + 2 + T * prods, 1 ) ;
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
  end

  %% Compare simulation and quadrature

    % Configure various integration options
    
    % Number of nodes for Gauss-Hermite Product Rule
%     vGHNodes = [ 3 ; 4 ; 5 ; 7 ; 9 ] ;      % Gauss-Hermite node spacings to try
    vGHNodes = [ 3 ; 4 ; 5 ; 7 ] ; 
%     vGHNodes = [ 7 ] ;
    nGHTypes = length( vGHNodes ) ;

    % Allocate storage for computed shares
    mSharesSim = zeros( prods * T, N_MC_INT, nParamValues ) ;
    
    mSharesNied = zeros( prods * T, N_MC_INT, nParamValues ) ;
    
    mSharesGH  = zeros( prods * T, nGHTypes, nParamValues ) ;
    mSharesStroud = zeros( prods * T, 2, nParamValues ) ;
% %     mResidStroud  = zeros( prods * T, 2, nParamValues ) ;
    mSharesSparse = zeros( prods * T, 1, nParamValues ) ; 
    
    mOwnElastSim  = mSharesSim ;
    mOwnSparse    = mSharesSparse ;
    
    % Generate different parameter values.  In all cases, we use the same
    % unobserved product-market shocks xi_jt, only varying theta1 & theta2
    %
    % xOpt( 1 : nTheta )                      - theta1 ( params on
    %                                           characteristics and price )
    % xOpt( nTheta + 1 : 2 * nTheta )         - theta2 ( variance for RCs )
    % xOpt( 2 * nTheta + ( 1 : T * prods ) )  - delta ( [ -price  x ]'theta1 + xi_jt )
    % xOpt( 2 * nTheta + T * prods + ( 1 : T * prods ) )  - g = Z'xi_jt
    
    dwParamStdErr = 0.25 ;     % Standard Error for perturbing parameter values
    mParams = repmat( xOpt( 1 : 2 * nTheta ), [ 1, nParamValues ] )         ...
            .* ( ones( 2 * nTheta, nParamValues ) + dwParamStdErr * randn( 2 * nTheta, nParamValues ) ) ;
    
    % Ensure first item is original parameter values
    mParams( :, 1 ) = xOpt( 1 : 2 * nTheta ) ;
  
  for ixParam = 1 : nParamValues
    
    % Get perturbed parameter value at which to evaluate market share
    % integrals
    xParam = xOpt ;  % XXX generate random paramter values about xOpt
    xParam( 1 : 2 * nTheta ) = mParams( :, ixParam ) ;
    
    
    %% Simulation

    bCalcSim = 1 ;

    if bCalcSim

      disp( 'Using pseudo-Monte Carlo to compute shares and residuals' ) ;

  % %     mSharesSim = zeros( prods * T, N_MC_INT, nParamValues ) ;
  % %     mResidSim  = zeros( prods * T, N_MC_INT ) ;
  % %     vIter      = zeros( N_MC_INT ) ;      % number of iterations for Berry's mapping to converge

      for ix = 1 : N_MC_INT

        % resimulate the shock, because the econometrician does not observe it.
        v_sim = randn( length( betatrue ), N_DRAWS_MC_INT ) ;                     % draws for share integrals during estimation

        v  = v_sim ;
        nn = size( v, 2 ) ;
        oo = ones( 1, nn ) ; 
        Q_WEIGHTS = ones( nn, 1 ) / nn ;

        % Computer shares via Monte Carlo integration
        mSharesSim( :, ix, ixParam ) = ComputeShares( xParam ) ;   % Market shares, sorted by product, then market

        mOwnSim( :, ix, ixParam )    = ComputeOwnPriceElasticity( xParam ) ;
        
% %         [ mResidSim( :, ix ), vIter( ix ) ]  = ComputeShock( xParam ) ;    % Compute xi_jt
        
        %% Check conditions for contraction mapping theorem to be satisfied
        bCheckMap = 0 ;        
        if bCheckMap
          CheckMap( xOpt ) ;          
        end
        
        
      end 

    end

    % mJ    = CSDJacobian( @(x) ComputeMap( x, xParam ), mResidSim(:,1), 1e-5 ) ;
    % vEigs = eigs( mJ, 5, 'sm' ) ;       % 5 smallest eigenvectors
    % cond( full( mJ ) )                  % Condition number
    % hist( log10( abs( eig( full( mJ ) ) ) ), 10 ) ;



    %% Quasi-Monte Carlo

    bUseQMC = 0 ;
    if bUseQMC

      disp( 'Using quasi-Monte Carlo to compute shares and residuals' ) ;

      nPoints = nQMCDraws ;     % 1000
      nBurnIn = nQMCBurnIn ;    % 1000
      nDim    = 5 ;

      mNied = NiedNormQuadInit( nPoints * N_MC_INT, nDim, nBurnIn ) ;

    % %   mSharesNied = zeros( prods * T, N_MC_INT, nParamValues ) ;

      ixPos = 0 ;
      for ix = 1 : nPoints : ( N_MC_INT * nPoints )

        ixPos = ixPos + 1 ;
        v  = mNied( :, ix : ix + nPoints - 1 ) ;
        nn = size( v, 2 ) ;
        oo = ones( 1, nn ) ;
        Q_WEIGHTS = ones( nn, 1 ) / nn ;

        mSharesNied( :, ixPos, ixParam ) = ComputeShares( xParam ) ;
      end

    end


    %% GH Quadrature with Product Rule

    bCalcGH = 1 ;
    if bCalcGH

      disp( 'Using GH Product Rule to compute shares and residuals' ) ;

  % %   %   vGHNodes = [ 3 ; 4 ; 5 ; 7 ; 9 ] ;      % Gauss-Hermite node spacings to try
  % %     vGHNodes = [ 7 ] ;
  % %     nGHTypes = length( vGHNodes ) ;
    % %   covrc = diag( [ .5 .5 .5 .5 .2] );                  % true co-variances in tastes
      theta2 = xParam( K + 2 : 2 * K + 2, 1 ) ;                  % st. deviation of tastes
      covrc = diag( theta2 ) ;

  % %     mSharesGH  = zeros( prods * T, nGHTypes, nParamValues ) ;

      for jx = 1 : nGHTypes

        % Compute Shares via Gaussian-Hermite Quadrature
        GHQuadInit( betatrue, covrc, vGHNodes( jx ) ) ;  % Change the last argument to change the number of nodes

        global Q_NODES ;

        bUseHatchet = 0 ;
        if bUseHatchet
          % now remove nodes whose weights are less than 1e-10

          nNodes    = length( Q_WEIGHTS ) ;
          vBig      = find( Q_WEIGHTS >= 1e-8 ) ;

          fprintf( 1, '\tFraction of nodes remaining after cutoff for grid spacing %d^5: %6.2f (Total Nodes = %10d)\n',     ...
                   vGHNodes( jx ), length( vBig ) / nNodes, length( vBig ) ) ;

          Q_WEIGHTS = Q_WEIGHTS( vBig ) ;
          Q_NODES   = Q_NODES( :, vBig ) ;
        end

        v  = Q_NODES ;    % Use Gaussian-Hermite nodes

        nn = size( v, 2 ) ;
        oo = ones( 1, nn ) ; 

        mSharesGH( :, jx, ixParam ) = ComputeShares( xParam ) ;

      %   for ix = 1 : T
      %     vIx = ix + ( 0 : 50 : ( T - 1 ) * prods ) ;
      %     totGH  = sum( mSharesGH( vIx, jx ) ) ;
      %     totSim = sum( mSharesSim( vIx, jx ) ) ;
      %     fprintf( 1, 'Market = %4d : (GH SIM diff): %f %f %15.7e\n', ix, totGH, totSim, ( totGH - totSim ) / totGH ) ;
      %   end

      end

    end


    %% Stroud Monomial Rule

    disp( 'Using Monomial Rule to compute shares and residuals' ) ;

    % Compute Shares via Stroud's 11-1 monomial rule
  % %   mSharesStroud = zeros( prods * T, 2, nParamValues ) ;
  % %   mResidStroud  = zeros( prods * T, 2, nParamValues ) ;

    global Q_NODES ;

    for ixRule = 1 : 2
      MonoQuadInit( ixRule ) ;

      v  = Q_NODES ;    % Use Gaussian-Hermite nodes

      nn = size( v, 2 ) ;
      oo = ones( 1, nn ) ; 

      mSharesStroud( :, ixRule, ixParam ) = ComputeShares( xParam ) ;

    %   [ mResidStroud( :, ixRule ), nIter ]  = ComputeShock( xParam ) ;    % Compute xi_jt
    end


    %% Sparse Grids

    bUseSparseGrids = 1 ;

    if bUseSparseGrids

      disp( 'Using Sparse Grids to compute shares and residuals' ) ;

      % Make sure you download matlab code from www.sparse-grids.de
      % Then put the code in your path

      addpath( '~/Tools/SparseGrids' ) ;    % Whereever the sparse grids code is installed

      % Konrad-Patterson Sparse Grid for Normal Kernel in 5 dim with K = 6
      % (so that it will be exact to same degree as monomial rule)
      % K = 7 => accuracy for total order 2*K-1
      [ Q_NODES, Q_WEIGHTS ] = nwspgr( 'KPN', 5, 6 ) ;      

  % %     mSharesSparse = zeros( prods * T, 1, nParamValues ) ; 

      v  = Q_NODES' ;
      nn = size( v, 2 ) ;
      oo = ones( 1, nn ) ;

      mSharesSparse( :, ixParam ) = ComputeShares( xParam ) ;
      
      mOwnSparse( :, ixParam )    = ComputeOwnPriceElasticity( xParam ) ;
      
% %       [ vResidSparse, nIterSparse ]  = ComputeShock( xParam ) ;    % Compute xi_jt

    end
    
  end     % End of computing shares for a variety of parameter values
  
  
  %% Save computed market share data
  
  % Create dummies for variables which may not have been computed but I 
  % want to save
  
  if ~exist( 'vIter' )
    vIter = NaN ;
  end
  
  if ~exist( 'mResidSim' )
    mResidSim = NaN ;
  end
  
  if ~exist( 'vResidSparse' )
    vResidSparse = NaN ;
  end
  
  if ~exist( 'nIterSparse' )
    nIterSparse = NaN ;
  end
  
  if ~exist( 'mSharesNied' )
    mSharesNied = NaN ;
  end
  
  save( szShareFile, 'mSharesSim', 'mSharesNied', 'mSharesGH', 'mSharesStroud',     ...
                     'mSharesSparse', 'xOpt', 'mParams', 'mOwnSparse', 'mOwnSim',   ...
                      'vIter', 'mResidSim', 'vResidSparse', 'nIterSparse' ) ;                  
  
                    
                    
  %% Compare elasticities
  
  mTmp  = repmat( mOwnSparse( :, 1 ), 1, N_MC_INT ) ;
  
  mDiff = 1 - mOwnSim( :, :, 1 ) ./ mTmp ;
  
% %   figure( 222 ) ;
% %   set( 222, 'Visible', 'off' ) ;
  h = newplot() ;
  
  hist( mDiff(:), 50 ) ;
  title( 'Histogram of Deviations of Own Price Elasticity (MC vs SGI)' ) ;
  xlabel( 'Percent Deviation' ) ;
  ylabel( 'Frequency' ) ;
  
  saveas( gcf, sprintf( '%s/Fig-DevOwnPriceElasticityN%05d.jpg', szDataDir, nMCDraws ), 'jpg' ) ;
  
  % Compute number with more than x% deviation
  dwCutoff = 0.01 ;   % percent deviation
  fprintf( 1, '\n\t\tNumber of deviations of more than %d%%: %d of %d observations\n\n', dwCutoff * 100,      ...
              sum( abs( mDiff(:) >= dwCutoff ) ), length( mDiff(:) ) ) ;
  
            
  
  %% Compute Residuals                    
  
  bComputeResiduals = 1 ;
  if bComputeResiduals
    
    jxTrial = 1 ;   % Which parameter value to use
    vRefShares = mSharesGH( :, ( find( 7 == vGHNodes ) ), jxTrial ) ;
    
    fprintf( 1, '\n\n============================================================\n\n' ) ;
    fprintf( 1, '\t\tTrial : %d\n', jxTrial ) ;
    fprintf( 1, '\tResiduals for pMC %d draws and qMC %d draws\n', N_DRAWS_MC_INT, nQMCDraws ) ;
    fprintf( 1, '\t\t\t\t\t\t(max abs residual) (min abs residual) (mean abs residual)\n' ) ;
    fprintf( 1, '\tpMC\t\t\t\t\t\t: %10.5e %10.5e %10.5e\n', max( max( abs( mSharesSim( :, :, jxTrial ) - repmat( vRefShares, 1, N_MC_INT ) ) ) ),    ...
                                            min( min( abs( mSharesSim( :, :, jxTrial ) - repmat( vRefShares, 1, N_MC_INT ) ) ) ),                     ...
                                            mean( mean( abs( mSharesSim( :, :, jxTrial ) - repmat( vRefShares, 1, N_MC_INT ) ) ) ) ) ;

    if bUseQMC                                          
      fprintf( 1, '\tqMC\t\t\t\t\t\t: %10.5e %10.5e %10.5e\n', max( max( abs( mSharesNied( :, :, jxTrial ) - repmat( vRefShares, 1, N_MC_INT ) ) ) ),   ...
                                            min( min( abs( mSharesNied( :, :, jxTrial ) - repmat( vRefShares, 1, N_MC_INT ) ) ) ),                    ...
                                            mean( mean( abs( mSharesNied( :, :, jxTrial ) - repmat( vRefShares, 1, N_MC_INT ) ) ) ) ) ;
    end
    
    for ix = 1 : length( vGHNodes )                                      
      fprintf( 1, '\tGH( %d)\t\t\t\t: %10.5e %10.5e %10.5e\n', vGHNodes( ix ),           ...
                                              ( max( abs( mSharesGH( :, ix, jxTrial ) - vRefShares ) ) ),   ...
                                              ( min( abs( mSharesGH( :, ix, jxTrial ) - vRefShares ) ) ),   ...
                                              ( mean( abs( mSharesGH( :, ix, jxTrial ) - vRefShares ) ) ) ) ;
    end       
    
    fprintf( 1, '\tStroud(Left)\t: %10.5e %10.5e %10.5e\n',           ...
                                        ( max( abs( mSharesStroud( :, 1, jxTrial ) - vRefShares ) ) ),    ...
                                        ( min( abs( mSharesStroud( :, 1, jxTrial ) - vRefShares ) ) ),    ...
                                        ( mean( abs( mSharesStroud( :, 1, jxTrial ) - vRefShares ) ) ) ) ;


    fprintf( 1, '\tStroud(Right)\t: %10.5e %10.5e %10.5e\n',           ...
                                        max( max( abs( mSharesStroud( :, 2, jxTrial ) - vRefShares ) ) ),   ...
                                        min( min( abs( mSharesStroud( :, 2, jxTrial ) - vRefShares ) ) ),   ...
                                        mean( mean( abs( mSharesStroud( :, 2, jxTrial ) - vRefShares ) ) ) ) ; 
                                      

    fprintf( 1, '\tSparse\t\t\t\t: %10.5e %10.5e %10.5e\n',           ...
                                        max( max( abs( mSharesSparse( :, 1, jxTrial ) - vRefShares ) ) ),   ...
                                        min( min( abs( mSharesSparse( :, 1, jxTrial ) - vRefShares ) ) ),   ...
                                        mean( mean( abs( mSharesSparse( :, 1, jxTrial ) - vRefShares ) ) ) ) ;  
  
  end
    
  
  bComputeRelDiff = 1 ;
  if bComputeRelDiff
    
    jxTrial = 1 ;   % Which parameter value to use
    vRefShares = mSharesGH( :, ( find( 7 == vGHNodes ) ), jxTrial ) ;
    
    fprintf( 1, '\n\n============================================================\n\n' ) ;
    fprintf( 1, '\t\tTrial : %d\n', jxTrial ) ;
    fprintf( 1, '\tPercentage Difference for pMC %d draws and qMC %d draws\n', N_DRAWS_MC_INT, nQMCDraws ) ;
    fprintf( 1, '\t\t\t\t\t\t(max abs residual) (min abs residual)\n' ) ;
    fprintf( 1, '\tpMC\t\t\t\t\t\t: %10.5e %10.5e\n', max( max( abs( 1 - mSharesSim( :, :, jxTrial ) ./ repmat( vRefShares, 1, N_MC_INT ) ) ) ),   ...
                                            min( min( abs( 1 - mSharesSim( :, :, jxTrial ) ./ repmat( vRefShares, 1, N_MC_INT ) ) ) ) ) ;
     
    if bUseQMC                                          
      fprintf( 1, '\tqMC\t\t\t\t\t\t: %10.5e %10.5e\n', max( max( abs( 1 - mSharesNied( :, :, jxTrial ) ./ repmat( vRefShares, 1, N_MC_INT ) ) ) ),   ...
                                            min( min( abs( 1 - mSharesNied( :, :, jxTrial ) ./ repmat( vRefShares, 1, N_MC_INT ) ) ) ) ) ;
    end
    
    for ix = 1 : length( vGHNodes )                                      
      fprintf( 1, '\tGH( %d)\t\t\t\t: %10.5e %10.5e\n', vGHNodes( ix ),           ...
                                              ( max( abs( 1 - mSharesGH( :, ix, jxTrial ) ./ vRefShares ) ) ),   ...
                                              ( min( abs( 1 - mSharesGH( :, ix, jxTrial ) ./ vRefShares ) ) ) ) ;
    end       
    
    fprintf( 1, '\tStroud(Left)\t: %10.5e %10.5e\n',           ...
                                        ( max( abs( 1 - mSharesStroud( :, 1, jxTrial ) ./ vRefShares ) ) ),   ...
                                        ( min( abs( 1 - mSharesStroud( :, 1, jxTrial ) ./ vRefShares ) ) ) ) ;


    fprintf( 1, '\tStroud(Right)\t: %10.5e %10.5e\n',           ...
                                        max( max( abs( 1 - mSharesStroud( :, 2, jxTrial ) ./ vRefShares ) ) ),   ...
                                        min( min( abs( 1 - mSharesStroud( :, 2, jxTrial ) ./ vRefShares ) ) ) ) ; 
                                      

    fprintf( 1, '\tSparse\t\t\t\t: %10.5e %10.5e\n',           ...
                                        max( max( abs( 1 - mSharesSparse( :, 1, jxTrial ) ./ vRefShares ) ) ),   ...
                                        min( min( abs( 1 - mSharesSparse( :, 1, jxTrial ) ./ vRefShares ) ) ) ) ;  
  
  end

  
  %% Compute statistics for simulated shares

  vMeans   = mean( mSharesSim( :, :, 1 ), 2 ) ;
  vMedians = median( mSharesSim( :, :, 1 ), 2 ) ;
  vStdDev  = std( mSharesSim( :, :, 1 ), 1, 2 ) ;
  mStdDev  = repmat( vStdDev, 1, N_MC_INT ) ;


  %% Make nice plots
  
  bVisible = 'off' ;   % Whether to display plots or only save to a file { on|off }

  % Create Index to grab all values for xParam == xOpt, i.e. one parameter
  % value.
  
  jxParam       = 1 ;    % Select which parameter to plot
  vIxSingleVec  = ( jxParam - 1 ) * T * prods + ( 1 : T * prods )' ;
  vIxStroud     = ( jxParam - 1 ) * T * prods * size( mSharesStroud, 2 ) + ( 1 : T * prods )' ;
  
  nGHOffset     = ( find( 7 == vGHNodes ) - 1 ) * T * prods ;               % Offset to get to the 7^5 rule
  vIxGH         = ( jxParam - 1 ) * T * prods * size( mSharesGH, 2 ) + ( 1 : T * prods )'         ...
                    + nGHOffset ;
  vIxSingleMat  = ( jxParam - 1 ) * T * prods * N_MC_INT + ( 1 : T * prods * N_MC_INT )' ;
  
  % Market Share vs. Variance( s_i )
% %   figure( 900 ) ;
% %   set( 900, 'Visible', bVisible ) ;
  h = newplot() ;
  hold on ;

  plot( mSharesSim( vIxSingleMat ), mStdDev(:), 'go' ) ;
    
% %   hold on ;
  
  if bUseQMC
    plot( mSharesNied( vIxSingleMat ), mStdDev(:), 'kh' ) ;
  end
  
  plot( mSharesStroud( vIxStroud ), vStdDev, 'r<' ) ;               % Left col of Stroud
  plot( mSharesStroud( vIxStroud + T * prods ), vStdDev, 'b>' ) ;   % Right Col of Stroud
  plot( mSharesGH( vIxGH ), vStdDev, 'mp' ) ;     % XXX uses first GH set of nodes
  plot( mSharesSparse( vIxSingleVec ), vStdDev, 'yd' ) ;     % Sparse Grids

  szNumDraws = sprintf(  'N_{Sim} = %d, R_{Draws} = %d', N_MC_INT, N_DRAWS_MC_INT ) ;
  % title( { 'Market Share vs. StdDev [ Market Share ]', szNumDraws } ) ;
  xlabel( 'Market Share ( s_{i} )' ) ;
  ylabel( 'Std. Err.[ pMC Market Share ]' ) ;

  if bUseQMC
    legend( 'Monte Carlo', 'Niederreiter', 'Stroud 11-1 Left', 'Stroud 11-1 Right',     ...
          'Product Rule (7^5 Nodes)', 'Sparse Grids' ) ;
  else
    legend( 'Monte Carlo', 'Stroud 11-1 Left', 'Stroud 11-1 Right',     ...
          'Product Rule (7^5 Nodes)', 'Sparse Grids' ) ;
  end
  
  hold off ;
  
  xlim( [ 0 1 ] ) ;     % Stroud Rule 11-1 Right produces six negative market shares < o(1e-9) so I nuke them

  print( '-dpsc2', szFigFile ) ;
  
  saveas( gcf, sprintf( '%s/Fig-SharesVsStdDevN%05d.new.jpg', szDataDir, nMCDraws ), 'jpg' ) ;
  
  
  % close up for one share
  
  [ junk, ixMaxShare ] = max( mSharesSparse ) ;
  
% %   figure( 110 ) ;
% %   set( 110, 'Visible', bVisible ) ;
  
  h = newplot() ;
  hold on ;
  
  plot( ones( N_MC_INT, 1 ), mSharesSim( ixMaxShare, : ) , 'go' ) ;
    
% %   hold on ;
  
  if bUseQMC 
    plot( 2 * ones( N_MC_INT, 1 ), mSharesNied( ixMaxShare, : ), 'kh' ) ;
    nQMCOffset = 0 ;
  else
    nQMCOffset = 1 ;
  end
  
  plot( 3 - nQMCOffset, mSharesStroud( ixMaxShare, 1 ), 'r<' ) ;               % Left col of Stroud
  plot( 4 - nQMCOffset, mSharesStroud( ixMaxShare, 2 ), 'b>' ) ;   % Right Col of Stroud
  plot( 5 - nQMCOffset, mSharesGH( ixMaxShare, nGHTypes ), 'mp' ) ;     % XXX uses first GH set of nodes
  plot( 6 - nQMCOffset, mSharesSparse( ixMaxShare ), 'yd' ) ;     % Sparse Grids

  szNumDraws = sprintf(  'N_{Sim} = %d, R_{Draws} = %d', N_MC_INT, N_DRAWS_MC_INT ) ;
  title( { 'Comparison of Rules for a Single Share Integral', szNumDraws } ) ;
  xlabel( 'Quadrature Rule' ) ;
  ylabel( 'Market Share ( s_{i} )' ) ;

  if bUseQMC
    legend( 'Monte Carlo', 'Niederreiter', 'Stroud 11-1 Left', 'Stroud 11-1 Right',     ...
          'Product Rule', 'Sparse Grids' ) ;
  else
    legend( 'Monte Carlo', 'Stroud 11-1 Left', 'Stroud 11-1 Right',     ...
          'Product Rule', 'Sparse Grids' ) ;    
  end
  
  hold off ;
  
  xlim( [ 0, 7 ] ) ;

  print( '-dpsc2', szFigFile ) ;
  
  saveas( gcf, sprintf( '%s/Fig-CloseUpSharesVsStdDevN%05d.new.jpg', szDataDir, nMCDraws ), 'jpg' ) ;
  
    
  % Market Share vs. Log Variance( s_i )
  
  vLogStdDev = log10( vStdDev ) ;
  mLogStdDev = log10( mStdDev ) ;
  
% %   figure( 105 ) ;
% %   set( 105, 'Visible', bVisible ) ;


  h = newplot() ;
  hold on ;
  
% %   xlim( [ 0 1 ] ) ;     % Stroud Rule 11-1 Right produces six negative market shares < o(1e-9) so I nuke them

  plot( log10( mSharesSim( vIxSingleMat ) ), mLogStdDev(:), 'go' ) ;
  
% %   hold on ;
  
  if bUseQMC
    plot( log10( mSharesNied( vIxSingleMat ) ), mLogStdDev(:), 'kh' ) ;
  end
  
  plot( log10( mSharesStroud( vIxStroud ) ), vLogStdDev, 'r<' ) ;               % Left col of Stroud
  plot( log10( mSharesStroud( vIxStroud + T * prods ) ), vLogStdDev, 'b>' ) ;   % Right Col of Stroud
  plot( log10( mSharesGH( vIxGH ) ), vLogStdDev, 'mp' ) ;     % XXX uses first GH set of nodes
  plot( log10( mSharesSparse( vIxSingleVec ) ), vLogStdDev, 'yd' ) ;     % Sparse Grids
  
% %   szNumDraws = sprintf(  'N_{Sim} = %d, R_{Draws} = %d', N_MC_INT, N_DRAWS_MC_INT ) ;
  title( { 'Log_{10}[ Market Share ] vs. Log_{10}[ Std Dev ( Market Share ) ]', szNumDraws } ) ;
  xlabel( 'Log_{10}[ Market Share ( s_{i} ) ]' ) ;
  ylabel( 'Log_{10}[ Std Dev( Market Share ) ]' ) ;

  if bUseQMC
    legend( 'Monte Carlo', 'Niederreiter', 'Stroud 11-1 Left', 'Stroud 11-1 Right',       ...
          'Product Rule (7^5 Nodes)', 'Sparse Grids', 'Location', 'SouthEast' ) ;
  else        
    legend( 'Monte Carlo', 'Stroud 11-1 Left', 'Stroud 11-1 Right',       ...
          'Product Rule (7^5 Nodes)', 'Sparse Grids', 'Location', 'SouthEast' ) ;
  end

  hold off ;

  print( '-append', '-dpsc2', szFigFile ) ;
  
  saveas( gcf, sprintf( '%s/Fig-LogSharesVsLogStdDevN%05d.jpg', szDataDir, nMCDraws ), 'jpg' ) ;
  
  
  %% Change in Elasticity vs. Share size
  
% %   figure( 600 ) ;
% %   set( 600, 'Visible', bVisible ) ;
    
  h = newplot() ;
  
  plot( mSharesGH( :, find( 7 == vGHNodes ) ), mDiff, 'bx' ) ; 
  title( 'Percentage difference in elasticity vs. computed shares' ) ;
  xlabel( 'Computed Shares (GH 5^7)' );
  ylabel( 'Percentage Difference in Elasticity' ) ;
  
  
  print( '-append', '-dpsc2', szFigFile ) ;
  
  saveas( gcf, sprintf( '%s/Fig-ElasticityVsSharesN%05d.jpg', szDataDir, nMCDraws ), 'jpg' ) ;
  
  
  %% Histogram of residuals
  
% %   figure( 777 ) ;
% %   set( 777, 'Visible', bVisible ) ;
  
  h = newplot() ;
  
  hist( log10( [ mSharesGH(:,4,1) mSharesStroud(:,:,1) mSharesSparse(:,1) ] ) - repmat( log10( vMeans ), 1, 4), 25) ;
  legend( 'Product Rule (7^5 Nodes)', 'Stroud 11-1 Left', 'Stroud 11-1 Right', 'Sparse Grids' ) ;
  xlabel( 'Market Share - Mean of Simulated Market Share' ) ;
  ylabel( 'Count' ) ;
  title( 'Histogram of Market Share Residuals - Polynomial Rules vs. Simulation' )
  
  print( '-append', '-dpsc2', szFigFile ) ;
    
  saveas( gcf, sprintf( '%s/Fig-HistResidualN%05d.jpg', szDataDir, nMCDraws ), 'jpg' ) ;
  
  %% histogram of market shares

% %   figure( 200 ) ;
% %   set( 200, 'Visible', bVisible ) ;
    
  h = newplot() ;
  
  hist( real( log10( [ mSharesGH( :, end - 1, 1 ), mSharesStroud( :, :, 1 ), mSharesSparse( :, 1 ), vMeans, vMedians ] ) ), 25 ) ;
  title( { 'Histogram of Market Shares by Rule', szNumDraws } ) ;
  xlabel( 'log_{10}( s_i )' );
  ylabel( 'Count' ) ;

  legend( 'Product Rule (9^5 Nodes)', 'Stroud 11-1 Left', 'Stroud 11-1 Right', 'Sparse Grid', 'Simulated Mean( s_i )', 'Simulated Median( s_i )' ) ;

  print( '-append', '-dpsc2', szFigFile ) ;
  
  saveas( gcf, sprintf( '%s/Fig-HistLogSharesN%05d.jpg', szDataDir, nMCDraws ), 'jpg' ) ;
  
  
  %% Histogram of bias

% %   figure( 300 ) ;
% %   set( 300, 'Visible', bVisible ) ;
  
  h = newplot() ;
  
  hist( real( log10( [ mSharesGH( :, end - 1, 1 ), mSharesStroud( :, :, 1 ) ] - repmat( vMeans, 1, 3 ) ) ), 25 ) ;
  title( { 'Histogram of Difference in Market Shares from Simulation Mean by Rule', szNumDraws } ) ;
  xlabel( 'log_{10}( s_i )' );
  ylabel( 'Count' ) ;

  legend( 'Product Rule (9^5 Nodes)', 'Stroud 11-1 Left', 'Stroud 11-1 Right' ) ;

  print( '-append', '-dpsc2', szFigFile ) ;

  %% Cummulative Histogram of Market Shares

% %   figure( 400 ) ;
% %   set( 400, 'Visible', bVisible ) ;
  
  h = newplot() ;
  
  [ mCount, vBins ] = hist( real( log10( [ mSharesGH( :, end - 1, 1 ), mSharesStroud( :, :, 1 ), vMeans, vMedians ] ) ), 25 ) ;
  plot( vBins, cumsum( mCount ) ./ ( T * prods ), '.' ) ;
  title( { 'Cummulative Fraction of Market Shares by Rule', szNumDraws } ) ;
  xlabel( 'log_{10}( s_i )' );
  ylabel( 'Cummulative Fraction of Market Shares \leq s_i' ) ;

  legend( 'Product Rule (9^5 Nodes)', 'Stroud 11-1 Left', 'Stroud 11-1 Right', 'Simulated Mean( s_i )',         ...
          'Simulated Median( s_i )', 'Location', 'SouthEast' ) ;

  print( '-append', '-dpsc2', szFigFile ) ;

  saveas( gcf, sprintf( '%s/Fig-SharesCDF.fig', szDataDir ), 'fig' ) ;
  
  %% Plots vs. xi_jt heterogeneity

  vDelta = xOpt( 2 * K + 3 : 2 * K + 2 + T * prods, 1 ) ;         % mean utilities
  vShock = vDelta - x * thetaMPEC1 ;  % xi

% %   figure( 500 ) ;
% %   set( 500, 'Visible', bVisible ) ;
  
  h = newplot() ;
  
  plot( vShock, x * thetaMPEC1, 'o' ) ;
  title( '\xi_{jt} vs. X^{T}\beta' ) ;
  xlabel( '\xi_{jt}' ) ;
  ylabel( 'X^{T}\beta' ) ;

  print( '-append', '-dpsc2', szFigFile ) ;

  saveas( gcf, sprintf( '%s/Fig-XiVsIndex.fig', szDataDir ), 'fig' ) ;
  
  %% 

  % Simulation
  for ix = 1:K+1
    nu( ix, : ) = thetaMPEC2( ix ) * v_sim( ix, : ) ;
  end

  MU = x * nu ;

% %   figure( 501 ) ;
% %   set( 501, 'Visible', bVisible ) ;
  
  h = newplot() ;
  
  plot( vShock, MU, 'x' ) ;
  title( { '\mu - Simulation', szNumDraws } ) ;
  xlabel( '\xi_{jt}' ) ;
  ylabel( '\mu_{sim}' ) ;

  print( '-append', '-dpsc2', szFigFile ) ;
  
  saveas( gcf, sprintf( '%s/Fig-MuSimulation.jpg', szDataDir ), 'jpg' ) ;

  % Stroud 11-1
  
  MonoQuadInit( 1 ) ;   % Should use right?

  v  = Q_NODES ;    % Use Gaussian-Hermite nodes

  nn = size( v, 2 ) ;
  oo = ones( 1, nn ) ; 

  for ix = 1:K+1
    nu2( ix, : ) = thetaMPEC2( ix ) * v( ix, : ) ;
  end

  MU2 = x * nu2 ;

% %   figure( 502 ) ;
% %   set( 502, 'Visible', bVisible ) ;
  
  h = newplot() ;
  
  plot( vShock, MU2, 'x' ) ;
  title( { '\mu - Stroud 11-1 Left', szNumDraws } ) ;
  xlabel( '\xi_{jt}' ) ;
  ylabel( '\mu_{Stroud}' ) ;

  print( '-append', '-dpsc2', szFigFile ) ;
  
  saveas( gcf, sprintf( '%s/Fig-MuStroud.jpg', szDataDir ), 'jpg' ) ;
  
  %% Plot of Error in Shares vs. Shock
  
  
  %% Plot histogram of iterations to compute Berry Map
  
  if exist( 'mResidSim', 'var' )
    
% %     figure( 505 ) ;
% %     set( 505, 'Visible', bVisible ) ; 
    
    h = newplot() ;
    
% %     vBrks = linspace( 0, 1000, 25 ) ;
    hist( vIter, 25 ) ;
    title( sprintf( 'Iterations Until Convergence of Berry''s Mapping: R = %d', N_DRAWS_MC_INT ) ) ;
    xlabel( 'Number of Iterations' ) ;
    
    print( '-append', '-dpsc2', szFigFile ) ;
  
    saveas( gcf, sprintf( '%s/Fig-HistBerryN%05d.jpg', szDataDir, N_DRAWS_MC_INT ), 'jpg' ) ;
  end
  
  
  %% Clean up
  
  diary off;
  
  

