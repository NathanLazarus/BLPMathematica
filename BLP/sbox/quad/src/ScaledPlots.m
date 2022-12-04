% ScaledPlots - plots Fig 1 in paper but rescaled

% Make sure you load the Shares.mat file for R=1000.....

%% Compute statistics for simulated shares

  vMeans   = mean( mSharesSim( :, :, 1 ), 2 ) ;
  vMedians = median( mSharesSim( :, :, 1 ), 2 ) ;
  vStdDev  = std( mSharesSim( :, :, 1 ), 1, 2 ) ;
  mStdDev  = repmat( vStdDev, 1, N_MC_INT ) ;

  mMeans   = repmat( vMeans, 1, N_MC_INT ) ;

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

  plot( mStdDev(:), mSharesSim( vIxSingleMat ) ./ mMeans( vIxSingleMat ) - 1, 'go' ) ;
  
  plot( vStdDev, mSharesStroud( vIxStroud ) ./ vMeans - 1, 'r<' ) ;               % Left col of Stroud
  plot( vStdDev, mSharesStroud( vIxStroud + T * prods ) ./ vMeans - 1, 'b>' ) ;   % Right Col of Stroud
  plot( vStdDev, mSharesGH( vIxGH ) ./ vMeans - 1, 'mp' ) ;     % XXX uses first GH set of nodes
  plot( vStdDev, mSharesSparse( vIxSingleVec ) ./ vMeans - 1, 'yd' ) ;     % Sparse Grids

  szNumDraws = sprintf(  'N_{Sim} = %d, R_{Draws} = %d', N_MC_INT, N_DRAWS_MC_INT ) ;
  % title( { 'Market Share vs. StdDev [ Market Share ]', szNumDraws } ) ;
  ylabel( 'Relative Market Share Error [ ( s_{i} - s^{pMC}_{i} ) / s^{pMC}_{i} ]' ) ;
  xlabel( 'Std. Err.[ pMC Market Share ]' ) ;

  legend( 'Monte Carlo', 'Stroud 11-1 Left', 'Stroud 11-1 Right',     ...
        'Product Rule (7^5 Nodes)', 'Sparse Grids' ) ;

  
  hold off ;
  
  % xlim( [ 0 1 ] ) ;     % Stroud Rule 11-1 Right produces six negative market shares < o(1e-9) so I nuke them


  %% Save

  saveas( gcf, sprintf( '%s/Fig-RelSharesErrorVsStdErrN%05d.new.jpg', szDataDir, nMCDraws ), 'jpg' ) ;
  
  
  %% Versus Market Share
  
  h = newplot() ;
  hold on ;

  plot( mMeans(vIxSingleMat), mSharesSim( vIxSingleMat ) ./ mMeans( vIxSingleMat ) - 1, 'go' ) ;
  
  plot( vMeans, mSharesStroud( vIxStroud ) ./ vMeans - 1, 'r<' ) ;               % Left col of Stroud
  plot( vMeans, mSharesStroud( vIxStroud + T * prods ) ./ vMeans - 1, 'b>' ) ;   % Right Col of Stroud
  plot( vMeans, mSharesGH( vIxGH ) ./ vMeans - 1, 'mp' ) ;     % XXX uses first GH set of nodes
  plot( vMeans, mSharesSparse( vIxSingleVec ) ./ vMeans - 1, 'yd' ) ;     % Sparse Grids

  szNumDraws = sprintf(  'N_{Sim} = %d, R_{Draws} = %d', N_MC_INT, N_DRAWS_MC_INT ) ;
  % title( { 'Market Share vs. StdDev [ Market Share ]', szNumDraws } ) ;
  ylabel( 'Relative Market Share Error [ ( s_{i} - s^{pMC}_{i} ) / s^{pMC}_{i} ]' ) ;
  xlabel( 'Mean [ pMC Market Share ]' ) ;

  legend( 'Monte Carlo', 'Stroud 11-1 Left', 'Stroud 11-1 Right',     ...
        'Product Rule (7^5 Nodes)', 'Sparse Grids' ) ;

  
  hold off ;
  
  %% Save
  
  saveas( gcf, sprintf( '%s/Fig-RelSharesErrorN%05d.new.jpg', szDataDir, nMCDraws ), 'jpg' ) ;
