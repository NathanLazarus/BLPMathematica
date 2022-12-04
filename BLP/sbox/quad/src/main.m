%%  XXX THIS FILE IS DEPRECATED AND NOT YET UPDATED


%%%%%%%%%%%%%%%%%%%%%%%
% MAIN SCRIPT
%
% 1) Simulate Data from Random Coefficients Logit Demand model
% 2) Estimate model using MPEC
%
% source: Dube, Fox and Su (2009)
% Code Revised: May 2009
%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
%
% modification history
% --------------------
% 20jan2010 bss added QMC integration.
% 24oct2009 bss changed diary name to current day's date.
% 23oct2009 bss changed so sample data is simulated.
% 22oct2009 bss added Gaussian-Hermite Quadrature.
%

%%
%%%%%%%%%%% 
% SETTINGS
%%%%%%%%%%%

error('main.m is deprecated. Do not use until updated for latest DFS code')

%------------------
% Choose type of quadrature
bUseGHQuadrature = 0 ;
bUseQuasiMC = 0 ;


%------------------
% Setup variables

global share W PX x IV rc expmeanval
global datasort denomexpand sharesum oo numProdsTotal
global T prods nn K v rc rctrue oo oo1 tol_inner tol_outer

N_DRAWS_GEN_DATA = 100 ;        % How many draws to generate the data
N_DRAWS_MC_INT   = 5000 ;      % How many draws to simulate the integral

szDate      = date ;
szDiaryName = [ 'Log' szDate(8:end) szDate(4:6) szDate(1:2) '.out' ] ;

diary( szDiaryName ) ;

% randn('seed',155)                                   % reset normal random number generator
% rand('seed',155)                                    % reset uniform random number generator
randn('seed',231)                                   % reset normal random number generator
rand('seed',159)                                    % reset uniform random number generator

nn = 100 ;                                           % # draws to simulate shares
tol_inner = 1.e-14;
tol_outer = 1.e-6;

% prods = 25;                                         % # products
% T = 50;                                             % # Markets
prods = 15;                                         % # products
T     = 20;                                         % # Markets

starts = 5;                                         % # random start values to check during estimation

%%
%%%%%%%%%%%%%%%%%%%%%%%%
% TRUE PARAMETER VALUES
%%%%%%%%%%%%%%%%%%%%%%%%
costparams = ones(6,1)/2;
betatrue = [2 1.5 1.5 .5 -3]';                     % true mean tastes
betatrue0 = betatrue;
K = size(betatrue,1)-1;                             % # attributes = length( beta ) - 1 (for constant)
covrc = diag( [ .5 .5 .5 .5 .2] );                  % true co-variances in tastes
rctrue = covrc(find(covrc));
thetatrue = [betatrue;sqrt(rctrue)];
v = randn(length(betatrue),nn);                     % draws for share integrals during estimation
rc = chol(covrc)'*v;                                % draws for share integrals for data creation
sigmaxi=1;
covX = -.8*ones(K-1)+1.8*eye(K-1);                  % allow for correlation in X's
covX(1,3) = .3;
covX(2,3) = .3;
covX(3,1) = .3;
covX(3,2) = .3;

% % %% Create Quadrature Nodes and Weights
% % 
% % bUseGHQuadrature = 1 ;
% % 
% % if bUseGHQuadrature 
% %   % Setup Gaussian Quadrature nodes and weights
% %   
% %   GHQuadInit( length(betatrue), 5 ) ;  % Change the last argument to change the number of nodes
% % 
% %   global Q_NODES ;
% %   v  = Q_NODES ;
% %   rc = chol(covrc)' * v ;
% %   nn = size( v, 2 ) ;
% % 
% % else
% %   % Default to basic simulation
% %   
  global Q_WEIGHTS ;
  Q_WEIGHTS = ones( nn, 1 ) / nn ;
% %   
% % end

%% Back to Dube's code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create indices to speed share calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oo = ones(1,nn);                                    % index for use in simulation of shares
oo1 = (0:T-1)*prods+prods;                          % indicates last observation for each market
sharesum = kron( ones(1,prods),speye(T));           % used to create denominators in predicted shares (i.e. sums numerators)

% Create a set of indexes which select by product index and then market
% index.  Necessary because data is stored by market and then product
%
% First create the one-dimensional index in column 1, then market index,
% and finally product index within a market.  Next, extract just the one
% dimensional index in the first column.
datasort = sortrows( [(1:numProdsTotal)' kron( (1:prods)',ones(T,1)) repmat((1:T)',prods,1)] , [3 2]); 
datasort = datasort(:,1);

denomexpand = repmat( (1:T)',prods,1);

%%
%%%%%%%%%%%%%%%%%%%%%
% SIMULATE DATA
%%%%%%%%%%%%%%%%%%%%%
% randn('seed',5000)                                   % reset normal random number generator
% rand('seed',5000)                                    % reset uniform random number generator

randn('seed',3837)                                   % reset normal random number generator
rand('seed',2234)                                    % reset uniform random number generator

xi = randn( numProdsTotal,1)*sigmaxi;                     % draw demand shocks
A = [kron(randn( prods,K-1)*chol(covX),ones(T,1))]; %product attributes
prand = rand( numProdsTotal,1)*5;
price = 3 +   xi*1.5 +  prand + sum(A,2);
z = rand( numProdsTotal,length(costparams)) + 1/4*repmat( abs( prand +  sum(A,2)*1.1 ) ,1,length(costparams));
x = [ones(numProdsTotal,1) A price];
[share,nopurch] = mksharesim(betatrue,x,xi,rc);
y = log(share) - log(repmat(nopurch,prods,1));      % log-odds ratio
iv = [ones(numProdsTotal,1) A z];                         % matrix of instruments


%% Create Quadrature Nodes and Weights

if bUseGHQuadrature 
  % Setup Gaussian Quadrature nodes and weights
  disp( '***Using GH Quadrature' ) ;
  
  GHQuadInit( length(betatrue), 5 ) ;  % Change the last argument to change the number of nodes

  global Q_NODES ;
  v  = Q_NODES ;
%   rc = chol(covrc)' * v ;
  nn = size( v, 2 ) ;
  oo = ones(1,nn);                    % Redefine index for use in simulation of shares
  
elseif bUseQuasiMC
  disp( '***Using qMC' ) ;
  
  nPoints = 1000 ;
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
  nn = N_DRAWS_MC_INT ;
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
beta2sls = PX2sls*y;
se2sls = sqrt(diag( mean((y-x*beta2sls).^2)*inv(xhat'*xhat) ));

expmeanval0 = exp(y);
IV = [ones(numProdsTotal,1) A z A.^2 A.^3 z.^2 z.^3 prod(A,2) prod(z,2) kron(A(:,1),ones(1,size(z,2))).*z  kron(A(:,2),ones(1,size(z,2))).*z];
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

for reps=1:starts,
    %%% NOTE: START VALUES SET TO GIVE SAME INITIAL GMM OBJECTIVE FUNCTIONAL
    %%% VALUE AS WITH NFP
    theta20 = startvalues(reps,:)';
    delta0 = invertshares(theta20);
    theta10 = PX*delta0;
    resid0 = delta0 - x*theta10;
    g0 = IV'*resid0;
    x0 = [theta10; theta20; delta0; g0];
    [theta1MPEC_rep,theta2MPEC_rep,GMPEC_rep,INFOMPEC_rep,CPUtMPEC_rep,FuncEvalMPEC_rep,GradEvalMPEC_rep,xOpt] = runGMMMPECKnitro(x0,1);
    
    CPUtMPEC = CPUtMPEC + CPUtMPEC_rep;
    FuncEvalMPEC = FuncEvalMPEC + FuncEvalMPEC_rep;
    GradEvalMPEC = GradEvalMPEC + GradEvalMPEC_rep;
    
    mTheta1( :, reps ) = theta1MPEC_rep ;
    mTheta2( :, reps ) = theta2MPEC_rep ;
    mXOpt( :, reps )   = xOpt ;
    
    if (GMPEC_rep < GMPEC && INFOMPEC_rep==0),
        thetaMPEC1 = theta1MPEC_rep;
        thetaMPEC2 = theta2MPEC_rep;
        GMPEC = GMPEC_rep;
        INFOMPEC = INFOMPEC_rep;
    end
end


% COVARIANCE MATRIX FOR MPEC STRUCTURAL PARAMETERS
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

save( './Data/SeedB/mpecVarsN5000.mat', 'share', 'W', 'PX', 'x', 'IV', 'rc', 'expmeanval',               ...
      'datasort', 'denomexpand', 'sharesum', 'oo', 'T', 'prods', 'nn', 'K',                 ...
      'v', 'rctrue', 'oo1', 'tol_inner', 'tol_outer', 'mXOpt', 'covMPEC',                   ...
      'betatrue', 'thetatrue', 'covrc' ) ;
    
%% Compare simulation and quadrature

mSharesSim = zeros( prods * T, starts ) ;
mSharesGH  = zeros( prods * T, starts ) ;

for ix = 1 : starts
  
  % Computer shares via Monte Carlo integration
  mSharesSim( :, ix ) = ComputeShares( mXOpt( :, ix ) ) ;   % Market shares, sorted by product, then market

end


% Compute Shares via Gaussian-Hermite Quadrature
GHQuadInit( length(betatrue), 5 ) ;  % Change the last argument to change the number of nodes

if ~bUseGHQuadrature
  global Q_NODES ;
end

v  = Q_NODES ;
%   rc = chol(covrc)' * v ;
nn = size( v, 2 ) ;
oo = ones(1,nn);  


for jx = 1 : starts
  
  mSharesGH( :, jx ) = ComputeShares( mXOpt( :, jx ) ) ;

  for ix = 1 : T
    vIx = ix + ( 0 : 50 : ( T - 1 ) * prods ) ;
    totGH  = sum( mSharesGH( vIx, jx ) ) ;
    totSim = sum( mSharesSim( vIx, jx ) ) ;
    fprintf( 1, 'Market = %4d : (GH SIM diff): %f %f %15.7e\n', ix, totGH, totSim, ( totGH - totSim ) / totGH ) ;
  end

end

diary off;
