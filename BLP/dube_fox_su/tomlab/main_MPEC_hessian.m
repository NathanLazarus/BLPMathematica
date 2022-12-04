%%%%%%%%%%%%%%%%%%%%%%%
% MAIN SCRIPT
%
% 1) Simulate Data from Random Coefficients Logit Demand model
% 2) Estimate model using MPEC
%
% source: Dube, Fox and Su (2010)
% Code Revised:  February 2012
%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%
% SETTINGS
%%%%%%%%%%%
global share W PX x IV rc expmeanval nIV % logshare
global datasort denomexpand sharesum oo
global T prods nn K v rc rctrue oo oo1 tol_inner tol_outer delta

% aa=maxNumCompThreads(1);                            % disable multi-thread processing b/c Knitro cannot exploit this so it slows the code down

diary output_MPEC_Hessian_Tomlab.out;

randn('seed',155)                                   % reset normal random number generator
rand('seed',155)                                    % reset uniform random number generator
nn = 500;                                           % # draws to simulate shares
tol_inner = 1.e-14;
tol_outer = 1.e-6;
prods = 25;                                         % # products
T = 50;                                             % # Markets
starts = 5;                                         % # random start values to check during estimation

%%%%%%%%%%%%%%%%%%%%%%%%
% TRUE PARAMETER VALUES
%%%%%%%%%%%%%%%%%%%%%%%%
costparams = ones(6,1)/2;
betatrue = [2 1.5 1.5 .5 -3]';                     % true mean tastes
betatrue0 = betatrue;
K = size(betatrue,1)-1;                             % # attributes
covrc = diag( [ .5 .5 .5 .5 .2] );                  % true co-variances in tastes
rctrue = covrc(find(covrc));
thetatrue = [betatrue;rctrue];
v = randn(length(betatrue),nn);                     % draws for share integrals during estimation
rc = chol(covrc)'*v;                                % draws for share integrals for data creation
sigmaxi=1;
covX = -.8*ones(K-1)+1.8*eye(K-1);                  % allow for correlation in X's
covX(1,3) = .3;
covX(2,3) = .3;
covX(3,1) = .3;
covX(3,2) = .3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create indices to speed share calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oo = ones(1,nn);                                    % index for use in simulation of shares
oo1 = (0:T-1)*prods+prods;                          % indicates last observation for each market
sharesum = kron(speye(T), ones(1,prods));           % used to create denominators in predicted shares (i.e. sums numerators)
datasort = [(1:T*prods)' repmat((1:T)',prods,1) kron( (1:prods)',ones(T,1))];
datasort = datasort(:,1);
denomexpand = kron((1:T)',ones(prods,1));

%%%%%%%%%%%%%%%%%%%%%
% SIMULATE DATA
%%%%%%%%%%%%%%%%%%%%%
randn('seed',5000)                                   % reset normal random number generator
rand('seed',5000)                                    % reset uniform random number generator
xi = randn( T*prods,1)*sigmaxi;                      % draw demand shocks
A = [kron(ones(T,1), randn( prods,K-1)*chol(covX))];  % product attributes
prand = rand( T*prods,1)*5;
price = 3 +   xi*1.5 +  prand + sum(A,2);
z = rand( T*prods,length(costparams)) + 1/4*repmat( abs( prand +  sum(A,2)*1.1 ) ,1,length(costparams));
x = [ones(T*prods,1) A price];
[share,nopurch] = mksharesim(betatrue,x,xi,rc);
y = log(share) - log(kron(nopurch,ones(prods,1)));      % log-odds ratio
iv = [ones(prods*T,1) A z];                         % matrix of instruments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2SLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xhat = iv*inv(iv'*iv)*iv'*x;
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
IterMPEC = 0;
expmeanval = expmeanval0;

for reps=1:starts,
    
    %%% NOTE: START VALUES SET TO GIVE SAME INITIAL GMM OBJECTIVE FUNCTIONAL
    %%% VALUE AS WITH NFP
    theta20 = startvalues(reps,:)';
    delta0 = invertshares(theta20);
    theta10 = PX*delta0;
    resid0 = delta0 - x*theta10;
    g0 = IV'*resid0;
    x0 = [theta10; theta20; delta0; g0];
    
    [theta1MPEC_rep,theta2MPEC_rep,GMPEC_rep,INFOMPEC_rep,CPUtMPEC_rep,FuncEvalMPEC_rep,GradEvalMPEC_rep, d2L, lambda, deltaMPEC, IterMPEC_rep] = runGMMMPECTomlab(x0,1);
    
    CPUtMPEC = CPUtMPEC + CPUtMPEC_rep;
    FuncEvalMPEC = FuncEvalMPEC + FuncEvalMPEC_rep;
    GradEvalMPEC = GradEvalMPEC + GradEvalMPEC_rep;
    IterMPEC = IterMPEC + IterMPEC_rep;
    
    if (GMPEC_rep < GMPEC && INFOMPEC_rep==0),
        theta1MPEC = theta1MPEC_rep;
        theta2MPEC = theta2MPEC_rep;
        GMPEC = GMPEC_rep;
        INFOMPEC = INFOMPEC_rep;
    end
    
    
end

% COVARIANCE MATRIX FOR MPEC STRUCTURAL PARAMETERS
deltaSE = invertshares(thetaMPEC2);    % mean utilities
resid = deltaSE - x*theta1MPEC;  % xi
Ddelta = jacobian(theta2MPEC);       % jacobian matrix of mean utilities
covg = zeros(size(IV,2));
for ii =1:length(IV),
    covg = covg + IV(ii,:)'*IV(ii,:)*(resid(ii)^2);
end
Dg = [x Ddelta]'*IV;            % gradients of moment conditions
covMPEC = inv( Dg*W*Dg')*Dg*W*covg*W*Dg'*inv( Dg*W*Dg');
stderr = sqrt(diag(covMPEC));

save(['MPEC_Hess_' num2str(T) 'mkts_' num2str(prods) 'products',num2str(nn),'draws']);

diary off;