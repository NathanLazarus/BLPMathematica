function [theta1,theta2,G,INFO,CPUt,FuncEval,GradEval,MPECx] = runGMMMPECKnitro(x0,gradient)

% source: Dube, Fox and Su (2009)
% Code Revised: March 2009
%
% modification history
% --------------------
% 13dec2020 bss add support for KNITRO 12.3.0
% 21aug2010 bss added box constraints to avoid bad regions
%

global iteration K T prods IV share
global prodsMarket marketStarts marketEnds numProdsTotal

disp('***Solving for local optimum using KNITRO 12.3.0');

if nargin<2,
    % gradient=0;
    error('Error: only gradient==1 is supported');
end

nBeta = K + 1 ;      % number of parameters
nIV = size(IV, 2);

iteration = 0;
INFO = -1;
reps =0;
CPUt = 0;
differ = 1;
while INFO == -1 && reps <= 1
    
    %% Setup problem specification for KNITRO

    % Setup bounds for optimum
    x_L = -Inf*ones(length(x0),1);   % Lower bounds for x.
    x_L(K+2:2*K+2) = 0;              % standard deviations are nonnegative
    x_U =  Inf*ones(length(x0),1);   % Upper bounds for x.
    

    nx0 = size(x0,1);   % number of MPEC optimization parameters
    
    % Sparsity patterns of constraints, 0 if derivative 0, 1 otherwise
    % Derivatives of market shares
    
    c11 = zeros(numProdsTotal, K+1);   % market shares with respect to mean parameters
    c12 = ones(numProdsTotal, K+1);  % market shares with respect to standard deviations
    
    c13 = zeros(numProdsTotal,numProdsTotal); % market shares with respect to mean utilities
    for t=1:T
        c13(marketStarts(t):marketEnds(t),marketStarts(t):marketEnds(t)) = 1;
    end
    
    c14 = zeros(numProdsTotal, nIV); % market shares with respect to moment values
    
    % Derivatives of moments
    
    c21 = ones(nIV, K+1);   % moments with respect to mean parameters 
    c22 = zeros(nIV, K+1);    % moments with respect to standard deviations
    c23 = ones(nIV, numProdsTotal);  % moments with respect to mean utilities
    c24 = eye(nIV);       % moments with respect to moment values 
    
    ConsPattern = [ c11 c12 c13 c14; c21 c22 c23 c24 ];   
        
    % Hessian pattern
    
    HessianPattern = zeros(nx0, nx0);
    HessianPattern(2*K+3+numProdsTotal:end, 2*K+3+numProdsTotal:end) = ones(nIV, nIV);
    HessianPattern(K+2:2*K+2,K+2:2*K+2+numProdsTotal) = ones(K+1, K+1+numProdsTotal);
    HessianPattern(2*K+3:2*K+2+numProdsTotal,K+2:2*K+2) = ones(numProdsTotal, K+1);
    for t=1:T
        range = marketStarts(t):marketEnds(t);
        index = 2*K+2+range;
        HessianPattern(index, index) = ones(prodsMarket(t), prodsMarket(t));
    end 


    %% Find optimum using KNITRO 12.3.0
    ktropts = knitro_options('algorithm',1,'outlev',4, ...
                     'gradopt',1,'hessopt',1,'maxit',1000, ...
                     'feastol',1e-8,'opttol',1e-8, ...
                     'derivcheck', 0, 'derivcheck_tol', 1E-4);      % XXX Enable dervicheck by setting to 1, 2, or 3 (all)
                     % 'par_numthreads', 12   % this option seems to make KNITRO slower
                    
    extendedFeatures = {};
    extendedFeatures.JacobPattern = ConsPattern;
    extendedFeatures.HessPattern = HessianPattern;
    extendedFeatures.HessFcn = @GMMMPEC_hess;
            
    t1 = cputime;   
       
    [X, fval_rep, exitflag, output, lambda, f_grad, f_hess] = knitro_nlp(@GMMMPEC_f, x0, ...
                                                                         [], [], [], [], x_L, x_U, ...
                                                                         @GMMMPEC_c, extendedFeatures, ...
                                                                         ktropts, 'knitroOptions2.opt') ;
                                                                     
                                                                     
    CPUt = cputime - t1;
    MPECx = X;
    G = fval_rep;
    INFO = exitflag;
    theta1 = MPECx(1:K+1);
    theta2 = MPECx(K+2:2*K+2);
    FuncEval = output.funcCount;
    GradEval = nan;
    
    x0 = MPECx;
    reps = reps+1;
end
