function [theta1,theta2,G,INFO,CPUt,FuncEval,GradEval, d2L, lambda, delta, Iteration] = runGMMMPECTomlab(x00,gradient)

% source: Dube, Fox and Su (2010)
% Code Revised:  February 2012

global theta1 K T prods IV tol_outer share nIV

if nargin<2,
    gradient=0;
end


CPUt = 0;

x_L = -Inf*ones(length(x00),1);   % Lower bounds for x.
x_L(K+2:2*K+2) = 0;
x_U =  Inf*ones(length(x00),1);   % Upper bounds for x.

neF    = prods*T + size(IV, 2);
c_L = zeros(neF,1);
c_U = zeros(neF,1);

c_L(1:prods*T) = share;
c_U(1:prods*T) = share;

nx0 = size(x00,1);

c11 = zeros(T*prods, K+1);
c12 = ones(T*prods, K+1);
c13 = kron(eye(T), ones(prods,prods));
c14 = zeros(T*prods, nIV);

c21 = ones(nIV, K+1);
c22 = zeros(nIV, K+1);
c23 = ones(nIV, T*prods);
c24 = eye(nIV);

ConsPattern = [ c11 c12 c13 c14; c21 c22 c23 c24 ];

ConsPattern = ceil(abs(ConsPattern))./max(ceil(abs(ConsPattern)),1);

HessPattern = zeros(nx0, nx0);
HessPattern(2*K+3+T*prods:end, 2*K+3+T*prods:end) = ones(nIV, nIV);
d2LPatt = zeros(nx0, nx0);
d2LPatt(2*K+3+T*prods:end, 2*K+3+T*prods:end) = ones(nIV, nIV);
d2LPatt(K+2:2*K+2,K+2:2*K+2+T*prods) = ones(K+1, K+1+T*prods);
d2LPatt(2*K+3:2*K+2+T*prods,K+2:2*K+2) = ones(T*prods, K+1);
for tt=1:T,
    index = 2*K+2+(1:prods)'+(tt-1)*prods;
    d2LPatt(index, index) = ones(prods, prods);
end        % specify sparsity pattern for constraints

Name = 'BLP MPEC Problem';
fLowBnd = 0;                                                % Lower bound on function.

if gradient==1,
    Prob = conAssign('GMMMPEC_f', 'GMMMPEC_grad', 'GMMMPEC_hess', HessPattern, x_L, x_U, Name, x00,...
        [], fLowBnd, [], [], [], 'GMMMPEC_c', 'GMMMPEC_dc', 'GMMMPEC_d2c', ConsPattern, c_L, c_U);
else
    Prob = conAssign('GMMMPEC_f', [], 'GMMMPEC_hess', HessPattern, x_L, x_U, Name, x00,...
        [], fLowBnd, [], [], [], 'GMMMPEC_c', [], [], ConsPattern, c_L, c_U);
end
Prob.d2LPattern = d2LPatt;
Prob.PriLevOpt = 3;
Prob.KNITRO.PrintFile = '';
Prob.KNITRO.options.ALG = 1;
Prob.KNITRO.options.MAXIT = 200;    % Setting maximum number of iterations
% Prob.KNITRO.options.BAR_MAXBACKTRACK = 5;
% Prob.KNITRO.options.BAR_MAXREFACTOR = 10;
% Prob.KNITRO.options.HESSOPT = 4;    % to be used along with first-order derivatives

Prob.KNITRO.options.OPTTOL = tol_outer ;    % set outer-loop tolerance

MPEC_Result = tomRun('knitro', Prob, 1);

MPECx = MPEC_Result.x_k;
G = MPEC_Result.f_k;
d2L = MPEC_Result.H_k;
lambda = MPEC_Result.v_k;
INFO = MPEC_Result.Inform;
CPUt = MPEC_Result.CPUtime;
Iteration = MPEC_Result.Iter;
FuncEval = MPEC_Result.FuncEv;
GradEval = MPEC_Result.GradEv;
theta1 = MPECx(1:K+1);
theta2 = MPECx(K+2:2*K+2);
delta = MPECx(2*K+3:2*K+2+T*prods);

x00 = MPECx;
reps = 1;


