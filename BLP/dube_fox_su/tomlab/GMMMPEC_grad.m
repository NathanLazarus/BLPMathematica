function grad_f = GMMMPEC_grad(x0, Prob)

% GMMMPEC_grad
% Gradient of the GMM objective function for the random coefficients Logit
% estimated via MPEC.
%
% source: Dube, Fox and Su (2010)
% Code Revised:  February 2012

global W K prods T

theta1 = x0(1:K+1, 1);                      % mean tastes
theta2 = x0(K+2:2*K+2, 1);                  % st. deviation of tastes
delta = x0(2*K+3:2*K+2+T*prods, 1);         % mean utilities
g = x0(2*K+3+T*prods:end, 1);               % moment condition values
nx0 = size(x0,1);
grad_f = zeros(nx0,1);
grad_f(2*K+3+T*prods:end,1) = 2*W*g;
