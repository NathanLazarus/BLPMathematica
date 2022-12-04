function [cineq, c, dcineq, dc] = GMMMPEC_c(x0) 


% GMMMPEC_c
% Constraints for the random coefficients Logit estimated via MPEC.
%
% Arguments:
%   cineq:  non-linear inequality constraints
%   c:      non-linear equality constraints
%   dcineq: gradient of cineq
%   dc:     gradient of c
%
% source: Dube, Fox and Su (2009)
% Code Revised: March 2009


global x IV K prods T v x
global share

theta1 = x0(1:K+1, 1);                      % mean tastes
theta2 = x0(K+2:2*K+2, 1);                  % st. deviation of tastes
delta = x0(2*K+3:2*K+2+T*prods, 1);         % mean utilities
g = x0(2*K+3+T*prods:end, 1);               % moment condition values

cong = g - IV'*(delta - x*theta1);          % constraints on moment conditions

for i = 1:K+1
    nu(i,:) = theta2(i) * v(i,:) ;          % Note: this uses the same draws for each (J,T) integral
end
MU = x*nu;
expmu = exp(MU);                            % Exponentiated deviations from the mean utility
expmeanval = exp(delta);

vProb = ind_shnorm( expmeanval, expmu ) ;   % Estimated share probability
logEstShares = log( vProb ); % constraints on predicted shares  
% % logEstShares = log(ind_shnorm_stable(delta,MU)); % constraints on predicted shares  

% Set NaNs to Inf to avoid bad regions
% logEstShares( isnan(logEstShares ) ) = 1e10 ; % Inf
% cong( isnan( cong ) )                = 1e10 ; % Inf

cineq = [];
dcineq = [];

% c = [logEstShares ;
%      cong ];
c = [vProb - share ;
     cong ];

if nargout > 2
    dc = GMMMPEC_dc(x0);
end

% % Validate( c ) ;


