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


global x IV K prods T v x nn ;
global share Q_WEIGHTS numProdsTotal;

% 1) Assemble parameters
theta1 = x0(1:K+1, 1);                      % mean tastes
theta2 = x0(K+2:2*K+2, 1);                  % st. deviation of tastes
delta = x0(2*K+3:2*K+2+numProdsTotal, 1);         % mean utilities
g = x0(2*K+3+numProdsTotal:end, 1);               % moment condition values

cong = g - IV'*(delta - x*theta1);          % constraints on moment conditions

% 2) Compute exponentiated deviations from the mean utilities
for i = 1:K+1
    nu(i,:) = theta2(i) * v(i,:) ;          % Note: this uses the same draws for each (J,T) integral
end

% 3) Compute predicted shares
MU = x*nu;
expmu = exp(MU);                            % Exponentiated deviations from the mean utility
% expmu = exp(x*diag(theta2)*v);      % exponentiated deviations from mean utilities
expmeanval = exp(delta);

[EstShare, simShare] = ind_shnormMPEC(expmeanval,expmu);   % Estimated share probability

% vProb = ind_shnorm( expmeanval, expmu ) ;   % Estimated share probability
% logEstShares = log( vProb ); % constraints on predicted shares  
% % logEstShares = log(ind_shnorm_stable(delta,MU)); % constraints on predicted shares  

% Set NaNs to Inf to avoid bad regions
% logEstShares( isnan(logEstShares ) ) = 1e10 ; % Inf
% cong( isnan( cong ) )                = 1e10 ; % Inf

cineq = [];
dcineq = [];

% c = [logEstShares ;
%      cong ];
c = [EstShare - share ;
     cong ];

if nargout > 2

    % Indices
    ng           = size(g,1);
    ooo          = ones(1,K+1);
    ooo1         = ones(prods,1);
    dSdtheta2    = zeros(numProdsTotal,K+1);
    dSddeltaDIAG = zeros(numProdsTotal,prods);

    % 4) Evaluate the Gradients

    for tt=1:T
        index = (0:prods-1)'*T+tt;
        for rr = 1:nn
    % %         dSddeltaDIAG(index,:) = dSddeltaDIAG(index,:) + (diag(simShare(index,rr)) - simShare(index,rr)*simShare(index,rr)')/nn;
    % %         dSdtheta2(index,:) = dSdtheta2(index,:) + (simShare(index,rr)*ooo).*(ooo1*v(:,rr)').*( x(index,:) - (ooo1*(simShare(index,rr)'*x(index,:))))/nn;
            dSddeltaDIAG(index,:) = dSddeltaDIAG(index,:) + (diag(simShare(index,rr))                               ...
                                  - simShare(index,rr)*simShare(index,rr)') * Q_WEIGHTS( rr ) ;
            dSdtheta2(index,:) = dSdtheta2(index,:) + (simShare(index,rr)*ooo).*(ooo1*v(:,rr)').*(x(index,:)       ...
                               - (ooo1*(simShare(index,rr)'*x(index,:)))) * Q_WEIGHTS( rr ) ;
        end
        
%          dSdtheta2(index,:) = dSdtheta2(index,:) ./ (EstShare(index)*ooo);
%          dSddeltaDIAG(index,:) = dSddeltaDIAG(index,:) ./ (EstShare(index)*ooo1');
    end
    
    % XXX dSdtheta2    = dSdtheta2    ./ ( EstShare * ooo ) ;
    % XXX dSddeltaDIAG = dSddeltaDIAG ./ ( EstShare * ooo1' ) ;
        
    dSddelta = zeros(numProdsTotal, numProdsTotal);
    
    for i = 1:prods
      for j = 1:prods
        dSddelta( (j-1)*T+1:j*T, (i-1)*T+1:i*T ) = diag(dSddeltaDIAG((j-1)*T+1:j*T, i));     
      end
    end
    
    dc1 = [zeros(numProdsTotal,K+1), dSdtheta2, dSddelta, zeros(numProdsTotal,ng)];
    dc2 = [IV'*x, zeros(ng, K+1), -IV', eye(ng)];
    dc = [dc1; dc2];
    
    % dc( isnan( dc ) ) = 1e10 ; % Inf
    % % Validate( dc ) ;

end

% % Validate( c ) ;
