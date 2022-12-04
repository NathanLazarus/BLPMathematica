function dc = GMMMPEC_dc(x0)

% GMMMPEC_dc
% Gradient of the constraints of the random coefficients Logit estimated
% via MPEC.
%
% source: Dube, Fox and Su (2009)
% Code Revised: March 2009

global x IV K prods v nn T

% 1) Assemble parameters
theta2 = x0(K+2:2*K+2, 1);                  % st. deviation of tastes
delta = x0(2*K+3:2*K+2+T*prods, 1);         % mean utilities
g = x0(2*K+3+T*prods:end, 1);               % moment condition values

% 2) Compute exponentiated deviations from the mean utilities
expmu = exp(x*diag(theta2)*v);      % exponentiated deviations from mean utilities

% 3) Compute predicted shares
expmeanval = exp(delta);
[EstShare, simShare] = ind_shnormMPEC(expmeanval,expmu);

% Indices
ng           = size(g,1);
ooo          = ones(1,K+1);
ooo1         = ones(prods,1);
dSdtheta2    = zeros(T*prods,K+1);
dSddeltaDIAG = zeros(T*prods,prods);

% 4) Evaluate the Gradients

global Q_WEIGHTS ;

for tt=1:T,
    index = (0:prods-1)'*T+tt;
    for rr = 1:nn
% %         dSddeltaDIAG(index,:) = dSddeltaDIAG(index,:) + (diag(simShare(index,rr)) - simShare(index,rr)*simShare(index,rr)')/nn;
% %         dSdtheta2(index,:) = dSdtheta2(index,:) + (simShare(index,rr)*ooo).*(ooo1*v(:,rr)').*( x(index,:) - (ooo1*(simShare(index,rr)'*x(index,:))))/nn;
        dSddeltaDIAG(index,:) = dSddeltaDIAG(index,:) + (diag(simShare(index,rr))                               ...
                              - simShare(index,rr)*simShare(index,rr)') * Q_WEIGHTS( rr ) ;
        dSdtheta2(index,:) = dSdtheta2(index,:) + (simShare(index,rr)*ooo).*(ooo1*v(:,rr)').*( x(index,:)       ...
                           - (ooo1*(simShare(index,rr)'*x(index,:)))) * Q_WEIGHTS( rr ) ;
    end
    
% %     dSdtheta2(index,:) = dSdtheta2(index,:) ./ (EstShare(index)*ooo);
% %     dSddeltaDIAG(index,:) = dSddeltaDIAG(index,:) ./ (EstShare(index)*ooo1');
end

dSdtheta2    = dSdtheta2    ./ ( EstShare * ooo ) ;
dSddeltaDIAG = dSddeltaDIAG ./ ( EstShare * ooo1' ) ;
    
dSddelta = zeros(T*prods, T*prods);

for i = 1:prods
  for j = 1:prods
    dSddelta( (j-1)*T+1:j*T, (i-1)*T+1:i*T ) = diag(dSddeltaDIAG((j-1)*T+1:j*T, i));     
  end
end

dc1 = [zeros(T*prods,K+1), dSdtheta2, dSddelta, zeros(T*prods,ng)];
dc2 = [IV'*x, zeros(ng, K+1), -IV', eye(ng)];
dc = [dc1; dc2];

dc( isnan( dc ) ) = 1e10 ; % Inf
% % Validate( dc ) ;

