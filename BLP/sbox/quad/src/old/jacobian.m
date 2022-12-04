function Ddelta = jacobian(theta2)

% JACOBIAN
% Computes the Jacobian of the mean utilities associated with
% the random coefficients Logit demand model.
% I.e., compute D delta_jt / D theta_2
%
% source: Dube, Fox and Su (2009)
% Code Revised: March 2009

global prods T nn v K x delta
global Q_WEIGHTS

expmu      = exp( x * diag( theta2 ) * v ) ;      % exponentiated deviations from mean utilities
expmeanval = exp( delta ) ;
[EstShare, simShare] = ind_shnormMPEC(expmeanval,expmu);

ooo          = ones(1,K+1);
ooo1         = ones(prods,1);
dSdtheta2    = zeros(T*prods,K+1);
dSddeltaDIAG = zeros(T*prods,prods);

for tt=1:T
  index = (0:prods-1)'*T+tt;

  % Process each quadrature node
  for rr = 1:nn
% %         dSddeltaDIAG(index,:) = dSddeltaDIAG(index,:) + (diag(simShare(index,rr)) - simShare(index,rr)*simShare(index,rr)')/nn;
% %         dSdtheta2(index,:) = dSdtheta2(index,:) + (simShare(index,rr)*ooo).*(ooo1*v(:,rr)').*( x(index,:) - (ooo1*(simShare(index,rr)'*x(index,:))))/nn;

    dSddeltaDIAG(index,:) = dSddeltaDIAG(index,:) + (diag(simShare(index,rr))                               ...
                          - simShare(index,rr)*simShare(index,rr)') * Q_WEIGHTS( rr ) ;
    dSdtheta2(index,:) = dSdtheta2(index,:) + (simShare(index,rr)*ooo).*(ooo1*v(:,rr)').*( x(index,:)       ...
                       - (ooo1*(simShare(index,rr)'*x(index,:)))) * Q_WEIGHTS( rr ) ;        
  end
end

%% Do I need to divide by EstShare?  It is not in JP's implementation...
%
% Apparently not because the EstShare terms cancel out on inversion...

% % dSdtheta2    = dSdtheta2    ./ ( EstShare * ooo ) ;
% % dSddeltaDIAG = dSddeltaDIAG ./ ( EstShare * ooo1' ) ;


%% Back to original code...

dSddelta = zeros(T*prods, T*prods);

for i = 1:prods
  for j = 1:prods
    dSddelta( (j-1)*T+1:j*T, (i-1)*T+1:i*T ) = diag(dSddeltaDIAG((j-1)*T+1:j*T, i)); 
  end
end

% This is D delta_jt / D theta_2
Ddelta = -inv(dSddelta)*dSdtheta2;
