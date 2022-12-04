function Ddelta = jacobian(theta2)

% JACOBIAN
% Computes the Jacobian of the mean utilities associated with
% the random coefficients Logit demand model.
%
% source: Dube, Fox and Su (2010)
% Code Revised: March 2009

global prods T nn v K x delta

expmu = exp(x*diag(theta2)*v);      % exponentiated deviations from mean utilities
expmeanval = exp(delta);
[EstShare, simShare] = ind_shnormMPEC(expmeanval,expmu);

ooo = ones(1,K+1);
ooo1 = ones(prods,1);
dSdtheta2 = zeros(T*prods,K+1);
dSddeltaDIAG = zeros(T*prods,prods);
dSddelta = zeros(T*prods, T*prods);

for tt=1:T,
    index = (1:prods)'+(tt-1)*prods;
    for rr = 1:nn,
        dSddeltaDIAG(index,:) = dSddeltaDIAG(index,:) + (diag(simShare(index,rr)) - simShare(index,rr)*simShare(index,rr)')/nn;
        dSdtheta2(index,:) = dSdtheta2(index,:) + (simShare(index,rr)*ooo).*(ooo1*v(:,rr)').*( x(index,:) - (ooo1*(simShare(index,rr)'*x(index,:))))/nn;
    end
    dSdtheta2(index,:) = dSdtheta2(index,:);
    dSddelta(index,index) = dSddeltaDIAG(index,:);
end


Ddelta = -inv(dSddelta)*dSdtheta2;



