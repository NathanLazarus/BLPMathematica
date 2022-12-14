function Ddelta = jacobian(theta2)

% JACOBIAN
% Computes the Jacobian of the mean utilities associated with
% the random coefficients Logit demand model.
% I.e., compute D delta_jt / D theta_2
%
% source: Dube, Fox and Su (2012)
% Code Revised: January 2012

global prods T nn v K x delta
global prodsMarket numProdsTotal marketStarts marketEnds
global Q_WEIGHTS

expmu = exp(x*diag(theta2)*v);      % exponentiated deviations from mean utilities
expmeanval = exp(delta);
[EstShare, simShare] = ind_shnormMPEC(expmeanval,expmu);

ooo = ones(1,K+1);

dSdtheta2 = zeros(numProdsTotal,K+1);
dSddeltaDIAG = zeros(numProdsTotal,prods);
dSddelta = zeros(numProdsTotal, numProdsTotal);

for t=1:T,
    index = marketStarts(t):marketEnds(t);
    ooo1 = ones(prodsMarket(t),1);
    for rr = 1:nn,
        dSddeltaDIAG(index,1:prodsMarket(t)) = dSddeltaDIAG(index,1:prodsMarket(t)) + (diag(simShare(index,rr)) - simShare(index,rr)*simShare(index,rr)') * Q_WEIGHTS(rr);
        dSdtheta2(index,:) = dSdtheta2(index,:) + (simShare(index,rr)*ooo).*(ooo1*v(:,rr)').*( x(index,:) - (ooo1*(simShare(index,rr)'*x(index,:)))) * Q_WEIGHTS(rr);
    end
    dSddelta(index,index) = dSddeltaDIAG(index,1:prodsMarket(t));
end


% This is D delta_jt / D theta_2
Ddelta = -inv(dSddelta)*dSdtheta2;
