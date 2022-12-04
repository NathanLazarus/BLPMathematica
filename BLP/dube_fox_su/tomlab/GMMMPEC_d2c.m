function [d2c] = GMMMPEC_d2c(x0,lambda) 

% GMMMPEC_c
% Constraints for the random coefficients Logit estimated via MPEC.
%
% source: Dube, Fox and Su (2010)
% Code Revised: February 2012

global x K prods T v x nn W

theta1 = x0(1:K+1, 1);                              % mean tastes
theta2 = x0(K+2:2*K+2, 1);                          % st. deviation of tastes
delta = x0(2*K+3:2*K+2+T*prods, 1);                 % mean utilities
g = x0(2*K+3+T*prods:end, 1);                       % moment condition values

nx0 = size(x0,1);
ng = size(g,1);
ooo = ones(1,K+1);
ooo1 = ones(prods,1);

d2c = zeros(nx0, nx0);

expmu = exp(x*diag(theta2)*v);      % exponentiated deviations from mean utilities
expmeanval = exp(delta);
[EstShare, simShare] = ind_shnormMPEC(expmeanval,expmu);

dL2dtheta22 = zeros(K+1,K+1);
dL2dtheta22rr = zeros(K+1,K+1);
dL2ddeltadtheta = zeros(T*prods, K+1);
dL2ddelta2DIAG = zeros(T*prods,prods);
dL2ddelta2 = zeros(T*prods, T*prods);

    
% Evaluate the hessian    
for tt=1:T,
    index = (1:prods)'+(tt-1)*prods;    
    multip = lambda(index);
    for rr = 1:nn,
        
        simS = simShare(index,rr);
        sumprod_mpsimS = multip'*simS;
        blk1 = sumprod_mpsimS*(-diag(simS) + 2*simS*simS');
        blk2 = -(multip.*simS)*simS';
        blk3 = diag(multip.*simS);   
        blk = blk1 + blk2 + blk2'+blk3;        
        dL2ddelta2DIAG(index,:) = dL2ddelta2DIAG(index,:) + blk/nn;
              
        xsimSx = x(index,:) - (ooo1*(simS'*x(index,:)));
        xsimSxv = xsimSx.*(ooo1*v(:,rr)');
        dSdtheta2rr = (simS*ooo).*xsimSxv;
        dL2ddeltadthetarr = -simS*multip'*dSdtheta2rr - sumprod_mpsimS*dSdtheta2rr + (multip*ooo).*dSdtheta2rr;
        dL2ddeltadtheta(index,:) = dL2ddeltadtheta(index,:) + dL2ddeltadthetarr/nn;     
        dL2dtheta22rr = ((multip*ooo).*dSdtheta2rr)'*xsimSxv + sumprod_mpsimS*(-dSdtheta2rr'*x(index,:).*(ones(K+1,1)*v(:,rr)')) ;
        dL2dtheta22 = dL2dtheta22 + dL2dtheta22rr/nn;    
    end
    dL2ddelta2(index,index) = dL2ddelta2DIAG(index,:);    
end

d2c(K+2:2*K+2,K+2:2*K+2) = dL2dtheta22;    
d2c(2*K+3:2*K+2+T*prods,K+2:2*K+2) = dL2ddeltadtheta;
d2c(K+2:2*K+2,2*K+3:2*K+2+T*prods) = dL2ddeltadtheta';        
d2c(2*K+3:2*K+2+T*prods, 2*K+3:2*K+2+T*prods) = dL2ddelta2;     
d2c = (d2c+d2c')/2;
    
    