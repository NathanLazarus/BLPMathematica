function [share,nopurch] = mksharesim(betatrue,x,xi,rc)

%%%%%%%%%%%%
% MKSHARESIM
% generates market share data.
%
% source: Dube, Fox and Su (2010)
% Code Revised: April 2008

global prods T datasort

delta = x*betatrue+xi;

if nargin==4,
    MU = x*rc;
    numer = exp( repmat(delta,1,size(rc,2)) + MU );
else
    numer = exp(delta);
end
temp = cumsum(numer);
temp1 = temp( (1:T)*prods,: );
temp1(2:end,:) = diff(temp1);
denom = 1+temp1;
share = mean(numer./kron(denom,ones(prods,1)),2);
nopurch = mean(1./denom,2);
