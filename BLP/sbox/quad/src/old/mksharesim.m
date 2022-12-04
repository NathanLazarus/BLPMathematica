function [share,nopurch] = mksharesim(betatrue,x,xi,rc)

%%%%%%%%%%%%
% MKSHARESIM
% generates market share data.
%
% source: Dube, Fox and Su (2008)
% Code Revised: April 2008

global prods T datasort

delta = x*betatrue+xi;

if nargin==4,
    MU = x*rc;
    numer = exp( repmat(delta,1,size(rc,2)) + MU );
else
    numer = exp(delta);
end
temp = cumsum(numer(datasort,:));
temp1 = temp( (1:T)*prods,: );
temp1(2:end,:) = diff(temp1);

denom = 1+temp1;  % This is sum( over j of exp( V_jt ) )  : ( T x nn )


% XXX multiply by quad
% weights numer./denom * Q_WEIGHTS share is ordered by data sort, so each
% product in a market is T slots apart

% % share = mean(numer./repmat(denom,prods,1),2);   
% %
% % nopurch = mean(1./denom,2);

%% Calculate Using Gaussian-Hermite Quadrature

global Q_WEIGHTS ;
share = ( numer ./ repmat( denom, prods, 1 ) ) * Q_WEIGHTS ;

nopurch = ( 1 ./ denom ) * Q_WEIGHTS ;
