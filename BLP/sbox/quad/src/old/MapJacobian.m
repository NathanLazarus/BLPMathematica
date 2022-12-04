function [ mJac ] = MapJacobian( x0 )

% MapJacobian
% 
% Jacobian of the BLP inner loop mapping to verify the conditions for the
% proof that it is a contraction mapping in BLP (1995)
%
% source: Dube, Fox and Su (2009)
% Code Revised: March 2009
% 
% modification history
% --------------------
% 04jun2010 bss hacked from the bones of GMMMPEC_dc.m.
%

global x K prods v nn T share

% 1) Assemble parameters
theta1 = x0( 1 : K + 1, 1 ) ;                  % mean tastes
theta2 = x0( K + 2 : 2 * K + 2, 1 ) ;          % st. deviation of tastes
delta = x0( 2*K+3 : 2*K+2+T*prods, 1 ) ;       % mean utilities
g = x0( 2*K+3+T*prods:end, 1 ) ;               % moment condition values

% 2) Compute exponentiated deviations from the mean utilities
expmu = exp( x * diag( theta2 ) * v ) ;      % exponentiated deviations from mean utilities

% 3) Compute predicted shares
expmeanval = exp( delta ) ;

% EstShare is market share for good j in market t
% simShare is market share for good j in market t conditional on a
%   node/draw
[ EstShare, simShare ] = ind_shnormMPEC( expmeanval, expmu ) ;    

% Indices
% % nx0  = size( x0, 1 ) ;
% % ng   = size( g, 1 ) ;
ooo  = ones( 1, K+1 ) ;
ooo1 = ones( prods, 1 ) ;

dSddeltaDIAG = zeros( T * prods, prods ) ;
dSddelta     = zeros( T * prods, T * prods ) ;

mJac         = zeros( T * prods, prods ) ;

% 4) Evaluate the Gradients

global Q_WEIGHTS ;

vTmp = expmeanval .* share ./ ( EstShare .^ 2 ) ;   % exp(delta_jt) S_jt / s_jt ^ 2
mTmpSimShare = simShare ./ repmat( expmeanval, 1, nn ) ;     % exp(mu_jt) / ( Sum exp( delta_kt + mu_kt ) for a draw/node

for  tt = 1 : T                             % tt is which market to process
    index = ( 0 : prods-1 )' * T + tt ;     % index grabs all products for market tt
    for rr = 1 : nn                         % processes each node/draw
% %         dSddeltaDIAG(index,:) = dSddeltaDIAG(index,:) + (diag(simShare(index,rr))                               ...
% %                               - simShare(index,rr)*simShare(index,rr)') * Q_WEIGHTS( rr ) ;
        mJac( index, : ) = mJac( index, : ) + simShare(index,rr)*mTmpSimShare(index,rr)' * Q_WEIGHTS( rr ) ;                            
    end
% %     dSddeltaDIAG(index,:) = dSddeltaDIAG(index,:) ./ (EstShare(index)*ooo1');
end

% % dSddelta = zeros(T*prods, T*prods);
% % 
% % for i = 1:prods
% %     for j = 1:prods
% %     dSddelta( (j-1)*T+1:j*T, (i-1)*T+1:i*T ) = diag(dSddeltaDIAG((j-1)*T+1:j*T, i)); 
% %     end
% % end

% % dc1 = [zeros(T*prods,K+1), dSdtheta2, dSddelta, zeros(T*prods,ng)];

% % mJac = [dc1; dc2];
