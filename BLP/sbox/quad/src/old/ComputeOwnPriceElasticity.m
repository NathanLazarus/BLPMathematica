% ComputeOwnPriceElasticity - computes own-price elasticity
%
% modification history
% --------------------
% 28mar2010 bss written.
%

function [ vOwn ] = ComputeOwnPriceElasticity( x0 )


%% Compute expmu and expmeanval

global x K prods T v numProdsTotal

theta1 = x0(1:K+1, 1);                      % mean tastes
theta2 = x0(K+2:2*K+2, 1);                  % st. deviation of tastes

dwPriceCoef = theta1( end ) ;
vOwnPrice   = x( :, end ) ;

delta = x0(2*K+3:2*K+2+numProdsTotal, 1);         % mean utilities

for i = 1:K+1,
    nu(i,:) = theta2(i)*v(i,:);
end

MU = x*nu;
expmu = exp(MU);                            % Exponentiated deviations from the mean utility
expmeanval = exp(delta);


%% Compute derivative of shares

global oo sharesum denomexpand 

numer = (expmeanval*oo ).*expmu;        % this is the numerator (oo speeds-up expanding mean utility by number of draws)
sum1 = sharesum*numer;
sum11 = 1./(1+sum1);                    % this is the denominator of the shares
denom1 = sum11(denomexpand,:);          % this expands the denominator
mCondShares = numer.*denom1;                     % simulated shares for each draw

% % PROB = mean(SS,2);                      % expected share (i.e. mean across draws)
                                        % XXX replace with quad weights

%% Integrate using Gaussian-Hermite Quadrature

global Q_WEIGHTS ;
vProb       = mCondShares * Q_WEIGHTS ;
vGradShares = dwPriceCoef( 1 ) * ( mCondShares .* ( 1 - mCondShares ) ) * Q_WEIGHTS ;    

vOwn = vOwnPrice .* vGradShares ./ vProb ;

