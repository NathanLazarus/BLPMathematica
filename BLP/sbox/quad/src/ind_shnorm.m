function PROB = ind_shnorm(expmeanval,expmu)

% IND_SHNORM
% This function computes the "individual" probabilities of choosing each brand.
% The probabilities are those associated with the normally-distributed r.c. logit.
%
% ARGUMENTS:
% expmeanval = vector of exponentiated mean utilities, sorted by product then market
% expmu = matrix of exponentiated draws of deviations from mean utilities, sorted by product then market
%
% OUTPUT:
% PROB = vector of expected market shares, sorted by product then market
%
% source: Dube, Fox and Su (2009)
% Code Revised: March 2009

global oo sharesum marketForProducts SS PROB
global Q_WEIGHTS ;


numer = (expmeanval*oo ).*expmu;        % this is the numerator (oo speeds-up expanding mean utility by number of draws)
sum1 = sharesum*numer;
sum11 = 1./(1+sum1);                    % this is the denominator of the shares
denom1 = sum11(marketForProducts,:);    % this expands the denominator

% XXX This is numerically unstable and needs to be fixed!
SS = numer.*denom1;                     % simulated shares for each draw

% % PROB = mean(SS,2);                      % expected share (i.e. mean across draws)
                                        % XXX replace with quad weights

%% Integrate using Quadrature
PROB = SS * Q_WEIGHTS ;                                        
