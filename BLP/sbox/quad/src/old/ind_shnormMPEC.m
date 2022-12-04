function [EstShare, simShare] = ind_shnormMPEC(expmeanval,expmu)

% IND_SHNORM
% This function computes the distribution of "individual" probabilities of choosing each brand.
% The probabilities are those associated with the normally-distributed r.c. logit.
%
% ARGUMENTS:
% expmeanval = vector of exponentiated mean utilities, sorted by product then market
% expmu = matrix of exponentiated draws of deviations from mean utilities, sorted by product then market
%
% OUTPUT:
% PROB = vector of expected market shares, sorted by product then market
%
% source: Dube, Fox and Su (2008)
% Code Revised: April 2008

global oo sharesum denomexpand

numer = (expmeanval*oo ).*expmu;        % this is the numerator (oo speeds-up expanding mean utility by number of draws)
sum1 = sharesum*numer;
sum11 = 1./(1+sum1);                    % this is the denominator of the shares
denom1 = sum11(denomexpand,:);          % this expands the denominator
simShare = numer.*denom1;               % simulated shares for each draw
                                        % XXX Do I need a quadrature weight
                                        % here?

% % EstShare = mean(simShare,2);            % expected share (i.e. mean across draws)

%% Use Gaussian-Hermite Quadrature

global Q_WEIGHTS ;
EstShare = simShare * Q_WEIGHTS ;

