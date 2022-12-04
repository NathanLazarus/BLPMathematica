function PROB = ind_shnorm_stable( vMeanUtil, mMu )

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
%
% modification history
% --------------------
% 20aug2010 bss hacked from ind_shnorm.m to improve numerical stability.
%

global oo sharesum denomexpand SS PROB

% NOTE:  ordering for expmeanval & expmu is 
% ( 
%   (product 1 for markets 1 .. T )
%   (product 2 for markets 1 .. T )
%   ...
%   (product J for markets 1 .. T ) 
% )

%{
expmeanval = exp( vMeanUtil ) ;
expmu      = exp( mMu ) ;

numer = (expmeanval*oo ).*expmu;        % this is the numerator (oo speeds-up expanding mean utility 
                                        % by number of draws -- alternative to repmat)
sum1 = sharesum*numer;
sum11 = 1./(1+sum1);                    % this is the denominator of the shares
denom1 = sum11(denomexpand,:);          % this duplicates the denominator with correct market ordering

% XXX This is numerically unstable and needs to be fixed!
SS = numer.*denom1;                     % simulated shares for each draw

% % PROB = mean(SS,2);                      % expected share (i.e. mean across draws)
                                        % XXX replace with quad weights

%% Integrate using Gaussian-Hermite Quadrature

global Q_WEIGHTS ;
PROB2 = SS * Q_WEIGHTS ;                                        

%}


%% Slow but numerically stable -- something is wrong here

%{

global Q_WEIGHTS T prods

nNodes = length( Q_WEIGHTS ) ;
PROB1   = zeros( length( vMeanUtil ), 1 ) ;

% process each quadrature node
for ixNode = 1 : nNodes
  
  % process each market
  for ixMarket = 1 : T
    vCurRange = ( ixMarket : T : T * prods )' ;  % because data is organized by ( prod 1 in all markets, then prod 2 in 1 .. T, etc.)
    
    vTmp = vMeanUtil( vCurRange ) + mMu( vCurRange, ixNode ) ;  % V_j for all prods in a market
    mTmp = repmat( vTmp, 1, prods ) ;
    
    mTmp = mTmp' - mTmp ;   % mTmp_ij = V_i - V_j, i.e. j is chosen alternative 
    mTmp = exp( mTmp ) ;
    vProb = 1 ./ ( exp( -vTmp ) + sum( mTmp, 2 ) ) ;   % exp( -vTmp ) is for the outside option
    
    PROB1( vCurRange ) = PROB( vCurRange ) + vProb * Q_WEIGHTS( ixNode ) ;
  end
  
end

%}

%% faster but numerically stable

% we proceed by computing choice probabilities for each product
% Remember:  Data is organized as ( prod 1 in markets 1 .. T ; prod 2 in
% markets 1 .. T ; etc. )

global Q_WEIGHTS T prods

nNodes = length( Q_WEIGHTS ) ;
PROB   = zeros( length( vMeanUtil ), 1 ) ;

mUtil = repmat( vMeanUtil, 1, nNodes ) + mMu ;

for ixProd = 1 : prods
  
  % Extract utility for current product
  vCurRange   = ( ixProd - 1 ) * T + ( 1 : T )' ;
  mChoiceUtil = repmat( mUtil( vCurRange, : ), prods, 1 ) ;
  mTmp        = mUtil - mChoiceUtil ;
  mSSj        = 1 ./ ( sharesum * exp( mTmp ) + exp( - mUtil( vCurRange, : ) ) ) ;    % Last term is outside good
  
  PROB( vCurRange ) = mSSj * Q_WEIGHTS ;
end


% % keyboard ;
