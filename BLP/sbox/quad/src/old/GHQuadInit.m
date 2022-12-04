% GHQuadInit - Create Gaussian-Hermite nodes and weights for integration
%
% Compute nodes and weights for numerical integration with a weight
% function Exp[ -x^2 ].  For use with the normal density, this means that
% the change of variables causes the det( 2 Var )^(1/2) to cancel the same
% factor whose reciprical multiplies the density.
%
%
% OUTPUT (in global variables)
%   Q_POS         Scaled nodes ( N_DIM x N_NODES^NDIM )
%   Q_WEIGHTS     Weights for quadrature node ( N_NODES^N_DIM ) scaled by
%                 pi^(N_DIM/2).
%   N_NODES       Number of nodes to use
%
% WARNING: Q_WEIGHTS is scaled by pi^(5/2) because I assume that you are
% integrating a normal distribution versus something.  I have also rescaled 
% the nodes by sqrt( 2 ) to convert from Gaussian to Normal density.
%
% modification history
% --------------------
% 19aug2010 bss updated paths for new tree structure.
% 29jan2010 bss fixed scaling bug in conversion from Gaussian to Normal density.
% 21oct2009 bss written.
%

function [ ] = GHQuadInit( vMean_, mVar_, nQuadNodes_ )


%% Get quadrature weights and nodes

global N_NODES   ;  % How many quadrature nodes
global N_DIM     ;  % How many dimensions we are integrating over
global Q_NODES   ;  % location of nodes in transformed coordinates
global Q_WEIGHTS ;  % weights for nodes
% % global Q_POS     ;  % location of nodes in untransformed coordinates


% Number of dimensions to integrate over
N_DIM = length( vMean_ ) ;

% How many nodes to use for integral approximation
if nargin < 3
  N_NODES = 5 ;
else
  N_NODES = nQuadNodes_ ;
end

addpath( '~/Tools/matlab' ) ;
tmp = gauher( N_NODES ) ;

% extract quadrature information for one dimension
vNodes   = tmp( :, 1 ) ;
vWeights = tmp( :, 2 ) ;

% calculate three dimensional nodes and weights
% Q_WEIGHTS = kron( vWeights, kron( vWeights, vWeights ) ) ;

Q_WEIGHTS = vWeights ;
for ix = 2 : N_DIM
  Q_WEIGHTS = kron( vWeights, Q_WEIGHTS ) ; 
end


% WARNING!!!
% Correct for integration versus normal density XXX
Q_WEIGHTS = Q_WEIGHTS / ( pi ^ ( N_DIM / 2 ) ) ;


% Calculate node vector for each location and log income level
% Q_NODES = [ kron( vNodes, ones( N_NODES * N_NODES, 1 ) )                    ... % first column: XLoc
%             kron( ones( N_NODES, 1 ), kron( vNodes, ones( N_NODES, 1 ) ) )  ... % second column: YLoc
%             kron( ones( N_NODES * N_NODES, 1 ), vNodes )                    ... % third column: log income
%           ] ;


% This is some slickness to make sure that the right-most dimension (ixDim = N_DIM) varies
% most quickly and the left-most (ixDim = 1), most slowly
Q_NODES = zeros( N_DIM, N_NODES^N_DIM ) ;
for ixDim = 1 : N_DIM
  Q_NODES( ixDim, : ) = kron( ones( N_NODES^(ixDim - 1), 1 ), kron( vNodes, ones( N_NODES^(N_DIM - ixDim), 1 ) ) ) ;
end

%% Calculate distance from stores to quadrature nodes
%
% Note the presence of the 'mysterious' factor of two.  This is
% required because the quadrature formula works for a weighting
% function exp( - \sum_{i=1}^{3} x_i^2 ), so when transforming
% to scaleless coordinates, we must get rid of the 2 in the normal
% pdfs.

% % U = chol( 2 * mVar_ ) ;

% converst normal mode dimensions to correlated position: 
% row is [xloc; yloc; income ], column is a different node ...
% % Q_POS = U * Q_NODES' + repmat( vMean_, 1, N_NODES^N_DIM ) ;
% Q_POS = Q_POS' ;  % now each row corresponds to a node and the columns are the correctly scaled integration nodes


% Correct Q_NODES for conversion from exp( -x^2 ) kernel to exp( -0.5 * x^2 ) Normal Density
Q_NODES = Q_NODES * sqrt( 2 ) ;


