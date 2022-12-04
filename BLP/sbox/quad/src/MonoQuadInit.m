% MonoQuadInit - Create nodes and weights for integration using monomial rule
%
% Computes nodes and weights based on monomial rule in Stroud, p. 322,
% 11-1.  This rule is exact for all degree 11 monomials in 5 dimensions.
%
% INPUT
%   nRule_        Which Stroud rule to use
% OUTPUT (in global variables)
%   Q_POS         Scaled nodes ( N_DIM x N_NODES^NDIM )
%   Q_WEIGHTS     Weights for quadrature node ( N_NODES^N_DIM ) scaled by
%                 pi^(N_DIM/2).
%   N_NODES       Number of nodes to use
%
% WARNING: I rescale by pi ^ ( N_DIM / 2 ) because I assume you are
% integrating a normal density!!!  I have also rescaled the nodes by sqrt( 2 )
% to convert from Gaussian to Normal density.
%
% modification history
% --------------------
% 29jan2010 bss fixed scaling bug in conversion from Gaussian to Normal density.
% 29oct2009 bss written.
%

function [ ] = MonoQuadInit( nWhichRule_ )


%% Setup

global N_NODES   ;  % How many quadrature nodes
global N_DIM     ;  % How many dimensions we are integrating over
global Q_NODES   ;  % location of nodes in transformed coordinates
global Q_WEIGHTS ;  % weights for nodes
% % global Q_POS     ;  % location of nodes in untransformed coordinates


% Number of dimensions to integrate over
N_DIM   = 5 ;   % Assume we are using 5 dimensions.  if not => trouble!

N_NODES = 983 ;    % Total number of nodes

% extract quadrature information for one dimension
Q_NODES   = zeros( N_NODES, N_DIM ) ;
Q_WEIGHTS = zeros( N_NODES, 1 ) ;

%% Stroud setup
%   I use Stroud's notation as much as possible to keep things clean
%
%----------------------------
%
% Monomial Rule 11-1, pp. 322-3
%
% (0, 0, 0, 0, 0, 0, ... ,0)        B0      0
% (u, 0, 0, 0, 0, 0, ... , 0)   FS  B1      1 u
% (v, 0, 0, 0, 0, 0, ... ,O)    FS  B2      1 v 
% (w, 0, 0, 0, 0, 0, ... , 0)   FS  B3      1 w 
% (u, u, 0, 0, 0, 0, ... , 0)   FS  B4      2 u
% (v, v, 0, 0, 0, 0, ... ,O)    FS  B5      2 v
% (w, w, 0, 0, 0, 0, ... ,O)    FS  B6      2 w 
% (u, v, 0, 0, 0, 0, ... , 0)   FS  B7      1u 1v
% (u, w, 0, 0, 0, 0, ... , O)   FS  B8      1u 1w
% (u, u, u, 0, 0, 0 , ... , 0)  FS  B9      3u
% (v, v, v, 0, 0, 0, ... ,O)    FS  B10     3v
% (w, W, W, 0, 0, 0, ... , 0)   FS  B11     3w
% (u, U, v, 0, 0, 0, ... , 0)   FS  B12     2u 1v
% (u, U, u, U, 0, 0, ... , 0)   FS  B13     4u
% (v, v, v, v, 0, 0, ... ,0)    FS  B14     4v
% (u, u, u, U, u, 0, ... ,0)    FS  B15     5u
%
% FS == fully symmetric, i.e. all permutations and changes of signs
%
%----------------------------
%
% Nodes and Weights
%
% I list both of Stroud's solutions below
%
% n = 5, N = 983 
%
%   (1)0.23506 04973 67449 
%      0.43607 74119 27617 
%   (1)0.13358 49074 01370 
%   (4)0.25588 52693 11763 
%     -0.43959 86774 91526 
%  -(4)0.10654 14061 44610 
%      0.45354 09090 54264 
% -(-1)0.13210 09056 23778 
%   (3)0.41860 65689 54203 
%  (-1)0.51139 45630 43680 
%  (-1)0.64558 10138 45604 
% -(-2)0.53341 72774 94500 
% -(-3)0.13798 16262 54496 
%  -(3)0.14743 69331 89884 
%  (-1)0.30425 38077 65057 
%  (-2)0.24810 86982 07828 
%  (-4)0.11365 20945 46015 
%   (2)0.39425 74071 60391 
%  (-5)0.33172 50113 58320 
% 
% 
% ------------------------------------------------------------
% Other solution
% 
%   (1)0.23506 04973 67449 
%   (1)0.13358 49074 01370 
%      0.43607 74119 27617 
%  -(3)0.76130 53475 48192 
%     -0.53636 08050 19297 
%      0.11066 98320 78736 
%   (3)0.24642 10889 23968 
% -(-3)0.77364 93279 68607 
%      0.16908 86412 05970 
%  -(2)0.67070 06802 43651 
% -(-2)0.85609 05602 29205 
%  (-1)0.79444 62327 70302 
% -(-3)0.22027 28632 63544 
% -(-2)0.37351 58122 28225 
%   (2)0.12360 65440 52884 
%  (-3)0.53778 88045 57843 
% -(-4)0.12210 18614 80881 
%  (-2)0.72511 70707 59373 
%  (-5)0.33172 50113 58320 

switch nWhichRule_

  case 1
    vStroudData = [                                 ...
                       0.235060497367449e1	;			...   % u
                       0.436077411927617		;			...   % v
                       0.133584907401370e1	;			...   % w
                       0.255885269311763e4	;			...   % B0
                      -0.439598677491526		;			...   % B1
                      -0.106541406144610e4	;			...
                       0.453540909054264		;			...
                      -0.132100905623778e-1	;			...
                       0.418606568954203e3	;			...
                       0.511394563043680e-1	;			...
                       0.645581013845604e-1	;			...
                      -0.533417277494500e-2	;			...
                      -0.137981626254496e-3	;			...
                      -0.147436933189884e3	;			...
                       0.304253807765057e-1	;			...
                       0.248108698207828e-2	;			...
                       0.113652094546015e-4	;			...
                       0.394257407160391e2	;			...
                       0.331725011358320e-5  			...   % B15
                  ] ;
              
  case 2
    vStroudData = [                                 ...
                       0.235060497367449e1			;	...
                       0.133584907401370e1			;	...
                       0.436077411927617				;	...
                      -0.761305347548192e3			;	...
                      -0.536360805019297				;	...
                       0.110669832078736				;	...
                       0.246421088923968e3			;	...
                      -0.773649327968607e-3			;	...
                       0.169088641205970				;	...
                      -0.670700680243651e2			;	...
                      -0.856090560229205e-2			;	...
                       0.794446232770302e-1			;	...
                      -0.220272863263544e-3			;	...
                      -0.373515812228225e-2			;	...
                       0.123606544052884e2			;	...
                       0.537788804557843e-3			;	...
                      -0.122101861480881e-4			;	...
                       0.725117070759373e-2			;	...
                       0.331725011358320e-5			 	...          
                  ] ;
  otherwise
    error( 'Unknown Stroud rule' ) ;    
end

% Setup Index for Stroud's data
IStroud.u   =  1 ;
IStroud.v   =  2 ;
IStroud.w   =  3 ;
IStroud.B0	=  4 ;
IStroud.B1	=  5 ;
IStroud.B2	=  6 ;
IStroud.B3	=  7 ;
IStroud.B4	=  8 ;
IStroud.B5	=  9 ;
IStroud.B6	= 10 ;
IStroud.B7	= 11 ;
IStroud.B8	= 12 ;
IStroud.B9	= 13 ;
IStroud.B10	= 14 ;
IStroud.B11	= 15 ;
IStroud.B12	= 16 ;
IStroud.B13	= 17 ;
IStroud.B14	= 18 ;
IStroud.B15	= 19 ;

ixStart = 1 ;
ixEnd   = ixStart ;


%% Begin computing nodes and weights

% B0
Q_WEIGHTS( ixStart ) = vStroudData( IStroud.B0 ) ;
Q_NODES( ixStart )   = 0 ;

% B1 ( u, 0, ... ) FS
nNodes  = nchoosek( N_DIM, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

mTmp = unique( perms( [ vStroudData( IStroud.u ) zeros( 1, N_DIM - 1 ) ] ), 'rows' ) ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B1 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B1 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;


% B2 ( v, 0, ... ) FS
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

mTmp = unique( perms( [ vStroudData( IStroud.v ) zeros( 1, N_DIM - 1 ) ] ), 'rows' ) ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B2 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B2 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;


% B3 ( w, 0, ... ) FS
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

mTmp = unique( perms( [ vStroudData( IStroud.w ) zeros( 1, N_DIM - 1 ) ] ), 'rows' ) ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B3 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B3 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;


% B4 ( u, u, 0, ... ) FS 
nNodes  = nchoosek( N_DIM, 2 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

mTmp = unique( perms( [ vStroudData( IStroud.u ) vStroudData( IStroud.u ) zeros( 1, N_DIM - 2) ] ), 'rows' ) ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B4 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B4 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;

mTmp = unique( perms( [ vStroudData( IStroud.u ) -vStroudData( IStroud.u ) zeros( 1, N_DIM - 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B4 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;


% B5 ( v, v, 0, ... ) FS
mTmp = unique( perms( [ vStroudData( IStroud.v ) vStroudData( IStroud.v ) zeros( 1, N_DIM - 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B5 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B5 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;

mTmp = unique( perms( [ vStroudData( IStroud.v ) -vStroudData( IStroud.v ) zeros( 1, N_DIM - 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B5 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;


% B6 ( w, w, 0, ... ) FS
mTmp = unique( perms( [ vStroudData( IStroud.w ) vStroudData( IStroud.w ) zeros( 1, N_DIM - 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B6 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B6 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;

mTmp = unique( perms( [ vStroudData( IStroud.w ) -vStroudData( IStroud.w ) zeros( 1, N_DIM - 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B6 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;


% B7 ( u, v, 0, ... ) FS
nNodes  = N_DIM * ( N_DIM - 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

mTmp = unique( perms( [ vStroudData( IStroud.u ) vStroudData( IStroud.v ) zeros( 1, N_DIM - 2) ] ), 'rows' ) ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B7 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B7 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;

mTmp = unique( perms( [ vStroudData( IStroud.u ) -vStroudData( IStroud.v ) zeros( 1, N_DIM - 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B7 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

mTmp = unique( perms( [ -vStroudData( IStroud.u ) vStroudData( IStroud.v ) zeros( 1, N_DIM - 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B7 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;


% B8 ( u, w, 0, ... ) FS
nNodes  = N_DIM * ( N_DIM - 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

mTmp = unique( perms( [ vStroudData( IStroud.u ) vStroudData( IStroud.w ) zeros( 1, N_DIM - 2) ] ), 'rows' ) ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B8 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B8 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;

mTmp = unique( perms( [ vStroudData( IStroud.u ) -vStroudData( IStroud.w ) zeros( 1, N_DIM - 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B8 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

mTmp = unique( perms( [ -vStroudData( IStroud.u ) vStroudData( IStroud.w ) zeros( 1, N_DIM - 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B8 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;


% B9 ( u, u, u, 0, ... ) FS
mTmp = unique( perms( [ vStroudData( IStroud.u ) vStroudData( IStroud.u ) vStroudData( IStroud.u ) zeros( 1, 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B9 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B9 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;

mTmp = unique( perms( [ -vStroudData( IStroud.u ) vStroudData( IStroud.u ) vStroudData( IStroud.u ) zeros( 1, 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B9 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B9 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;


% B10 ( v, v, v, 0, ... ) FS
mTmp = unique( perms( [ vStroudData( IStroud.v ) vStroudData( IStroud.v ) vStroudData( IStroud.v ) zeros( 1, 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B10 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B10 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;

mTmp = unique( perms( [ -vStroudData( IStroud.v ) vStroudData( IStroud.v ) vStroudData( IStroud.v ) zeros( 1, 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B10 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B10 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;


% B11 ( w, w, w, 0, ... ) FS
mTmp = unique( perms( [ vStroudData( IStroud.w ) vStroudData( IStroud.w ) vStroudData( IStroud.w ) zeros( 1, 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B11 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B11 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;

mTmp = unique( perms( [ -vStroudData( IStroud.w ) vStroudData( IStroud.w ) vStroudData( IStroud.w ) zeros( 1, 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B11 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B11 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;


% B12 ( u, u, v, 0, ... ) FS
mTmp = unique( perms( [ vStroudData( IStroud.u ) vStroudData( IStroud.u ) vStroudData( IStroud.v ) zeros( 1, 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B12 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B12 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;

mTmp = unique( perms( [ -vStroudData( IStroud.u ) vStroudData( IStroud.u ) vStroudData( IStroud.v ) zeros( 1, 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B12 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;


mTmp = unique( perms( [ -vStroudData( IStroud.u ) -vStroudData( IStroud.u ) vStroudData( IStroud.v ) zeros( 1, 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B12 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

mTmp = unique( perms( [ -vStroudData( IStroud.u ) vStroudData( IStroud.u ) -vStroudData( IStroud.v ) zeros( 1, 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B12 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

mTmp = unique( perms( [ vStroudData( IStroud.u ) vStroudData( IStroud.u ) -vStroudData( IStroud.v ) zeros( 1, 2) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B12 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;


% B13 ( u, u, u, u, 0 ) FS
mTmp = unique( perms( [ vStroudData( IStroud.u ) * ones( 1, 4 ),  0 ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B13 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B13 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;

mTmp = unique( perms( [ -vStroudData( IStroud.u ) vStroudData( IStroud.u ) * ones( 1, 3 ),  0 ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B13 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B13 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;


mTmp = unique( perms( [ -vStroudData( IStroud.u ) * ones( 1, 2 )  vStroudData( IStroud.u ) * ones( 1, 2 ),  0 ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B13 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;


% B14 ( v, v, v, v, 0 ) FS
mTmp = unique( perms( [ vStroudData( IStroud.v ) * ones( 1, 4 ),  0 ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B14 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B14 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;

mTmp = unique( perms( [ -vStroudData( IStroud.v )  vStroudData( IStroud.v ) * ones( 1, 3 ),  0 ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B14 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B14 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;

mTmp = unique( perms( [ -vStroudData( IStroud.v ) * ones( 1, 2 )  vStroudData( IStroud.v ) * ones( 1, 2 ),  0 ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B14 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;


% B15 ( u, u, u, u, u ) FS
mTmp = vStroudData( IStroud.u ) * ones( 1, 5 ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B15 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B15 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;

mTmp = unique( perms( [ -vStroudData( IStroud.u )  vStroudData( IStroud.u ) * ones( 1, 4 ) ] ), 'rows' )  ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B15 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B15 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;

mTmp = unique( perms( [ -vStroudData( IStroud.u ) * ones( 1, 2 )  vStroudData( IStroud.u ) * ones( 1, 3 ) ] ), 'rows' ) ;

nNodes  = size( mTmp, 1 ) ;
ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B15 ) ;
Q_NODES( ixStart : ixEnd, : ) = mTmp ;

ixStart = ixEnd + 1 ;
ixEnd   = ixStart + nNodes - 1 ;

Q_WEIGHTS( ixStart : ixEnd )  = vStroudData( IStroud.B15 ) ;
Q_NODES( ixStart : ixEnd, : ) = -mTmp ;


%% Rescale for use with Normal Density

Q_WEIGHTS = Q_WEIGHTS / (  pi ^ ( N_DIM / 2 ) ) ;    % Rescale for normal density

Q_NODES   = Q_NODES' * sqrt( 2 ) ;    % To agree with JP's code and convert to normal density

