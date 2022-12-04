# RulesTable - generate LaTeX table comparing rules.
#
# modification history
# --------------------
# 02dec2011 bss added column label f(x) to monomial rule comparison and
#               hline after exact rules.
# 27jan2011 bss written.
#

#---------------------------------------------------------------------
# LoadRule - loads comparison results
#

LoadRule <- function( szDataFile ) {

  require( R.matlab )

  mCmp <- readMat( szDataFile )

  return( mCmp ) 
}


#---------------------------------------------------------------------
# ComputeNodes -


#---------------------------------------------------------------------
# Setup

setwd('/Users/bss/sbox/quad/doc/scripts/')

library( xtable )

# szSGI   <- '~/sbox/blpcpp/data/quad/SGIRule.5-6.csv'
# szMonoL <- '~/sbox/blpcpp/data/quad/Mono.11-1.L.csv'
# szMonoR <- '~/sbox/blpcpp/data/quad/Mono.11-1.R.csv'
# szGH3   <- '~/sbox/blpcpp/data/quad/GHProdRule.3-5.csv'
# szGH5   <- '~/sbox/blpcpp/data/quad/GHProdRule.5-5.csv'
# szGH7   <- '~/sbox/blpcpp/data/quad/GHProdRule.7-5.csv'
# szPMC   <- '~/sbox/blpcpp/data/quad/pMCRule.10000.csv'

szCmpDataFile <- '/Users/bss/sbox/quad/doc/tables/CmpRules.mat' 
szOutDir      <- '/Users/bss/sbox/quad/doc/tables'

nDigits       <- 2

#---------------------------------------------------------------------
# Load data and label

mCmp <- LoadRule( szCmpDataFile )

mCmpErr <- mCmp$mCmpErr
mQCmpErr <- mCmp$mQCmpErr
mCases  <- mCmp$Cases

#colnames( mCmpErr ) <- c( "GH $3^5$", "GH $5^5$", "GH $7^5$", 
#                       "Sparse Grid", "11-1 L", "11-1 R", 
#                       "$\\overline{|pMC|}$", "$q25 |pMC|$", "$q50 |pMC|$", 
#                       "$q75 |pMC|$" )
#
#rownames( mCmpErr ) <- c( "1", "$x_1$", "$x_1*x_2$", "$x_1*x_2*x_3*x_4*x_5$",
#                       "$x_1^2$", "$x_1^4$", "$x_1^6$", "$x_2^6*x_4^4$" )

colnames( mCmpErr ) <- c( "GH $3^5$", "GH $5^5$", "GH $7^5$", 
                       "Sparse Grid", "11-1 L", "11-1 R", 
                       "$\\overline{|pMC|}$", "$\\sigma(pMC)$" )

colnames( mQCmpErr ) <- c( "GH $3^5$", "GH $5^5$", "GH $7^5$", 
                       "Sparse Grid", "11-1 L", "11-1 R", 
                       "$\\overline{|qMC|}$", "$\\sigma(qMC)$" )

#rownames( mCmpErr ) <- c( "1", "$x_1$", "$x_1*x_2$", "$x_1*x_2*x_3*x_4*x_5$",
#                       "$x_1^2$", "$x_1^4$", "$x_1^6$", "$x_2^6*x_4^4$" )

#-------------------
# Build rownames

nCases <- nrow( mCases )
nDim   <- ncol( mCases )

vRowNames <- list()
for( ixCase in seq( 1, nCases ) ) {

  szName <- "$"
  for( ixCol in seq( 1, nDim ) ) {

    # Get exponent for current term in monomial
    nExp <- mCases[ ixCase, ixCol ] 
    if( nExp > 0 ) {
      szName <- paste( szName, " x_{", ixCol, "}^{", nExp, "}", sep='' )
    }

  }
  szName <- paste( szName, "$", sep='' )

  if( "$$" == szName ) {
    szName <- "$1$"
  }

  vRowNames[ ixCase ] <- szName
}

rownames( mCmpErr ) <- vRowNames
rownames( mQCmpErr ) <- vRowNames

# ADDED to ensure f(x) label over rownames and will hopefully work!
mCmpErr <- data.frame( cbind( vRowNames, mCmpErr ) )
colnames( mCmpErr ) <- c( "$f\\left( x \\right)$", 
                       "GH $3^5$", "GH $5^5$", "GH $7^5$", 
                       "Sparse Grid", "11-1 L", "11-1 R", 
                       "$\\overline{|pMC|}$", "$\\sigma(pMC)$" )

#-------------------
# Create LaTeX tables

#----------
# pMC

szCap   <- 'Monomial Integration Errors: Polynomial Rules vs. pMC ($R=10,000$)'
szLabel <- 'Tab.Cmp.Rule'
vDisp   <- c( 's', rep( 'g', ncol( mCmpErr ) ) )

xt <- xtable( mCmpErr, caption=szCap, label=szLabel, 
              display=vDisp, digits=nDigits )

szOut <- paste( szOutDir, 'Cmp.Rule.inc', sep='/' )

# Commented out per ADDED comment above
#print( xt, file=szOut, type='latex', math.style.negative=T,
#      include.rownames=T, floating=T, caption.placement='bottom',
#      floating.environment='sidewaystable',
#      sanitize.text.function=function(x){ x } )

# Put lines around header, after exact monomials, and last row
vHLine <- c( -1, 0, 10, nrow( mCmpErr ) )

print( xt, file=szOut, type='latex', math.style.negative=T,
      include.rownames=F, floating=T, caption.placement='bottom',
      floating.environment='sidewaystable',
      hline.after=vHLine,
      sanitize.text.function=function(x){ x } )


#----------
# qMC

szCap   <- 'Monomial Integration Errors: Polynomial Rules vs. qMC ($R=10,000$)'
szLabel <- 'Tab.Cmp.Rule.qMC'
vDisp   <- c( 's', rep( 'g', ncol( mCmpErr ) ) )

xt <- xtable( mQCmpErr, caption=szCap, label=szLabel, 
              display=vDisp, digits=nDigits )

szOut <- paste( szOutDir, 'Cmp.Rule.qMC.inc', sep='/' )

print( xt, file=szOut, type='latex', math.style.negative=T,
      include.rownames=T, floating=T, caption.placement='bottom',
      floating.environment='sidewaystable',
      sanitize.text.function=function(x){ x } )


#----------
# pMC & qMC

szCap   <- 'Monomial Integration Errors: Polynomial Rules vs. pMC and qMC ($R=10,000$)'
szLabel <- 'Tab.Cmp.Rule.pMC.qMC'
vDisp   <- c( 's', rep( 'g', ncol( mCmpErr ) + 2 ) )

xt <- xtable( cbind( mCmpErr, mQCmpErr[,7:8] ), caption=szCap, label=szLabel, 
              display=vDisp, digits=nDigits )

szOut <- paste( szOutDir, 'Cmp.Rule.pMC.qMC.inc', sep='/' )

print( xt, file=szOut, type='latex', math.style.negative=T,
      include.rownames=T, floating=T, caption.placement='bottom',
      floating.environment='sidewaystable',
      sanitize.text.function=function(x){ x } )

#----------
# SGI vs. pMC - short for talk

szCap   <- 'Integration Errors: SGI vs. pMC ($R=10,000$)'
szLabel <- 'Tab.Cmp.Rule'

# Take a subset of rows & columns
vIxCol <- c( 4, 5, 7, 8 )
vIxRow <- seq( 1, 18 )

vHLine <- c( -1, 0, 10, max( vIxRow ) )

vDisp   <- c( 's', rep( 'g', ncol( mCmpErr[ vIxRow, vIxCol ] ) ) )

xt <- xtable( mCmpErr[ vIxRow, vIxCol ], caption=szCap, label=szLabel, 
              display=vDisp, digits=nDigits )

szOut <- paste( szOutDir, 'Cmp.Rule.Talk.inc', sep='/' )

print( xt, file=szOut, type='latex', math.style.negative=F,
      include.rownames=T, floating=F, caption.placement='bottom',
      hline.after=vHLine,
      sanitize.text.function=function(x){ x } )


