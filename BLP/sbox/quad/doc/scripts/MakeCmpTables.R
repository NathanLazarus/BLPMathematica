# MakeCmpTables - makes tables to compare results for different integration rules
#
# modification history
# --------------------
# 26oct2010 bss written.
#


#-----------------------------------------------------------
# Setup
library( xtable )

setwd('/Users/bss/sbox/quad/doc/tables')


#===============================================================================
# Create nice LaTeX long tables to compare solver results
# These tables are for the paper.  Shorter versions for the
# talk are created below.
#===============================================================================

# Setup pathnames to data

szIn.pMC.R01000 <- 'SolverResults-pMC5Good-R01000.txt'
szIn.pMC.R10000 <- 'SolverResults-pMC-R10000.txt'
szIn.GH         <- 'SolverResults-GH7-5Good.txt'
szIn.SGI        <- 'SolverResults-SGI5Good.txt'
szIn.Mono       <- 'SolverResults-Mono5Good.txt'  


szLabel <- 'TabCmpSolResults.f_k'
szCap   <- "Objective function values (f\\_k) at optimum by rule."


#---------------------------------------
# Load data
df.pMC.R01000 <- read.table( file = szIn.pMC.R01000 )
df.pMC.R10000 <- read.table( file = szIn.pMC.R10000 )
df.GH         <- read.table( file = szIn.GH )
df.SGI        <- read.table( file = szIn.SGI )
df.Mono       <- read.table( file = szIn.Mono )

# reorder rows
#sol.df <- cbind( df[,1], df[,c(3,4)], df[,2], df[,5] )
#colnames( sol.df ) <- c( 'Dataset', 'EXIT', 'INFORM', 'f_k', 'CPU Time' )

#-----------------------------
#-----------------------------
# Compare f_k
df.sol <- cbind( df.GH[,1], df.pMC.R01000[,2], df.pMC.R10000[,2],
                 df.GH[,2], df.SGI[,2], df.Mono[,2] )
colnames( df.sol ) <- c( 'Dataset', 'pMC R=1,000', 'pMC R=10,000', 
                         "GH $7^5$", "SGI", "Monomial" )

xt <- xtable( df.sol, caption=szCap, label=szLabel, digits=c(0,0,5,5,5,5,5) )

szOut <- 'CmpSolResults.f_k.inc'

print( xt, file=szOut, type='latex', math.style.negative=T, 
      include.rownames=F, floating=T, caption.placement='bottom',
      hline.after=c( -1, 0, seq(5,25,5) ), 
      sanitize.text.function=function(x){ x } )

# Make Short Tables by rule and dataset 

#-------------------
# f_k

# Dataset 1
# GH vs. pMC
szCap <- "Obj. function (f\\_k) at optimum for Dataset 1 by rule."
szOut <- 'CmpSolResults.f_k.Short.pMC.1.inc'

xt <- xtable( df.sol[seq(1,5),c(1,4,2,3)], caption=szCap, digits=c(0,0,5,5,5) )

print( xt, file=szOut, type='latex', math.style.negative=T, 
      include.rownames=F, floating=T, caption.placement='bottom',
      hline.after=c( -1, 0, seq(5,5,5) ), 
      sanitize.text.function=function(x){ x } )

# GH vs. SGI vs. Monomial
szCap <- "Obj. function (f\\_k) at optimum for Dataset 1 by rule."
szOut <- 'CmpSolResults.f_k.Short.poly.1.inc'

xt <- xtable( df.sol[seq(1,5),c(1,4,5,6)], caption=szCap, digits=c(0,0,5,5,5) )

print( xt, file=szOut, type='latex', math.style.negative=T, 
      include.rownames=F, floating=T, caption.placement='bottom',
      hline.after=c( -1, 0, seq(5,5,5) ), 
      sanitize.text.function=function(x){ x } )


# Dataset 5
# GH vs. pMC
szCap <- "Obj. function (f\\_k) at optimum for Dataset 5 by rule."
szOut <- 'CmpSolResults.f_k.Short.pMC.5.inc'

xt <- xtable( df.sol[seq(21,25),c(1,4,2,3)], caption=szCap, digits=c(0,0,5,5,5) )

print( xt, file=szOut, type='latex', math.style.negative=T, 
      include.rownames=F, floating=T, caption.placement='bottom',
      hline.after=c( -1, 0, seq(5,5,5) ), 
      sanitize.text.function=function(x){ x } )

# GH vs. SGI vs. Monomial
szCap <- "Obj. function (f\\_k) at optimum for Dataset 5 by rule."
szOut <- 'CmpSolResults.f_k.Short.poly.5.inc'

xt <- xtable( df.sol[seq(21,25),c(1,4,5,6)], caption=szCap, digits=c(0,0,5,5,5) )

print( xt, file=szOut, type='latex', math.style.negative=T, 
      include.rownames=F, floating=T, caption.placement='bottom',
      hline.after=c( -1, 0, seq(5,5,5) ), 
      sanitize.text.function=function(x){ x } )





#-----------------------------
#-----------------------------
# Compare CPU Time
df.sol <- cbind( df.GH[,1], df.pMC.R01000[,5], df.pMC.R10000[,5],
                 df.GH[,5], df.SGI[,5], df.Mono[,5] )
colnames( df.sol ) <- c( 'Dataset', 'pMC R=1,000', 'pMC R=10,000', 
                         "GH $7^5$", "SGI", "Monomial" )
szCap   <- "CPU time required by SNOPT by rule."
szLabel <- 'TabCmpSolResults.CPUTime'
szOut   <- 'CmpSolResults.CPUTime.inc'

xt <- xtable( df.sol, caption=szCap, label=szLabel, digits=c(0,0,2,2,2,2,2) )

print( xt, file=szOut, type='latex', math.style.negative=T, 
      include.rownames=F, floating=T, caption.placement='bottom',
      hline.after=c( -1, 0, seq(5,25,5) ), 
      sanitize.text.function=function(x){ x } )
  

# Make Short Tables by rule and dataset 

#-------------------
# CPU Time

# Dataset 1
# GH vs. pMC
szCap <- "CPU time required by SNOPT for Dataset 1 by rule."
szOut <- 'CmpSolResults.CPUTime.Short.pMC.1.inc'

xt <- xtable( df.sol[seq(1,5),c(1,4,2,3)], caption=szCap, digits=c(0,0,2,2,2) )

print( xt, file=szOut, type='latex', math.style.negative=T, 
      include.rownames=F, floating=T, caption.placement='bottom',
      hline.after=c( -1, 0, seq(5,5,5) ), 
      sanitize.text.function=function(x){ x } )

# GH vs. SGI vs. Monomial
szCap <- "CPU time required by SNOPT for Dataset 1 by rule."
szOut <- 'CmpSolResults.CPUTime.Short.poly.1.inc'

xt <- xtable( df.sol[seq(1,5),c(1,4,5,6)], caption=szCap, digits=c(0,0,2,2,2) )

print( xt, file=szOut, type='latex', math.style.negative=T, 
      include.rownames=F, floating=T, caption.placement='bottom',
      hline.after=c( -1, 0, seq(5,5,5) ), 
      sanitize.text.function=function(x){ x } )


# Dataset 5
# GH vs. pMC
szCap <- "CPU time required by SNOPT for Dataset 5 by rule."
szOut <- 'CmpSolResults.CPUTime.Short.pMC.5.inc'

xt <- xtable( df.sol[seq(21,25),c(1,4,2,3)], caption=szCap, digits=c(0,0,2,2,2) )

print( xt, file=szOut, type='latex', math.style.negative=T, 
      include.rownames=F, floating=T, caption.placement='bottom',
      hline.after=c( -1, 0, seq(5,5,5) ), 
      sanitize.text.function=function(x){ x } )

# GH vs. SGI vs. Monomial
szCap <- "CPU time required by SNOPT for Dataset 5 by rule."
szOut <- 'CmpSolResults.CPUTime.Short.poly.5.inc'

xt <- xtable( df.sol[seq(21,25),c(1,4,5,6)], caption=szCap, digits=c(0,0,2,2,2) )

print( xt, file=szOut, type='latex', math.style.negative=T, 
      include.rownames=F, floating=T, caption.placement='bottom',
      hline.after=c( -1, 0, seq(5,5,5) ), 
      sanitize.text.function=function(x){ x } )

#===============================================================================
# New Plots to compare rules
#===============================================================================

df.cpu <- cbind( df.Mono[,5], df.SGI[,5], df.pMC.R01000[,5], df.pMC.R10000[,5], df.GH[,5] )
colnames( df.cpu ) <- c( 'Monomial', 'SGI', 'pMC R=1,000', 'pMC R=10,000', 'GH 7^5' )

szOut <- 'Cmp.CPU.Time.Box.png'

quartz( file=szOut, type='png', bg='white', height=8, width=10 )
boxplot( df.cpu, xlab='Rule', ylab='CPU Time', main='CPU Time vs. Rule' )
dev.off()


