#!/bin/bash
#
# Script to convert solver logs into raw input for R to turn into LaTeX tables
#
# modification history
# --------------------
# 16aug2010 bss updated for refactored directory layout.
# 18jun2010 bss added ND for case of one start, multiple draws.
# 18may2010 bss written.
#

OUTDIR=/Users/bss/sbox/quad/doc/tables

# List of files
for file in  GH5        \
             SGI        \
             pMC        \
             pMC-R10000 \ 
             ND         \
             Monomial
do
  # Create input for solver results
  sed -n -f Extract.sed ${OUTDIR}/SolverLog-${file}.txt           \
    | gawk -f Join.awk > ${OUTDIR}/SolverResults-${file}.txt

  # Create input for point estimates
  sed -n -f ExtractPointEst.sed ${OUTDIR}/SolverLog-${file}.txt    \
    | gawk -f JoinPointEst.awk > ${OUTDIR}/PointEst-${file}.txt
done
