# Sed script to extract necessary results from solver log
#
# This script will extract the Seed ID followd by n pairs of lines where
# the first line is the f_k value of the objective function at the optimum
# and the second line is the EXIT and INFORM codes from SNOPT.
#
# modification history
# --------------------
# 10may2010 bss written.
#

s/Data.*Seed.*\([0-9]\{4\}\)/\1/p
s/.*f_k\(.*\)/\1/p
s/Solver: snopt.*EXIT=\([0-9]*\)\..*INFORM=\([0-9]*\).*/\1 \2/p
s/CPU time: \([0-9]*\.[0-9]*\).*/\1/p
