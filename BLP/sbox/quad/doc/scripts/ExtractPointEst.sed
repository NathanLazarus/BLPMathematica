# ExtractPointEst - extracts point estimates from solver log

s/Data.*Seed.*\([0-9]\{4\}\)/\1/p
#/Optimal x/,/State vector hs/p 
/Optimal x/,/State vector hs/s/\(x:\)\{0,1\}\(.*\)/\2/p 
