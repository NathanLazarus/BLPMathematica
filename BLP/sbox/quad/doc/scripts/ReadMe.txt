To produce tables for paper:

1) Extract data from solver log (Extract.sed ExtractPointEst.sed):

sed -n -f Extract*.sed SolverLog.txt > Out.txt

2) Edit Out.txt to remove extraneous text and prepend seed number to 
   each line

For Extract.sed, then run

  gawk -f Join.awk Out.txt

to perform necessary edits.  For ExtractPointEst.sed, use JoinPointEst.awk.

Note:  You can string it all together with:

sed -n -f Extract*.sed SolverLog.txt | gawk -f Join*.awk > Out.txt

3) Create tables using R scripts.

> source( 'MakePretty.R' )
> MakePretty( 'pMC' )
> MakePretty( 'SGI' )
> MakePretty( 'GH' )
> MakePointEst( 'pMC' )
> MakePointEst( 'SGI' )
> MakePointEst( 'GH' )
>
> mStdErr <- LoadStdErr( szBaseDir, nQuadType, nNodes )
> MakePointEstStdErr( 'GH', 'Some Caption', mStdErr )

4) Make short versions for talk

> MakePointEstShort( <rule>, <dataset> )
