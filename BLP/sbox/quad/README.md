# Skrainka-Judd Quadrature

This repository contains the code and paper for Skrainka-Judd 
[High performance quadrature rules: How numerical integration affects a popular model of product differentiation](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1870703)

##  Directory Layout

There are two main directories, `doc` and `src`. The former contains the resources to build the paper, including
data analysis and preparation scripts. The latter contains the source code to run simulations. Start with `Driver.m`
which controls and runs the simulation.

The simulation code is currently based on the first version of the 
[Dubé, Fox, and Su](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA8585)
MPEC [code](https://www.jp-dube.com/research/MPECcode.html). Consequently, our code also uses their 
original ordering of *products* then *markets*, which produces a Hessian which is not block diagonal.
Our code needs to be updated to their latest version which does use the correct *markets* then *products*
ordering.

When you run the simulation, it will store results in `data`.

```
$ tree -L 2 -d
.
├── data
│   └── DataJ25T50
├── doc
│   ├── Publication
│   ├── fig
│   ├── scripts
│   └── tables
└── src
```

The code `BLPSetup.m` will create the necessary synthetic datasets, storing them under the `data`
directory using names like `DataJjjTtt` where `jj` and `tt` refer to the number of products and markets, 
respectively.

##  Requirements

This code requires:

* MATLAB
* Tomlab + KNITRO (and, possibly, SNOPT)
* R
* LyX

Tomlab could easily be replaced with the native interfaces for KNITRO and SNOPT.

## Configuration

The code is highly configurable to support a variety of quadrature rules and analyses. To adjust
the data generating process, see the file `BLPSetup.m`. The solver is tuned and run using 
`runGMMMPECKnitro.m`; `knitroOptions2.opt` controls addition KNITRO tuning parameters. You
should adjust these settings for best performance.

`Driver.m` contains is written using MATLAB cells so that you can easily run independent parts of
estimation and analysis. In addition, it contains additional configuration settings such as
the number of replications, products,  and markets, and the quadrature rule to use.

To select different quatrature rules, set the variable `nQuadType`. Supported values are:

| `nQuadType`   |  Quadrature Rule  |
|-------------- |-------------------|
| 0 (default)   | Psuedo-Monte Carlo  |
| 1             | Gauss-Hermite Product Rule  |
| 2             | Stroud Monomial Rule 11-1 (five dimensional problems with Gaussian Kernel)  |
| 3             | Quasi-Monte Carlo |
| 4             | Sparse grids |


Perform additional tuning for these rules `ComputeGMMEstimates.m`.

##  Running the models

To run the models and much of the analysis, run `Driver.m` after you have adjusted the settings for
your situation. Because solving large BLP problems is very time consuming, start with smaller problems with
only one or two replications so that you can fail fast and debug any issues more quickly.
