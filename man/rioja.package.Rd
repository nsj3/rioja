\name{rioja-package}
\alias{rioja-package}
\alias{rioja}
\docType{package}
\title{
Analysis of Quaternary Science Data
}
\description{
\code{rioja}: An R package for the analysis of Quaternary science data. Contains functions for constrained clustering, transfer functions, and plotting stratigraphic data.
}

\details{
\tabular{ll}{
Package: \tab rioja\cr
Type: \tab Package\cr
Version: \tab 0.9-15\cr
Date: \tab 2017-06-16\cr
License: \tab GPL-2\cr
}

Main functions:
\preformatted{
WA                      Weighted averaging transfer functions
WAPLS                   Partial least squares (PLS) and Weighted 
                        Averaging PLS transfer functions
MAT                     Modern analogue technique transfer functions
IKFA                    Imbrie and Kipp Factor Analysis
MLRC                    Maximum likelihood response curves
LWR                     Locally weighted transfer functions
MR                      Multiple regression transfer functions
crossval                Cross validate a transfer function model
predict                 Predict using a transfer function model, 
                        with sample-specific errors
randomPTF               Variable important via random palaeoecological 
                        transfer functions                        
read.CEP                Read and write data in Cornell Ecology 
                        Program format
read.Tilia              Read old format Tilia (*.til) files
strat.plot              Plot a stratigraphic diagram
chclust                 Constrained heirarchical clustering (zonation)
bstick                  Broken stick model to determine number of 
                        significant groups in a constrained cluster 
                        analysis
inkspot                 Plot a two-way ordered diagram of species 
                        abundances
Merge                   Merge two or more data frames by column names
interp.dataset          Interpolate one dataset to the depth or ages 
                        of another
compare.datasets        Assess conformity of species occurrence and 
                        abundance in two datasets
}

Datasets:
\preformatted{
aber                    Pollen stratigraphic data from the Abernethy 
                        forest
IK                      Imbrie and Kipp foraminifera data
Ponds                   Southeast England ponds and pools diatom 
                        dataset
SWAP                    Surface Waters Acidification Programme lake 
                        diatom and pH dataset
RLGH                    Diatom stratigraphic data from the Round Loch 
                        of Glenhead
}
}

\author{
Steve Juggins

Maintainer: Steve Juggins <Stephen.Juggins@ncl.ac.uk>
}
\keyword{ package }

