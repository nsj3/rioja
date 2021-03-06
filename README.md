rioja
========

Functions for the analysis of Quaternary science data, including constrained clustering, WA, WAPLS, IKFA, MLRC and MAT transfer functions, and stratigraphic diagrams.

``rioja`` implements a number of numerical methods for inferring the value of an environmental variable from a set of sepecies abundances, given a modern training set of species data and  associated environmental values.  In palaeoecology these are known as "transfer functions" or "inference models" and are used to hindcast or "reconstruct" past environmental conditions from sub-fossil species assemblages preserved in sediment cores. The techniques included are weighted averaging ``WA``, partial least squares (PLS) and weighted average partial least squared ``WAPLS``, Imbrie and Kipp Factor Analysis ``IKFA`` a form of principal components regression, Maximum Likelihood Response Curves ``MLRC``, and the Modern Analogue Technique ``MAT``, a form of k-NN non-parametric regression (see Juggins & Birks (2010) for a review). 

The techniques are implemented in a consistent way and include functions for fitting a model to a training set of species and environmental data, with the function named after the technique: that is, ``WA`` fits a weighted averaging model. Any model can be cross-validated using the ``crossval`` function, which allows internal cross-validation using leave-one-out, leave-n-out, bootstrapping or h-block cross-validation. There are a number of generic functions that can be used to summarise and diagnose the models: ``print``, ``summary``, ``performance`` and ``plot``. Some techniques have additional diagnostic functions such as ``screeplot`` and ``rand.t.test`` to help estimate the approproate number of components (WAPLS), factors (IKFA) or number of analogues (MAT).

