# sstvars 1.0.0

* Initial CRAN submission.

# sstvars 1.0.1

* Updated configure script to fix an issue with the installation on Mac OS X.

# sstvars 1.0.2

* Updated readme.
* Updated documentation.

# sstvars 1.1.0

* MAJOR: Implemented independent skewed t distribution as a new conditional distribution.
* MAJOR: Implemented a three phase estimation for TVAR models to enhance computational efficiency.
* MAJOR: Implemented a possibility to maximize penalized log-likelihood function that penalizes from unstable and close-to-unstable estimates.
  Significantly improves the performance of the estimation algorithm in some cases, particularly when the time series are very persistent. 
* Estimates not satisfying the usual stability condition for the regimes can now be allowed. 
* Adjusted the step sizes in finite difference numerical differentiation. 
* The step size in finite difference numerical differentiation can now be adjusted in the function iterate_more.
* Changed the random parameter generation for ind_Student models (estimation results with specific seeds are not backward compatible).
* A new function: filter_estimates, which can be used considers includes estimates that are not deemed inappropriate).
* A new function: plot_struct_shocks, which plots the structural shock time series. 
* A new function: stvar_to_sstvars110, which makes STVAR models estimated with package versions <1.1.0 compatible with package versions >=1.1.0.
* Some adjustments to estimation with fitSTVAR. NOTE: estimation results with a particular seed may be different to the earlier version. 
* Removed the argument "filter_estimates" from fitSTVAR as a redundancy (it is now always applied), since the function alt_stvar can in
  any case be used to browse the estimates from any estimation round. 
* Added a new functionality to fitSSTVAR: structural models identified by non-Gaussianity can be estimated based on different orderings
  or signs of the columns of any of B_1,...,B_M (to conveniently examine models corresponding to various orderings and signs in the presence
  of weak identification with respect to ordering or signs of the columns of B_2,...,B_M)
* FIXED A BUG in the simulation algorithm for models incorporating independent Student's t conditional distributions
  (the variance of each structural shock was not scaled to one). 
* FIXED A BUG in the GIRF simulation algorithm: the transition weights were not necessarily high for 'init_regime' at impact (but
  the initial values were generated from the correct regimes).
* Made the function profile_logliks more user friendly. 
* Added a simplified table of contents to the vignette. 
* The argument standard_error_print can now be used directly in the summary-function to obtain printout of standard errors. 
* Updated the documentation. 
