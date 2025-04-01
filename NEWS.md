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

# sstvars 1.1.1

* Now also the NLS step in the three-phase estimation estimation checks that there are enough observations from each regime
  (previously only LS estimation checked this).
* Added the argument min_obs_coef to fitSTVAR to let the user to control the smallest accepted number of observations from each regime in
  the LS/NLS step of the three-phase estimation. Also increased its default value.
* Now alt_stvar, iterate_more, and filter_estimates retain LS_estimates if the original model contains them. 
* Now summary printout of class sstvar objects tells if the log-likelihood function is penalized. 
* Fixed CRAN check issues. 

# sstvars 1.1.2

* A new feature in GFEVD: initval_type = "data" and use_data_shocks = TRUE now allows to filter the histories based on the dominance of
  a specific regime.
* A new feature in GIRF: use_data_shocks, which allows to estimate the GIRF using the length p histories in the data, using the shocks recovered
  from the fitted model, with the possibility to filter the histories based on the dominance of a specific regime as well as on the sign and
  size of the shocks.
* GFEVDs are now calculated as the average over the GFEVDs based on the different initial values (previously calculated based on the average of the GIRFs
  based on the different initial values).
* Fixed the function stvar_to_sstvars110.
* Fixed a bug that caused GFEVD with data shocks to result in error when using model estimate with package versions <1.1.0.
* Various small adjustments to printouts, documentation, etc. 

# sstvars 1.1.3

* Fixed a bug that prevented printing models that impose restrictions on the AR parameters. 
* Now fitSTVAR normalizes the first row of the impact matrix B_1 to be positive and in a decreasing order also for skewed t models. 
* Now swap_B_signs also swaps the signs of the appropriate skewness parameter values for skewed t models, so that the resulting model
  is observationally equivalent to the original model.
* Increased the default maxit from 1000 to 2000 in fitSTVAR.

# sstvars 1.1.4

* Now also for Gaussian models, the functions GIRF and GFEVD use simulation procedure to draw initial values from a specific regime when init_regime is used 
  and not use_data_shocks. This ensures that the transition weights are high for the initial regime at impact, while with stationary distribution that was not checked.
* In simulate.stvar the default for Gaussian models is now also use simulation procedure to draw initial values from a specific regime, ensuring that the transition weights
  are high for the initial regime based on the initial values. The newly argument use_stat_for_Gaus can be set to TRUE to use stationary distribution instead (as was the 
  old functionality).
* Added the possibility to use more sparse grid in LS/NLS estimation by adding the argument sparse_grid to fitSTVAR.
* Fixed a bug causing the argument min_obs_coef in fitSTVAR to be ignored in the LS/NLS estimation step.
* Fixed a bug in plot.girf in which the y-axis was not wide enough when the point estimate is outside the "confidence bounds".
* The object returned by fitSTVAR now includes also the seeds used in the estimation (if specified).
* Updated the datafile usacpu to include observations until the end of 2024.
* Updated the references.
* Updates to the documentation.

# sstvars 1.1.5

* Updated the vignette to align with latest version of the working paper \url{https://arxiv.org/abs/2404.19707} of Virolainen, 
  establishing identification by non-Gaussianity for TVAR and STVAR models.
* In the documentation, simplified the notation of the intercept parameters from \phi_{m,0} to \phi_{m}.

